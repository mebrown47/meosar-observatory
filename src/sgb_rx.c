// SPDX-License-Identifier: GPL-3.0-or-later
#include "sgb_rx.h"

#include <fftw3.h>
#include <math.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>

#include "fgb.h"   // MID table

#define TC (1.0 / SGB_CHIP_RATE)
#define NSEG SGB_PRE_BITS
#define MDFT 32
#define ACQ_DEC 2
#define TN 4096          // tile FFT length
#define KTILE 2          // lags kept per tile and bin for the coherent stage
#define BIN_STEP 10      // coarse frequency bin = 10 tile-FFT bins (152.6 Hz at 62.5 kS/s)

struct sgb_acq {
	double fs_in, fs, seg_T, span;
	int ntaps, Ls, T, threads, nbins, npre;
	float *h;
	int o[NSEG];             // segment start, samples from the preamble start
	fftwf_complex *S;        // NSEG x TN conjugate spectra of the short segment replicas
	fftwf_plan fwd, inv;
	int *sh;                 // bin shift in tile-FFT bins
	double *fb;              // bin frequency, Hz
	double complex *corr;    // nbins x NSEG: exp(-j 2pi fb o_j / fs), absolute phase reference
	float *xd;               // decimated block
	size_t nd, xcap;
};

sgb_acq_t *sgb_acq_new(double fs_in, double span_hz, int threads) {
	sgb_init();
	sgb_acq_t *a = calloc(1, sizeof *a);
	a->fs_in = fs_in;
	a->fs = fs_in / ACQ_DEC;
	a->threads = threads < 1 ? 1 : threads;
	a->seg_T = SGB_CHIPS_PER_BIT * TC;
	a->ntaps = 63;
	a->h = fir_lowpass(a->ntaps, 0.45 * a->fs, fs_in, 8.0);
	double S = a->seg_T * a->fs;                       // samples per segment (fractional)
	a->Ls = (int)ceil(S) + 2;
	a->T = TN - a->Ls + 1;
	a->npre = (int)ceil(NSEG * S) + 2;
	const float *ci = sgb_chips_i(SGB_NORMAL), *cq = sgb_chips_q(SGB_NORMAL);
	a->S = fftwf_malloc(sizeof(fftwf_complex) * (size_t)NSEG * TN);
	fftwf_complex *tmp = fftwf_malloc(sizeof(fftwf_complex) * TN), *tmp2 = fftwf_malloc(sizeof(fftwf_complex) * TN);
	a->fwd = fftwf_plan_dft_1d(TN, tmp, tmp2, FFTW_FORWARD, FFTW_MEASURE);
	a->inv = fftwf_plan_dft_1d(TN, tmp, tmp2, FFTW_BACKWARD, FFTW_MEASURE);
	for (int j = 0; j < NSEG; j++) {
		a->o[j] = (int)floor(j * S);
		float *g = (float *)tmp;
		memset(tmp, 0, sizeof(fftwf_complex) * TN);
		int lo = j * SGB_CHIPS_PER_BIT, hi = lo + SGB_CHIPS_PER_BIT;
		for (int m = 0; m < a->Ls; m++) {
			double t = (a->o[j] + m) / a->fs;
			int ki = (int)(t / TC), kq = (int)floor(t / TC - 0.5);
			if (ki >= lo && ki < hi) g[2 * m] += ci[ki];
			if (kq >= lo && kq < hi) g[2 * m + 1] += cq[kq];
		}
		fftwf_complex *Sj = a->S + (size_t)j * TN;
		fftwf_execute_dft(a->fwd, tmp, Sj);
		float *p = (float *)Sj;
		for (int k = 0; k < TN; k++) p[2 * k + 1] = -p[2 * k + 1];
	}
	fftwf_free(tmp); fftwf_free(tmp2);
	double fres = a->fs / TN;
	int B = (int)floor(span_hz / (BIN_STEP * fres));
	a->nbins = 2 * B + 1;
	a->sh = malloc(sizeof(int) * a->nbins);
	a->fb = malloc(sizeof(double) * a->nbins);
	a->corr = malloc(sizeof(double complex) * a->nbins * NSEG);
	for (int b = 0; b < a->nbins; b++) {
		a->sh[b] = (b - B) * BIN_STEP;
		a->fb[b] = a->sh[b] * fres;
		for (int j = 0; j < NSEG; j++) a->corr[b * NSEG + j] = cexp(-I * 2 * M_PI * a->fb[b] * a->o[j] / a->fs);
	}
	a->span = span_hz;
	return a;
}

double sgb_acq_span(const sgb_acq_t *a) { return 1.9; }

void sgb_acq_free(sgb_acq_t *a) {
	if (!a) return;
	fftwf_free(a->S);
	fftwf_destroy_plan(a->fwd); fftwf_destroy_plan(a->inv);
	free(a->h); free(a->sh); free(a->fb); free(a->corr); free(a->xd); free(a);
}

typedef struct { double p, sum; long cnt; int lag; double f; } binbest_t;

typedef struct {
	sgb_acq_t *a;
	int nlag, ntiles, next;
	pthread_mutex_t mu;
	binbest_t *res;          // merged per bin
} job_t;

static void *worker(void *arg) {
	job_t *J = arg;
	sgb_acq_t *a = J->a;
	int T = a->T, nb = a->nbins;
	fftwf_complex *X = fftwf_malloc(sizeof(fftwf_complex) * (size_t)NSEG * TN);
	fftwf_complex *in = fftwf_malloc(sizeof(fftwf_complex) * TN), *z = fftwf_malloc(sizeof(fftwf_complex) * TN);
	float *Zt = malloc(sizeof(float) * 2 * (size_t)NSEG * T), *Nc = malloc(sizeof(float) * T);
	binbest_t *loc = calloc((size_t)nb, sizeof(binbest_t));
	double complex tw[MDFT][NSEG];
	for (int m = 0; m < MDFT; m++)
		for (int j = 0; j < NSEG; j++) tw[m][j] = cexp(-I * 2 * M_PI * m * j / MDFT);
	for (;;) {
		pthread_mutex_lock(&J->mu);
		int ti = J->next++;
		pthread_mutex_unlock(&J->mu);
		if (ti >= J->ntiles) break;
		int L0 = ti * T, nt = J->nlag - L0 < T ? J->nlag - L0 : T;
		for (int j = 0; j < NSEG; j++) {
			float *ip = (float *)in;
			long s0 = (long)L0 + a->o[j];
			for (int k = 0; k < TN; k++) {
				long idx = s0 + k;
				if (idx < (long)a->nd) { ip[2 * k] = a->xd[2 * idx]; ip[2 * k + 1] = a->xd[2 * idx + 1]; }
				else ip[2 * k] = ip[2 * k + 1] = 0;
			}
			fftwf_execute_dft(a->fwd, in, X + (size_t)j * TN);
		}
		for (int b = 0; b < nb; b++) {
			int sh = a->sh[b];
			memset(Nc, 0, sizeof(float) * nt);
			for (int j = 0; j < NSEG; j++) {
				const float *Xj = (const float *)(X + (size_t)j * TN), *Sj = (const float *)(a->S + (size_t)j * TN);
				float *P = (float *)in;
				for (int k = 0; k < TN; k++) {
					int ks = (k + sh) & (TN - 1);
					float xr = Xj[2 * ks], xi = Xj[2 * ks + 1], gr = Sj[2 * k], gi = Sj[2 * k + 1];
					P[2 * k] = xr * gr - xi * gi;
					P[2 * k + 1] = xr * gi + xi * gr;
				}
				fftwf_execute_dft(a->inv, in, z);
				const float *zz = (const float *)z;
				float cr = (float)creal(a->corr[b * NSEG + j]), cim = (float)cimag(a->corr[b * NSEG + j]);
				float *Zj = Zt + 2 * (size_t)j * T;
				for (int l = 0; l < nt; l++) {
					float r = zz[2 * l], i = zz[2 * l + 1];
					Zj[2 * l] = r * cr - i * cim;
					Zj[2 * l + 1] = r * cim + i * cr;
					Nc[l] += r * r + i * i;
				}
			}
			// best KTILE lags of the non-coherent sum, at least 3 lags apart
			int pick[KTILE], np = 0;
			for (int q = 0; q < KTILE; q++) {
				int bi = -1;
				for (int l = 0; l < nt; l++) {
					int near = 0;
					for (int r = 0; r < np; r++) if (abs(l - pick[r]) < 3) near = 1;
					if (!near && (bi < 0 || Nc[l] > Nc[bi])) bi = l;
				}
				if (bi >= 0) pick[np++] = bi;
			}
			double sum = 0;
			for (int l = 0; l < nt; l++) sum += Nc[l];
			loc[b].sum += sum;
			loc[b].cnt += nt;
			for (int q = 0; q < np; q++) {
				int l = pick[q];
				for (int m = 0; m < MDFT; m++) {
					double complex y = 0;
					for (int j = 0; j < NSEG; j++) {
						const float *Zj = Zt + 2 * (size_t)j * T;
						y += (Zj[2 * l] + I * Zj[2 * l + 1]) * tw[m][j];
					}
					double p = creal(y) * creal(y) + cimag(y) * cimag(y);
					if (p > loc[b].p) {
						int mf = m < MDFT / 2 ? m : m - MDFT;
						loc[b].p = p; loc[b].lag = L0 + l; loc[b].f = a->fb[b] + mf / (MDFT * a->seg_T);
					}
				}
			}
		}
	}
	pthread_mutex_lock(&J->mu);
	for (int b = 0; b < nb; b++) {
		J->res[b].sum += loc[b].sum;
		J->res[b].cnt += loc[b].cnt;
		if (loc[b].p > J->res[b].p) { J->res[b].p = loc[b].p; J->res[b].lag = loc[b].lag; J->res[b].f = loc[b].f; }
	}
	pthread_mutex_unlock(&J->mu);
	fftwf_free(X); fftwf_free(in); fftwf_free(z); free(Zt); free(Nc); free(loc);
	return NULL;
}

static int by_z(const void *x, const void *y) {
	double d = ((const sgb_peak_t *)y)->z - ((const sgb_peak_t *)x)->z;
	return (d > 0) - (d < 0);
}

int sgb_acq_search(sgb_acq_t *a, const cf_t *x, size_t n, double t_max, sgb_peak_t *out, int max) {
	size_t nd = n / ACQ_DEC;
	if (nd > a->xcap) { a->xcap = nd; a->xd = realloc(a->xd, sizeof(float) * 2 * nd); }
	a->nd = nd;
	int c = (a->ntaps - 1) / 2;
	for (size_t k = 0; k < nd; k++) {
		float re = 0, im = 0;
		for (int j = 0; j < a->ntaps; j++) {
			long idx = (long)k * ACQ_DEC + c - j;
			if (idx < 0 || idx >= (long)n) continue;
			re += a->h[j] * crealf(x[idx]);
			im += a->h[j] * cimagf(x[idx]);
		}
		a->xd[2 * k] = re; a->xd[2 * k + 1] = im;
	}
	int nlag = (int)fmin(t_max * a->fs, (double)(nd > (size_t)a->npre ? nd - a->npre : 0));
	if (nlag <= 0) return 0;
	job_t J = {.a = a, .nlag = nlag, .ntiles = (nlag + a->T - 1) / a->T, .next = 0};
	pthread_mutex_init(&J.mu, NULL);
	J.res = calloc((size_t)a->nbins, sizeof(binbest_t));
	pthread_t th[32];
	int nt = a->threads > 32 ? 32 : a->threads;
	for (int i = 0; i < nt; i++) pthread_create(&th[i], NULL, worker, &J);
	for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
	sgb_peak_t *pk = malloc(sizeof(sgb_peak_t) * a->nbins);
	for (int b = 0; b < a->nbins; b++) {
		double mean = J.res[b].cnt ? J.res[b].sum / J.res[b].cnt : 1;
		pk[b].z = J.res[b].p / mean;
		pk[b].t = J.res[b].lag / a->fs;
		pk[b].f = J.res[b].f;
	}
	qsort(pk, (size_t)a->nbins, sizeof(sgb_peak_t), by_z);
	int k = 0;
	for (int i = 0; i < a->nbins && k < max; i++) out[k++] = pk[i];
	free(pk); free(J.res);
	pthread_mutex_destroy(&J.mu);
	return k;
}

// ---------------- demodulation ----------------

// chip integrals of x*exp(-j2pi f t) over chips starting at t0 + off*Tc, nchips of them
static void chip_int(const cf_t *x, size_t n, double fs, double f, double t0, int nchips, double off,
                     double complex *out) {
	// mix and integrate on the fly over the needed sample span
	double a0 = (t0 + off * TC) * fs, step = TC * fs;
	long s0 = (long)floor(a0) - 1, s1 = (long)ceil(a0 + nchips * step) + 1;
	if (s0 < 0) s0 = 0;
	if (s1 > (long)n) s1 = (long)n;
	long m = s1 - s0;
	double complex *cs = malloc(sizeof(double complex) * (size_t)(m + 1));
	cs[0] = 0;
	double complex rot = cexp(-I * 2 * M_PI * f * s0 / fs), stp = cexp(-I * 2 * M_PI * f / fs);
	for (long i = 0; i < m; i++) {
		cs[i + 1] = cs[i] + x[s0 + i] * rot;
		rot *= stp;
	}
	double complex prev = 0;
	for (int k = 0; k <= nchips; k++) {
		double e = a0 + k * step - s0;
		double complex v;
		if (e <= 0) v = 0;
		else if (e >= m) v = cs[m];
		else { long i = (long)e; double fr = e - i; v = cs[i] * (1 - fr) + cs[i + 1] * fr; }
		if (k) out[k - 1] = v - prev;
		prev = v;
	}
	free(cs);
}

static double complex pre_sum_bits(const cf_t *x, size_t n, double fs, double f, double t0, const float *ci,
                                   const float *cq, double complex *zb) {
	int npc = SGB_PRE_BITS * SGB_CHIPS_PER_BIT;
	double complex *a = malloc(sizeof(double complex) * npc), *b = malloc(sizeof(double complex) * npc);
	chip_int(x, n, fs, f, t0, npc, 0.0, a);
	chip_int(x, n, fs, f, t0, npc, 0.5, b);
	double complex tot = 0;
	for (int k = 0; k < SGB_PRE_BITS; k++) {
		double complex s = 0;
		for (int c = 0; c < SGB_CHIPS_PER_BIT; c++) {
			int i = k * SGB_CHIPS_PER_BIT + c;
			s += a[i] * ci[i] - I * b[i] * cq[i];
		}
		zb[k] = s;
		tot += s;
	}
	free(a); free(b);
	return tot;
}

static double fft64_peak(const double complex *z) {   // max |DFT_64| of 25 values
	double best = 0;
	for (int m = 0; m < 64; m++) {
		double complex s = 0;
		for (int k = 0; k < SGB_PRE_BITS; k++) s += z[k] * cexp(-I * 2 * M_PI * m * k / 64);
		if (cabs(s) > best) best = cabs(s);
	}
	return best;
}

static void movavg_nearest(const double *x, int n, int w, double *y) {
	int h = w / 2;
	for (int i = 0; i < n; i++) {
		double s = 0;
		for (int j = i - h; j <= i - h + w - 1; j++) s += x[j < 0 ? 0 : j >= n ? n - 1 : j];
		y[i] = s / w;
	}
}

void sgb_demodulate(const cf_t *x, size_t n, double fs, double t0g, double fg, int mode, sgb_demod_t *d) {
	memset(d, 0, sizeof *d);
	const float *ci = sgb_chips_i(mode), *cq = sgb_chips_q(mode);
	double complex zb[SGB_PRE_BITS];
	double best = -1, t0 = t0g, f = fg;
	for (double dtc = -1.5; dtc <= 1.501; dtc += 0.25)
		for (double dfh = -20; dfh <= 20.01; dfh += 4) {
			pre_sum_bits(x, n, fs, fg + dfh, t0g + dtc * TC, ci, cq, zb);
			double p = fft64_peak(zb);
			if (p > best) { best = p; t0 = t0g + dtc * TC; f = fg + dfh; }
		}
	pre_sum_bits(x, n, fs, f, t0, ci, cq, zb);
	// residual frequency from the phase slope over the preamble bits
	double ph[SGB_PRE_BITS], sx = 0, sy = 0, sxx = 0, sxy = 0;
	for (int k = 0; k < SGB_PRE_BITS; k++) {
		ph[k] = carg(zb[k]);
		if (k) { while (ph[k] - ph[k - 1] > M_PI) ph[k] -= 2 * M_PI; while (ph[k] - ph[k - 1] < -M_PI) ph[k] += 2 * M_PI; }
		sx += k; sy += ph[k]; sxx += (double)k * k; sxy += k * ph[k];
	}
	double slope = (SGB_PRE_BITS * sxy - sx * sy) / (SGB_PRE_BITS * sxx - sx * sx);
	f += slope / (2 * M_PI * SGB_CHIPS_PER_BIT * TC);
	double complex tot = pre_sum_bits(x, n, fs, f, t0, ci, cq, zb);
	for (double dtc = -0.15; dtc <= 0.151; dtc += 0.05) {
		double complex zz[SGB_PRE_BITS];
		double complex tt = pre_sum_bits(x, n, fs, f, t0 + dtc * TC, ci, cq, zz);
		if (cabs(tt) > cabs(tot)) { tot = tt; t0 += dtc * TC; memcpy(zb, zz, sizeof zz); }
	}
	d->t0 = t0;
	// whole burst
	double complex *a = malloc(sizeof(double complex) * SGB_NCHIPS), *b = malloc(sizeof(double complex) * SGB_NCHIPS);
	chip_int(x, n, fs, f, t0, SGB_NCHIPS, 0.0, a);
	chip_int(x, n, fs, f, t0, SGB_NCHIPS, 0.5, b);
	double complex zi[SGB_BITS_PER_CH], zq[SGB_BITS_PER_CH];
	for (int k = 0; k < SGB_BITS_PER_CH; k++) {
		double complex si = 0, sq = 0;
		for (int c = 0; c < SGB_CHIPS_PER_BIT; c++) {
			int i = k * SGB_CHIPS_PER_BIT + c;
			si += a[i] * ci[i];
			sq += -I * b[i] * cq[i];
		}
		zi[k] = si; zq[k] = sq;
	}
	free(a); free(b);
	// frequency residual and drift from the squared bit correlations
	double tb[SGB_BITS_PER_CH], tm = 0;
	for (int k = 0; k < SGB_BITS_PER_CH; k++) { tb[k] = (k + 0.5) * SGB_CHIPS_PER_BIT * TC; tm += tb[k]; }
	tm /= SGB_BITS_PER_CH;
	for (int k = 0; k < SGB_BITS_PER_CH; k++) tb[k] -= tm;
	double complex w[SGB_BITS_PER_CH];
	for (int k = 0; k < SGB_BITS_PER_CH; k++) w[k] = zi[k] * zi[k] + zq[k] * zq[k];
	double bs = -1, bdf = 0, bfd = 0;
	for (double df = -3.0; df <= 3.001; df += 0.1)
		for (double fd = -10.0; fd <= 10.001; fd += 0.25) {
			double complex s = 0;
			for (int k = 0; k < SGB_BITS_PER_CH; k++) {
				double th = 2 * M_PI * (df * tb[k] + 0.5 * fd * tb[k] * tb[k]);
				s += w[k] * cexp(-2 * I * th);
			}
			if (cabs(s) > bs) { bs = cabs(s); bdf = df; bfd = fd; }
		}
	double theta[SGB_BITS_PER_CH];
	for (int k = 0; k < SGB_BITS_PER_CH; k++) theta[k] = 2 * M_PI * (bdf * tb[k] + 0.5 * bfd * tb[k] * tb[k]);
	double complex p0 = 0;
	for (int k = 0; k < SGB_PRE_BITS; k++) p0 += (zi[k] + zq[k]) * cexp(-I * theta[k]);
	double phi0 = carg(p0), phase[SGB_BITS_PER_CH];
	for (int k = 0; k < SGB_BITS_PER_CH; k++) phase[k] = phi0 + theta[k];
	d->f_hz = f + bdf;
	d->fdot = bfd;
	double re[SGB_BITS_PER_CH], im[SGB_BITS_PER_CH], sre[SGB_BITS_PER_CH], sim[SGB_BITS_PER_CH];
	for (int it = 0; it < 3; it++) {
		for (int k = 0; k < SGB_BITS_PER_CH; k++) {
			double complex rot = cexp(-I * phase[k]);
			double si = k < SGB_PRE_BITS ? 1 : (creal(zi[k] * rot) >= 0 ? 1 : -1);
			double sq = k < SGB_PRE_BITS ? 1 : (creal(zq[k] * rot) >= 0 ? 1 : -1);
			double complex res = (zi[k] * si + zq[k] * sq) * cexp(-I * (phi0 + theta[k]));
			re[k] = creal(res); im[k] = cimag(res);
		}
		movavg_nearest(re, SGB_BITS_PER_CH, 31, sre);
		movavg_nearest(im, SGB_BITS_PER_CH, 31, sim);
		for (int k = 0; k < SGB_BITS_PER_CH; k++) phase[k] = phi0 + theta[k] + atan2(sim[k], sre[k]);
	}
	double sig2 = 0, mu = 0;
	double ri[SGB_BITS_PER_CH], rq[SGB_BITS_PER_CH];
	for (int k = 0; k < SGB_BITS_PER_CH; k++) {
		double complex rot = cexp(-I * phase[k]);
		double complex a1 = zi[k] * rot, b1 = zq[k] * rot;
		ri[k] = creal(a1); rq[k] = creal(b1);
		if (k >= SGB_PRE_BITS) {
			sig2 += cimag(a1) * cimag(a1) + cimag(b1) * cimag(b1);
			mu += fabs(ri[k]) + fabs(rq[k]);
		}
	}
	int nm = 2 * (SGB_BITS_PER_CH - SGB_PRE_BITS);
	sig2 /= nm;
	mu /= nm;
	for (int k = 0; k < SGB_MSG_BITS / 2; k++) {
		d->llr[2 * k] = -2 * mu * ri[SGB_PRE_BITS + k] / sig2;
		d->llr[2 * k + 1] = -2 * mu * rq[SGB_PRE_BITS + k] / sig2;
	}
	d->ebn0_db = 10 * log10(fmax(mu * mu / (2 * sig2), 1e-6));
	d->pre_z = cabs(tot) / sqrt(sig2 * SGB_PRE_BITS * 2);
	d->ok = 1;
}

static int frame_checks(const uint8_t *bits, char *why, size_t wl) {
	int mid = 0, any = 0;
	for (int i = 30; i < 40; i++) mid = mid << 1 | bits[i];
	for (int i = 0; i < 202; i++) any |= bits[i];
	why[0] = 0;
	if (!any) { snprintf(why, wl, "all-zero frame"); return 0; }
	if (!fgb_mid_valid(mid)) { snprintf(why, wl, "country %d not assigned", mid); return 0; }
	return 1;
}

int sgb_process_block(sgb_acq_t *a, const cf_t *y, size_t n, double f1, double mid, double t_y0,
                      double t_lo, double t_hi, double acq_thr, double prez_thr, sgb_rec_t *out, int max) {
	int nout = 0;
	double span = sgb_acq_span(a);
	for (double ts = t_lo; ts < t_hi; ts += span) {
		double tmax = fmin(span, t_hi - ts);
		long i0 = lround((ts - t_y0) * f1);
		if (i0 < 0) continue;
		size_t need = (size_t)((tmax + 1.2) * f1);
		if ((size_t)i0 + need > n) need = n - (size_t)i0;
		sgb_peak_t pk[12];
		int np = sgb_acq_search(a, y + i0, need, tmax, pk, 12);
		sgb_demod_t dm[12];
		int nk = 0;
		for (int p = 0; p < np && nk < 12; p++) {
			if (pk[p].z < acq_thr) continue;
			sgb_demodulate(y + i0, need, f1, pk[p].t, pk[p].f, SGB_NORMAL, &dm[nk]);
			if (!dm[nk].ok || dm[nk].pre_z < prez_thr) continue;
			// same relay found twice (neighbouring bins)
			int dup = 0;
			for (int q = 0; q < nk; q++)
				if (fabs(dm[q].t0 - dm[nk].t0) < 0.002 && fabs(dm[q].f_hz - dm[nk].f_hz) < 20) dup = 1;
			if (dup) continue;
			nk++;
		}
		char used[12] = {0};
		for (int i = 0; i < nk && nout < max; i++) {
			if (used[i]) continue;
			int grp[8], ng = 0;
			grp[ng++] = i;
			used[i] = 1;
			for (int j = i + 1; j < nk && ng < 8; j++)
				if (!used[j] && fabs(dm[j].t0 - dm[i].t0) < 0.040) { grp[ng++] = j; used[j] = 1; }
			sgb_rec_t *r = &out[nout];
			memset(r, 0, sizeof *r);
			double llr[SGB_MSG_BITS] = {0};
			for (int k = 0; k < ng; k++) {
				sgb_demod_t *q = &dm[grp[k]];
				for (int b = 0; b < SGB_MSG_BITS; b++) llr[b] += q->llr[b];
				uint8_t hb[SGB_MSG_BITS];
				for (int b = 0; b < SGB_MSG_BITS; b++) hb[b] = q->llr[b] > 0;
				r->ralone[k] = sgb_bm_decode(hb) >= 0;
				r->rf[k] = mid + q->f_hz; r->rdt_ms[k] = 1e3 * (q->t0 - dm[i].t0);
				r->reb[k] = q->ebn0_db; r->rpz[k] = q->pre_z;
			}
			r->nrelay = ng;
			r->t = ts + dm[i].t0;
			r->f_hz = mid + dm[i].f_hz;
			if (r->t < t_lo || r->t >= t_hi) continue;   // owned by another block
			for (int b = 0; b < SGB_MSG_BITS; b++) r->bits[b] = llr[b] > 0;
			r->nerr = sgb_bm_decode(r->bits);
			if (r->nerr >= 0) { r->ok = 1; r->hard = 1; }
			else {
				double cost;
				r->flips = sgb_chase(llr, 8, r->bits, &cost);
				r->ok = r->flips >= 0 && r->flips <= 8;
			}
			if (r->ok && !frame_checks(r->bits, r->reasons, sizeof r->reasons)) r->ok = 0;
			r->accept = !r->ok ? 0 : r->hard ? 2 : 1;
			nout++;
		}
	}
	return nout;
}
