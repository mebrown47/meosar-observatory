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
#define KLAGS 32
#define ACQ_DEC 2
#define ACQ_NFFT (1 << 17)

struct sgb_acq {
	double fs_in, fs, span, seg_T;
	int nfft, npre, nbins, threads;
	float *h;                 // decimation filter
	int ntaps;
	fftwf_complex *G;         // NSEG x nfft conj replica spectra
	fftwf_plan fwd, inv;
	fftwf_complex *X;
	double *bins;
};

sgb_acq_t *sgb_acq_new(double fs_in, double span_hz, int threads) {
	sgb_init();
	sgb_acq_t *a = calloc(1, sizeof *a);
	a->fs_in = fs_in;
	a->fs = fs_in / ACQ_DEC;
	a->nfft = ACQ_NFFT;
	a->threads = threads < 1 ? 1 : threads;
	a->seg_T = SGB_CHIPS_PER_BIT * TC;
	a->ntaps = 63;
	a->h = fir_lowpass(a->ntaps, 0.45 * a->fs, fs_in, 8.0);
	a->npre = (int)ceil(NSEG * a->seg_T * a->fs) + 2;
	const float *ci = sgb_chips_i(SGB_NORMAL), *cq = sgb_chips_q(SGB_NORMAL);
	a->G = fftwf_malloc(sizeof(fftwf_complex) * (size_t)NSEG * a->nfft);
	a->X = fftwf_malloc(sizeof(fftwf_complex) * a->nfft);
	fftwf_complex *tmp = fftwf_malloc(sizeof(fftwf_complex) * a->nfft);
	a->fwd = fftwf_plan_dft_1d(a->nfft, tmp, a->X, FFTW_FORWARD, FFTW_MEASURE);
	a->inv = fftwf_plan_dft_1d(a->nfft, a->X, tmp, FFTW_BACKWARD, FFTW_MEASURE);
	float *g = (float *)tmp;
	for (int j = 0; j < NSEG; j++) {
		memset(tmp, 0, sizeof(fftwf_complex) * a->nfft);
		for (int i = 0; i < a->npre; i++) {
			double t = i / a->fs;
			int ki = (int)(t / TC), kq = (int)floor(t / TC - 0.5);
			int lo = j * SGB_CHIPS_PER_BIT, hi = lo + SGB_CHIPS_PER_BIT;
			if (ki >= lo && ki < hi) g[2 * i] += ci[ki];
			if (kq >= lo && kq < hi) g[2 * i + 1] += cq[kq];
		}
		fftwf_complex *Gj = a->G + (size_t)j * a->nfft;
		fftwf_execute_dft(a->fwd, tmp, Gj);
		float *p = (float *)Gj;
		for (int k = 0; k < a->nfft; k++) p[2 * k + 1] = -p[2 * k + 1];   // conjugate
	}
	fftwf_free(tmp);
	double df = 1.0 / a->seg_T;
	a->nbins = 2 * (int)floor(span_hz / df) + 1;
	a->bins = malloc(sizeof(double) * a->nbins);
	for (int b = 0; b < a->nbins; b++) a->bins[b] = (b - a->nbins / 2) * df;
	a->span = span_hz;
	return a;
}

double sgb_acq_span(const sgb_acq_t *a) { return (a->nfft - a->npre) / a->fs; }

void sgb_acq_free(sgb_acq_t *a) {
	if (!a) return;
	fftwf_free(a->G); fftwf_free(a->X);
	fftwf_destroy_plan(a->fwd); fftwf_destroy_plan(a->inv);
	free(a->h); free(a->bins); free(a);
}

typedef struct {
	sgb_acq_t *a;
	int nlag;
	int next;                    // shared bin counter
	pthread_mutex_t mu;
	sgb_peak_t *best;            // one per bin
} job_t;

static int cmp_desc(const void *x, const void *y) {
	float a = *(const float *)x, b = *(const float *)y;
	return (a < b) - (a > b);
}

static void *worker(void *arg) {
	job_t *J = arg;
	sgb_acq_t *a = J->a;
	int N = a->nfft, nlag = J->nlag;
	fftwf_complex *prod = fftwf_malloc(sizeof(fftwf_complex) * N);
	fftwf_complex *zt = fftwf_malloc(sizeof(fftwf_complex) * N);
	float *Z = malloc(sizeof(float) * 2 * (size_t)NSEG * nlag);
	float *Nc = malloc(sizeof(float) * nlag), *srt = malloc(sizeof(float) * nlag);
	double fres = a->fs / N;
	for (;;) {
		pthread_mutex_lock(&J->mu);
		int b = J->next++;
		pthread_mutex_unlock(&J->mu);
		if (b >= a->nbins) break;
		int sh = (int)lround(a->bins[b] / fres);
		memset(Nc, 0, sizeof(float) * nlag);
		const float *X = (const float *)a->X;
		for (int j = 0; j < NSEG; j++) {
			const float *G = (const float *)(a->G + (size_t)j * N);
			float *P = (float *)prod;
			for (int k = 0; k < N; k++) {
				int ks = (k + sh) & (N - 1);           // roll(X, -sh)
				float xr = X[2 * ks], xi = X[2 * ks + 1], gr = G[2 * k], gi = G[2 * k + 1];
				P[2 * k] = xr * gr - xi * gi;
				P[2 * k + 1] = xr * gi + xi * gr;
			}
			fftwf_execute_dft(a->inv, prod, zt);
			const float *zz = (const float *)zt;
			float *Zj = Z + 2 * (size_t)j * nlag;
			for (int l = 0; l < nlag; l++) {
				Zj[2 * l] = zz[2 * l]; Zj[2 * l + 1] = zz[2 * l + 1];
				Nc[l] += zz[2 * l] * zz[2 * l] + zz[2 * l + 1] * zz[2 * l + 1];
			}
		}
		double mean = 0;
		for (int l = 0; l < nlag; l++) mean += Nc[l];
		mean /= nlag;
		// K-th largest non-coherent value as the screen
		memcpy(srt, Nc, sizeof(float) * nlag);
		int K = KLAGS < nlag ? KLAGS : nlag;
		// partial selection: nth element via qsort of a sample is simpler but slow; use a heap-free pass
		float thr;
		{
			// quickselect for the K-th largest
			int lo = 0, hi = nlag - 1, k = K - 1;
			while (lo < hi) {
				float p = srt[(lo + hi) / 2];
				int i = lo, jx = hi;
				while (i <= jx) {
					while (srt[i] > p) i++;
					while (srt[jx] < p) jx--;
					if (i <= jx) { float t = srt[i]; srt[i] = srt[jx]; srt[jx] = t; i++; jx--; }
				}
				if (k <= jx) hi = jx; else if (k >= i) lo = i; else break;
			}
			thr = srt[k];
		}
		sgb_peak_t pk = {0, 0, 0};
		for (int l = 0; l < nlag; l++) {
			if (Nc[l] < thr) continue;
			for (int m = 0; m < MDFT; m++) {
				double re = 0, im = 0;
				for (int j = 0; j < NSEG; j++) {
					double ph = -2 * M_PI * m * j / MDFT, c = cos(ph), s = sin(ph);
					const float *Zj = Z + 2 * (size_t)j * nlag;
					re += Zj[2 * l] * c - Zj[2 * l + 1] * s;
					im += Zj[2 * l] * s + Zj[2 * l + 1] * c;
				}
				double p = (re * re + im * im) / mean;
				if (p > pk.z) {
					int mf = m < MDFT / 2 ? m : m - MDFT;
					pk.z = p; pk.t = l / a->fs; pk.f = a->bins[b] + mf / (MDFT * a->seg_T);
				}
			}
		}
		J->best[b] = pk;
	}
	fftwf_free(prod); fftwf_free(zt); free(Z); free(Nc); free(srt);
	return NULL;
}

static int by_z(const void *x, const void *y) {
	double d = ((const sgb_peak_t *)y)->z - ((const sgb_peak_t *)x)->z;
	return (d > 0) - (d < 0);
}

int sgb_acq_search(sgb_acq_t *a, const cf_t *x, size_t n, double t_max, sgb_peak_t *out, int max) {
	(void)cmp_desc;
	// decimate by 2 with the acquisition filter (output k <-> input 2k)
	size_t nd = n / ACQ_DEC;
	if (nd > (size_t)a->nfft) nd = (size_t)a->nfft;
	float *X = (float *)a->X;
	fftwf_complex *tmp = fftwf_malloc(sizeof(fftwf_complex) * a->nfft);
	float *t = (float *)tmp;
	memset(tmp, 0, sizeof(fftwf_complex) * a->nfft);
	int c = (a->ntaps - 1) / 2;
	for (size_t k = 0; k < nd; k++) {
		float re = 0, im = 0;
		for (int j = 0; j < a->ntaps; j++) {
			long idx = (long)k * ACQ_DEC + c - j;
			if (idx < 0 || idx >= (long)n) continue;
			re += a->h[j] * crealf(x[idx]);
			im += a->h[j] * cimagf(x[idx]);
		}
		t[2 * k] = re; t[2 * k + 1] = im;
	}
	fftwf_execute_dft(a->fwd, tmp, a->X);
	fftwf_free(tmp);
	(void)X;
	int nlag = (int)fmin(t_max * a->fs, (double)(nd > (size_t)a->npre ? nd - a->npre : 0));
	if (nlag <= 0) return 0;
	job_t J = {.a = a, .nlag = nlag, .next = 0};
	pthread_mutex_init(&J.mu, NULL);
	J.best = calloc((size_t)a->nbins, sizeof(sgb_peak_t));
	pthread_t th[32];
	int nt = a->threads > 32 ? 32 : a->threads;
	for (int i = 0; i < nt; i++) pthread_create(&th[i], NULL, worker, &J);
	for (int i = 0; i < nt; i++) pthread_join(th[i], NULL);
	qsort(J.best, (size_t)a->nbins, sizeof(sgb_peak_t), by_z);
	int k = 0;
	for (int i = 0; i < a->nbins && k < max; i++) out[k++] = J.best[i];
	free(J.best);
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
