#!/usr/bin/env python3
"""
skycov.py - where in the sky this station hears the SAR satellites (the Coverage page).

For every burst that satellite attribution (satattr.py) placed clearly, this asks which SAR
satellites could have relayed it: above the station's horizon (satattr.MIN_EL, as attribution
only places copies above it), and high enough above the
beacon's (BEACON_MIN_EL) for the beacon's uplink not to be the weak link. Each of them is an
expected copy, counted in the sky cell (10 degrees of azimuth by 5 of elevation, as the station
sees it) where that satellite was; the ones the burst was actually relayed by are heard copies,
with their Eb/N0. A cell's heard / expected is the station's hear rate in that part of the sky:
trees, roofs, the mounting and the antenna's pattern show as a horizon that rises.

A SAR satellite is one that has relayed bursts here (copies.sat): CelesTrak's GLONASS group
holds GLONASS-M satellites too, which carry no SAR transponder. A stray copy attributed to a
satellite that relays nothing (one copy in three weeks, against hundreds for the others) doesn't
make it one: it needs SAR_SHARE of the typical satellite's copies, and at least SAR_MIN.

Counts are hourly rows in cov_hourly, kept with the observatory database. Called by rollup.py
after satattr.run; the server's /api/coverage reads them.

Copyright (c) 2026 Mike Brown
SPDX-License-Identifier: GPL-3.0-or-later
"""

import satattr as SA

AZ_STEP, EL_STEP = 10, 5
BEACON_MIN_EL = 10.0    # deg: below it a miss may be the beacon's uplink, not the station
PER_RUN = 4000          # bursts per rollup run: a backlog is worked off over a few minutes
SAR_SHARE, SAR_MIN = 0.01, 3
VERSION = "2"           # change with the counting: every burst is counted again


def sar_sats(con):
    """The satellites that relay bursts here."""
    n = dict(con.execute("SELECT sat, COUNT(*) FROM copies WHERE sat IS NOT NULL GROUP BY sat"))
    if not n:
        return set()
    typical = sorted(n.values())[len(n) // 2]
    return {s for s, k in n.items() if k >= max(SAR_MIN, SAR_SHARE * typical)}


def init_tables(con):
    con.executescript(
        """
        -- Per hour, satellite and sky cell (az = azimuth // 10 deg, el = elevation // 5 deg, as
        -- the station sees the satellite): copies expected, copies heard, and the heard copies'
        -- Eb/N0 sum and count.
        CREATE TABLE IF NOT EXISTS cov_hourly (
            hour INTEGER, sat TEXT, az INTEGER, el INTEGER, exp INTEGER, heard INTEGER,
            eb_sum REAL, eb_n INTEGER, PRIMARY KEY (hour, sat, az, el));
        """
    )
    # attr.cov: the burst has been counted (NULL: not yet). satattr rewrites attr rows when it
    # attributes again, which leaves them uncounted.
    if "cov" not in [c[1] for c in con.execute("PRAGMA table_info(attr)")]:
        con.execute("ALTER TABLE attr ADD COLUMN cov INTEGER")


def run(con, where, tle_dir, verbose=False):
    """Count the clearly attributed bursts not counted yet. Returns how many, or None when off."""
    if SA.Satrec is None or not where:
        return None
    if con.execute("SELECT 1 FROM sqlite_master WHERE name='attr'").fetchone() is None:
        return None
    init_tables(con)
    # A new attribution model attributes every burst again: count them again from scratch.
    attr_model = (con.execute("SELECT value FROM meta WHERE key='attr_model'").fetchone() or [None])[0]
    model = f"{attr_model}/{VERSION}"
    if (con.execute("SELECT value FROM meta WHERE key='cov_model'").fetchone() or [None])[0] != model:
        con.execute("DELETE FROM cov_hourly")
        con.execute("UPDATE attr SET cov=NULL")
        con.execute("INSERT INTO meta(key,value) VALUES('cov_model',?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
                    (model,))
    tles = SA.Tles(tle_dir)
    if not tles.snaps:
        return None
    sar = sar_sats(con)
    station = SA.ecef(where[0], where[1], 0.1)
    rows = con.execute("SELECT gen, burst, t, sats FROM attr WHERE state='clear' AND cov IS NULL ORDER BY t LIMIT ?",
                       (PER_RUN,)).fetchall()
    acc = {}
    done = []
    for gen, burst, t, sats in rows:
        done.append((gen, burst))
        table = "fgb_bursts" if gen == 1 else "sgb_bursts"
        pos = con.execute(f"SELECT f.lat, f.lon FROM {table} b JOIN frames f ON f.hex=b.hex WHERE b.id=?", (burst,)).fetchone()
        if not pos or pos[0] is None:
            continue
        beacon = SA.ecef(pos[0], pos[1])
        eb = {s: e for s, e in con.execute("SELECT sat, ebn0 FROM copies WHERE gen=? AND burst=? AND sat IS NOT NULL",
                                          (gen, burst))}
        heard = set(sats.split(",")) if sats else set()
        hour = int(t) // 3600 * 3600
        for band in SA.GROUPS:
            for label, _, sat in tles.at(t, band):
                if label not in sar:
                    continue
                st = SA.sat_state(sat, t)
                if st is None:
                    continue
                r, v = st
                _, _, el = SA.look(r, v, station)
                if el < SA.MIN_EL:     # attribution never places a copy below it
                    continue
                if SA.look(r, v, beacon)[2] < BEACON_MIN_EL:
                    continue
                k = (hour, label, int(SA.azimuth(r, station) // AZ_STEP), min(int(el // EL_STEP), 90 // EL_STEP - 1))
                a = acc.setdefault(k, [0, 0, 0.0, 0])
                a[0] += 1
                if label in heard:
                    a[1] += 1
                    if eb.get(label) is not None:
                        a[2] += eb[label]
                        a[3] += 1
    for (hour, sat, az, el), (n, h, es, en) in acc.items():
        con.execute("INSERT INTO cov_hourly VALUES(?,?,?,?,?,?,?,?) ON CONFLICT(hour, sat, az, el) DO UPDATE SET "
                    "exp=exp+excluded.exp, heard=heard+excluded.heard, eb_sum=eb_sum+excluded.eb_sum, eb_n=eb_n+excluded.eb_n",
                    (hour, sat, az, el, n, h, es, en))
    con.executemany("UPDATE attr SET cov=1 WHERE gen=? AND burst=?", done)
    con.commit()
    if verbose and done:
        print(f"skycov: {len(done)} bursts counted")
    return len(done)
