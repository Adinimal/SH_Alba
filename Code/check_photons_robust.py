#!/usr/bin/env python
"""
Robust high-energy photon checker for Fermi-LAT data.
Reports per-file properties (energy range, time range, size),
PSF type (handling Pass 8 32X bit-column format), angular separations,
and MJD for each photon above a given energy threshold.
"""
import astropy.io.fits as fits
import numpy as np
import os, sys

# ---- Configuration ----
RA_SRC   = 80.737
DEC_SRC  = -36.4686
EMIN_MEV = 300000.0     # 300 GeV in MeV
SEP_DEG  = 0.5
EVTLIST  = '/sdf/data/fermi/u/adinesh/DATA/J0522/events_j0522.txt'
FERMI_MJD_REF = 51910.0

PSF_NAMES = {4: 'PSF0', 8: 'PSF1', 16: 'PSF2', 32: 'PSF3'}

def angsep_deg(ra1, dec1, ra2, dec2):
    """Great-circle angular separation (degrees), vectorized."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    cosc = np.clip(np.sin(dec1)*np.sin(dec2) +
                   np.cos(dec1)*np.cos(dec2)*np.cos(ra1-ra2), -1.0, 1.0)
    return np.degrees(np.arccos(cosc))

def get_psf_type(event_type_val):
    """
    Decode PSF type from EVENT_TYPE.
    Pass 8 stores EVENT_TYPE as a 32X bit-column: astropy reads it
    as a boolean array of 32 bits (LSB first). Convert to int first.
    Bit 2 (value 4)=PSF0, bit 3 (8)=PSF1, bit 4 (16)=PSF2, bit 5 (32)=PSF3.
    """
    try:
        etv = np.asarray(event_type_val).ravel()
        if len(etv) == 32:
            # 32X boolean bit array: reconstruct integer (LSB first)
            val = int(np.sum([int(b) << i for i, b in enumerate(etv)]))
        elif len(etv) == 1:
            val = int(etv[0])
        else:
            val = int(etv[0])
        for bit_val, name in PSF_NAMES.items():
            if val & bit_val:
                return name
        return 'UNKNOWN'
    except Exception as e:
        return 'ERR(%s)' % str(e)[:20]

def met_to_mjd(met):
    return met / 86400.0 + FERMI_MJD_REF

def sizeof_fmt(num_bytes):
    for unit in ['B','KB','MB','GB']:
        if num_bytes < 1024.0:
            return "%.1f %s" % (num_bytes, unit)
        num_bytes /= 1024.0
    return "%.1f TB" % num_bytes

# ---- Read file list ----
if not os.path.exists(EVTLIST):
    print("ERROR: Event list not found: %s" % EVTLIST)
    sys.exit(1)

with open(EVTLIST) as f:
    filelist = [l.strip() for l in f if l.strip() and not l.startswith('#')]

print("=" * 75)
print("Event list : %s" % EVTLIST)
print("Files listed: %d" % len(filelist))
print("Energy threshold : %.0f MeV = %.0f GeV" % (EMIN_MEV, EMIN_MEV/1000.))
print("Search radius    : %.2f deg" % SEP_DEG)
print("Source position  : RA=%.4f, Dec=%.4f" % (RA_SRC, DEC_SRC))
print("=" * 75)

n_ok, n_fail = 0, 0
total_events = 0
all_high_e = []

print("\n--- File-by-file properties ---\n")
print("%-45s  %8s  %10s  %10s  %12s  %12s  %6s  %8s  %8s" % (
    'Filename', 'Size', 'Emin(MeV)', 'Emax(MeV)',
    'Start(MJD)', 'Stop(MJD)', 'Events', 'N>300GeV', 'Within0.5'))
print("-" * 130)

for fpath in filelist:
    fname = os.path.basename(fpath)

    # --- Check path correctness ---
    expected_prefix = '/sdf/data/fermi/u/adinesh'
    path_ok = fpath.startswith(expected_prefix)
    path_flag = '' if path_ok else '  <-- WRONG PATH (home dir?)'

    if not os.path.exists(fpath):
        print("  MISSING: %s%s" % (fpath, path_flag))
        n_fail += 1
        continue

    try:
        fsize = sizeof_fmt(os.path.getsize(fpath))
        h = fits.open(fpath, memmap=True)

        # --- File-level properties from header ---
        hdr = h['EVENTS'].header
        emin_file = hdr.get('DSVAL3', 'N/A')   # energy range stored in DSS keywords
        emax_file = hdr.get('DSVAL4', 'N/A')
        # Fallback: use actual min/max of ENERGY column
        evts = h['EVENTS'].data
        n_evts = len(evts)
        total_events += n_evts

        emin_actual = float(evts['ENERGY'].min()) if n_evts > 0 else 0.0
        emax_actual = float(evts['ENERGY'].max()) if n_evts > 0 else 0.0
        tstart_mjd  = met_to_mjd(float(hdr.get('TSTART', 0)))
        tstop_mjd   = met_to_mjd(float(hdr.get('TSTOP',  0)))

        # --- High-energy selection ---
        emask = evts['ENERGY'] > EMIN_MEV
        n_he  = emask.sum()
        n_within = 0

        if n_he > 0:
            he  = evts[emask]
            sep = angsep_deg(he['RA'], he['DEC'], RA_SRC, DEC_SRC)
            within = sep < SEP_DEG
            n_within = within.sum()

            if n_within > 0:
                for idx in np.where(within)[0]:
                    evt = he[idx]
                    psf = get_psf_type(evt['EVENT_TYPE'])
                    mjd = met_to_mjd(float(evt['TIME']))
                    all_high_e.append({
                        'energy_MeV' : float(evt['ENERGY']),
                        'ra'         : float(evt['RA']),
                        'dec'        : float(evt['DEC']),
                        'sep_deg'    : float(sep[idx]),
                        'psf_type'   : psf,
                        'mjd'        : mjd,
                        'file'       : fname,
                    })

        flag = '  <-- HAS CANDIDATES' if n_within > 0 else ''
        print("%-45s  %8s  %10.1f  %10.1f  %12.2f  %12.2f  %6d  %8d  %8d%s%s" % (
            fname, fsize,
            emin_actual, emax_actual,
            tstart_mjd, tstop_mjd,
            n_evts, n_he, n_within,
            flag, path_flag))

        h.close()
        n_ok += 1

    except Exception as e:
        print("  ERROR reading %s: %s" % (fpath, str(e)))
        n_fail += 1

# ---- Summary ----
print("\n" + "=" * 75)
print("FILES READ OK : %d / %d" % (n_ok, len(filelist)))
print("FILES FAILED  : %d / %d" % (n_fail, len(filelist)))
print("TOTAL EVENTS  : %d" % total_events)
print("=" * 75)

if not all_high_e:
    print("\nNo photons above %.0f GeV found within %.2f deg." % (EMIN_MEV/1000., SEP_DEG))
else:
    all_high_e.sort(key=lambda x: x['energy_MeV'], reverse=True)

    print("\n--- Photons above %.0f GeV within %.2f deg of source ---" % (EMIN_MEV/1000., SEP_DEG))
    print("%-10s  %-8s  %-9s  %-9s  %-10s  %-6s  %-45s" % (
        'E(GeV)', 'E(TeV)', 'RA', 'Dec', 'Sep(deg)', 'PSF', 'File'))
    print("%-10s  %-8s  %-10s  %-10s  %-12s  %-12s  %-6s  %-45s" % (
        '', '', '', '', '', 'MJD', '', ''))
    print("-" * 105)
    for e in all_high_e:
        print("%-10.1f  %-8.4f  %-9.4f  %-9.4f  %-10.5f  %-12.2f  %-6s  %-45s" % (
            e['energy_MeV']/1000., e['energy_MeV']/1e6,
            e['ra'], e['dec'], e['sep_deg'],
            e['mjd'], e['psf_type'], e['file']))

    print("\nTotal: %d photon(s) found" % len(all_high_e))

    print("\n--- Counts at different search radii ---")
    for r in [0.05, 0.10, 0.20, 0.30, 0.50]:
        n = sum(1 for e in all_high_e if e['sep_deg'] < r)
        print("  Within %.2f deg: %d" % (r, n))

    print("\n--- PSF type breakdown ---")
    for psf in ['PSF0', 'PSF1', 'PSF2', 'PSF3', 'UNKNOWN']:
        n = sum(1 for e in all_high_e if e['psf_type'] == psf)
        if n > 0:
            print("  %s: %d" % (psf, n))

    # --- Time coverage check ---
    print("\n--- Temporal coverage of all 6 files ---")
