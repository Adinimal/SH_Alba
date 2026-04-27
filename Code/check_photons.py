import astropy.io.fits as fits
import numpy as np

# Point to your event files
ra_src, dec_src = 80.737, -36.4686  # J0522 coordinates
emin = 300000  # MeV

all_events = []
for f in open('/sdf/data/fermi/u/adinesh/DATA/J0522/events_j0522.txt').read().splitlines():
    h = fits.open(f)
    evts = h['EVENTS'].data
    mask = (evts['ENERGY'] > emin)
    if mask.sum() > 0:
        # Angular separation cut: ~0.2 deg at these energies is generous
        dra  = (evts['RA'][mask]  - ra_src) * np.cos(np.radians(dec_src))
        ddec =  evts['DEC'][mask] - dec_src
        sep  = np.sqrt(dra**2 + ddec**2)
        all_events.append(evts[mask][sep < 0.5])
    h.close()

if all_events:
    combined = np.concatenate(all_events)
    print("Events above 300 GeV within 0.5 deg: %d" % len(combined))
    print("Max energy event: %.1f MeV = %.2f TeV" % (combined['ENERGY'].max(), combined['ENERGY'].max()/1e6))
    print("Energies (GeV):", sorted(combined['ENERGY']/1000., reverse=True)[:20])
