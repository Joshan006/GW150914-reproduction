"""
GW150914 Reproduction — recovering the first LIGO binary-black-hole detection
from public GWOSC strain data.

This is a REPRODUCTION of the published result (Abbott et al. 2016), not an
original detection. It recovers the known signal from public data to confirm
the analysis method end to end.

Needs internet access to download the data from GWOSC on the first run.
Requires: GWpy (see environment.yml; on Windows install with conda, not pip)

Author: Joshanraj C V
Data:   Gravitational Wave Open Science Center (https://gwosc.org)
Based on the GWOSC "Signal Processing with GW150914" tutorial approach.
"""

import os
import matplotlib
matplotlib.use("Agg")            # so it saves figures without needing a screen
import matplotlib.pyplot as plt
import numpy as np
from astropy.utils.data import download_file
from gwosc.locate import get_event_urls
from gwpy.timeseries import TimeSeries

os.makedirs("figures", exist_ok=True)

# ----------------------------------------------------------------------
# 0. Event parameters
# ----------------------------------------------------------------------
# GW150914 merger time, in GPS seconds (published merger time ~1126259462.42;
# the rounded catalogue value .4 puts the peak ~0.02 s after t = 0).
GPS_EVENT = 1126259462.42
DETECTORS = ["H1", "L1"]         # Hanford, Livingston
TIMEOUT   = 120                  # seconds to wait on a slow download
WINDOW    = 16                   # seconds of data to keep on each side

# ----------------------------------------------------------------------
# 1. Fetch the public strain data
# ----------------------------------------------------------------------
# We ask GWOSC for the GW150914 data URLs and download with a generous
# timeout. GWOSC may serve the full 4096-second O1 file (~120 MB) rather than
# a short event file; cache=True means it is only downloaded once.
# 4096 Hz is plenty for GW150914 (the signal is below ~350 Hz).
#
# IMPORTANT: a 4096 s file can contain stretches where a detector was not
# recording science-quality data. GWOSC stores those stretches as NaN. A single
# NaN poisons the spectrum, whitening and Q-transform, so we keep only a short
# clean window (+/- 16 s) around the event, when both detectors were observing.
strain = {}
for det in DETECTORS:
    print(f"Fetching {det} ...")
    urls = get_event_urls("GW150914", detector=det, sample_rate=4096)
    url = next((u for u in urls if u.endswith("-32.hdf5")), urls[0])
    print(f"  {url}")
    path = download_file(url, cache=True, timeout=TIMEOUT, show_progress=True)
    ts = TimeSeries.read(path, format="hdf5.gwosc")
    ts = ts.crop(GPS_EVENT - WINDOW, GPS_EVENT + WINDOW)
    n_bad = int(np.count_nonzero(~np.isfinite(ts.value)))
    if n_bad:
        raise SystemExit(f"{det}: {n_bad} NaN/inf samples remain in the "
                         f"+/-{WINDOW} s window - cannot continue.")
    strain[det] = ts
    print(f"  {det}: {ts.duration} at {ts.sample_rate} (no NaNs)")

# ----------------------------------------------------------------------
# 2. Amplitude spectral density (the detector noise curve)
# ----------------------------------------------------------------------
# The ASD shows how much noise sits at each frequency. It is what "whitening"
# in step 3 divides out. Reproducing this plot is itself a good sanity check:
# you should see the seismic wall at low f, the thermal/quantum noise, and the
# narrow "violin" resonance spikes from the mirror suspensions.
plt.figure(figsize=(9, 5))
for det in DETECTORS:
    asd = strain[det].asd(fftlength=4)
    plt.loglog(asd.frequencies.value, asd.value, label=det)
plt.xlim(10, 2000)
plt.ylim(1e-24, 1e-19)
plt.xlabel("Frequency [Hz]")
plt.ylabel(r"ASD [strain / $\sqrt{\mathrm{Hz}}$]")
plt.title("GW150914 — detector amplitude spectral density")
plt.legend()
plt.grid(True, which="both", alpha=0.3)
plt.tight_layout()
plt.savefig("figures/fig1_asd.png", dpi=130)
print("saved fig1_asd.png")

# ----------------------------------------------------------------------
# 3. Whiten + band-pass  ->  the chirp
# ----------------------------------------------------------------------
# Whitening flattens the noise so every frequency contributes comparably.
# Band-passing 35-350 Hz keeps the band where the chirp lives and removes the
# low-frequency seismic wall and high-frequency junk.
BANDPASS = (35, 350)

whitened = {}
for det in DETECTORS:
    w = strain[det].whiten().bandpass(*BANDPASS)
    whitened[det] = w

# Zoom to a ~0.2 s window around the merger, where the chirp is visible.
t1, t2 = GPS_EVENT - 0.15, GPS_EVENT + 0.05

plt.figure(figsize=(9, 5))
for det in DETECTORS:
    seg = whitened[det].crop(t1, t2)
    plt.plot(seg.times.value - GPS_EVENT, seg.value, label=det, alpha=0.8)
plt.xlabel("Time since merger [s]")
plt.ylabel("Whitened strain")
plt.title("GW150914 — whitened, band-passed strain (the chirp)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figures/fig2_chirp.png", dpi=130)
print("saved fig2_chirp.png")

# ----------------------------------------------------------------------
# 3b. Overlay H1 and L1 aligned  (the "aha" comparison)
# ----------------------------------------------------------------------
# The two detectors are oriented differently and are ~10 ms apart in light
# travel time. For GW150914 the signal reached L1 ~6.9 ms BEFORE H1, so L1
# lines up with H1 after shifting it LATER by ~6.9 ms and flipping the sign.
# Seeing them overlap is the core of the reproduction.
L1_SHIFT = 0.0069   # seconds; L1 arrived ~6.9 ms before H1
h1 = whitened["H1"].crop(t1, t2)
l1 = whitened["L1"].crop(t1, t2)

plt.figure(figsize=(9, 5))
plt.plot(h1.times.value - GPS_EVENT, h1.value, label="H1")
plt.plot(l1.times.value - GPS_EVENT + L1_SHIFT, -l1.value,
         label="L1 (shifted +6.9 ms, inverted)", alpha=0.8)
plt.xlabel("Time since merger [s]")
plt.ylabel("Whitened strain")
plt.title("GW150914 — H1 vs L1 aligned")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("figures/fig3_h1_l1_aligned.png", dpi=130)
print("saved fig3_h1_l1_aligned.png")

# ----------------------------------------------------------------------
# 4. Q-transform spectrogram  (the iconic figure)
# ----------------------------------------------------------------------
# The time-frequency map shows the upward-sweeping arc: frequency rising from
# ~35 Hz to ~250 Hz over the final ~0.2 s as the black holes spiral in.
for det in DETECTORS:
    q = strain[det].q_transform(outseg=(GPS_EVENT - 0.3, GPS_EVENT + 0.1),
                                qrange=(4, 64), frange=(20, 400))
    plt.figure(figsize=(9, 5))
    plt.pcolormesh(q.times.value - GPS_EVENT, q.frequencies.value, q.value.T,
                   shading="auto")
    plt.yscale("log")
    plt.ylim(20, 400)
    plt.xlabel("Time since merger [s]")
    plt.ylabel("Frequency [Hz]")
    plt.title(f"GW150914 — {det} Q-transform spectrogram")
    plt.colorbar(label="Normalised energy")
    plt.tight_layout()
    plt.savefig(f"figures/fig4_qscan_{det}.png", dpi=130)
    print(f"saved fig4_qscan_{det}.png")

print("\nDone. You should have: fig1_asd, fig2_chirp, fig3_h1_l1_aligned, "
      "fig4_qscan_H1, fig4_qscan_L1.")
print("Compare fig2/fig3 (the chirp) and fig4 (the arc) against the published "
      "GW150914 figures.")
