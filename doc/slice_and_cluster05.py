import wfdb
import numpy as np
import matplotlib.pyplot as plt
from rspt_module import detect_peaks
from sklearn.metrics.pairwise import cosine_similarity

# Parameters
record_name = "100"
base_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
rec = wfdb.rdrecord(base_path + record_name)
sig = rec.p_signal[:,0]
fs = rec.fs
wb = int(0.2 * fs)
wa = int(0.4 * fs)

# 1. High-PPV segments
peaks_ppv = detect_peaks(sig, fs, mode="high_ppv")
segments_ppv = []
for peak in peaks_ppv:
    start = peak - wb
    end = peak + wa
    if 0 <= start and end < len(sig):
        segments_ppv.append(sig[start:end])
segments_ppv = np.array(segments_ppv)

# 2. Centering and SVD
mean_seg = np.mean(segments_ppv, axis=0)
centered_ppv = segments_ppv - mean_seg
U, S, Vt = np.linalg.svd(centered_ppv, full_matrices=False)
coords_ppv = (U[:, :3] * S[:3])

# 3. High-sensitivity segments
peaks_sens = detect_peaks(sig, fs, mode="high_sensitivity")
segments_sens = []
idx_sens = []
for peak in peaks_sens:
    start = peak - wb
    end = peak + wa
    if 0 <= start and end < len(sig):
        segments_sens.append(sig[start:end])
        idx_sens.append(peak)
segments_sens = np.array(segments_sens)
idx_sens = np.array(idx_sens)
centered_sens = segments_sens - mean_seg
coords_sens = centered_sens @ Vt[:3].T

# 4. Threshold sweep
thresholds = np.linspace(0.5, 0.99, 50)
counts = []
for thr in thresholds:
    sim = cosine_similarity(coords_sens, coords_ppv)
    max_sim = sim.max(axis=1)
    keep_mask = max_sim >= thr
    counts.append(np.sum(keep_mask))

# Plot threshold vs number of peaks kept
plt.figure()
plt.plot(thresholds, counts)
plt.xlabel("Cosine similarity threshold")
plt.ylabel("Number of R-peaks kept")
plt.title("Threshold sweep for cosine similarity")
plt.grid(True)
plt.tight_layout()
plt.show()
