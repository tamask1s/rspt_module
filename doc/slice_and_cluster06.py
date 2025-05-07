import wfdb
import numpy as np
import matplotlib.pyplot as plt
from rspt_module import detect_peaks
import neurokit2 as nk

def extract_segments(sig, peaks, wb, wa):
    segs, idx = [], []
    for p in peaks:
        if p-wb>=0 and p+wa<len(sig):
            segs.append(sig[p-wb:p+wa])
            idx.append(p)
    return np.array(segs), np.array(idx)

def slice_and_cluster(record_name, k=5):
    # 1) Betöltés
    rec = wfdb.rdrecord("/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"+record_name)
    sig, fs = rec.p_signal[:,0], rec.fs
    wb, wa = int(0.2*fs), int(0.4*fs)

    # 2) PPV szegmensek + SVD
    ppv = detect_peaks(sig, fs, mode="high_ppv")
    seg_ppv, _ = extract_segments(sig, ppv, wb, wa)
    mean_seg = seg_ppv.mean(axis=0)
    centered = seg_ppv - mean_seg
    U,S,Vt = np.linalg.svd(centered, full_matrices=False)
    comps = Vt[:k]   # k×L

    # 3) Sens szegmensek
    sens = detect_peaks(sig, fs, mode="high_sensitivity")
    seg_sens, idx_sens = extract_segments(sig, sens, wb, wa)

    final_peaks = []
    for seg, orig_idx in zip(seg_sens, idx_sens):
        # 4) reconstruct
        cen = seg - mean_seg
        coords = cen @ comps.T          # k-dim feature
        recon = coords @ comps + mean_seg

        # 5) újracsúcs detektálás a rekonstruált jelben
        r_new = detect_peaks(recon, fs, mode="default")
        # ha legalább egy csúcs van, megtartjuk az eredetit:
        if len(r_new)>0:
            final_peaks.append(orig_idx)

    final_peaks = np.array(sorted(set(final_peaks)))

    # 6) ábra
    #plt.figure(figsize=(20,4))
    #plt.plot(sig, label="ECG")
    #plt.scatter(final_peaks, sig[final_peaks], c="red", s=10, label="Filtered R")
    #plt.legend(); plt.grid(True); plt.tight_layout(); plt.show()

    return final_peaks

# Futtatás:
# peaks = slice_svd_reconstruct("100", k=5)
if __name__ == "__main__":
    slice_and_cluster("108", k=5)

