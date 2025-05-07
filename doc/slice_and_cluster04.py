import wfdb
import numpy as np
import matplotlib.pyplot as plt
from rspt_module import detect_peaks
from sklearn.metrics.pairwise import cosine_similarity

def extract_segments(signal, peaks, window_before, window_after):
    segments = []
    indices = []
    for peak in peaks:
        start = peak - window_before
        end   = peak + window_after
        if start >= 0 and end < len(signal):
            segments.append(signal[start:end])
            indices.append(peak)
    return np.array(segments), np.array(indices)

def slice_and_cluster(record_name, n_components=3, sim_threshold=0.99):
    # 1) Betöltés
    base_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
    rec = wfdb.rdrecord(base_path + record_name)
    sig = rec.p_signal[:,0]
    fs  = rec.fs

    # 2) High-PPV R-csúcsok és szegmensek
    peaks_ppv = detect_peaks(sig, fs, mode="high_ppv")
    wb = int(0.2 * fs)
    wa = int(0.4 * fs)
    seg_ppv, idx_ppv = extract_segments(sig, peaks_ppv, wb, wa)

    # 3) Centering: átlaglevonás
    mean_seg = np.mean(seg_ppv, axis=0)
    centered_ppv = seg_ppv - mean_seg

    # 4) SVD
    U, S, Vt = np.linalg.svd(centered_ppv, full_matrices=False)

    display = False
    # 5) PC hullámformák kirajzolása
    if display:
        plt.figure(figsize=(8,5))
        t = np.arange(seg_ppv.shape[1]) / fs * 1000  # ms tengely
        for i in range(n_components):
            pc = Vt[i] * S[i]
            plt.plot(t, pc, label=f'PC{i+1} (σ={S[i]:.1f})')
        plt.title("Első 3 főkomponens hullámformái")
        plt.xlabel("Idő (ms)")
        plt.ylabel("Amplitúdó")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        #plt.show()

    # 6) PPV-s koordináták (feature-tér)
    coords_ppv = (U[:, :n_components] * S[:n_components])

    # 7) High-sensitivity R-csúcsok és szegmensek
    peaks_sens = detect_peaks(sig, fs, mode="high_sensitivity")
    seg_sens, idx_sens = extract_segments(sig, peaks_sens, wb, wa)
    centered_sens = seg_sens - mean_seg

    # 8) Projekció a PC-térbe
    coords_sens = centered_sens @ Vt[:n_components].T

    # 9) Cosine-similarity és szűrés
    sim = cosine_similarity(coords_sens, coords_ppv)
    max_sim = sim.max(axis=1)
    keep_mask = max_sim >= sim_threshold
    final_peaks = idx_sens[keep_mask]

    # 10) Végső csúcsok kirajzolása
    if display:
        plt.figure(figsize=(16,4))
        plt.plot(sig, color='black', linewidth=0.8, label="ECG")
        plt.scatter(final_peaks, sig[final_peaks], color='red', s=10, label="Final R peaks")
        plt.title(f"Végső R-csúcsok (SVD+cosine ≥ {sim_threshold})")
        plt.xlabel("Mintavételi index")
        plt.ylabel("Feszültség (mV)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return np.sort(final_peaks)

if __name__ == "__main__":
    slice_and_cluster("108")
