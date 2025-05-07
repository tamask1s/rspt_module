import wfdb
import numpy as np
import matplotlib.pyplot as plt
from rspt_module import detect_peaks
import neurokit2 as nk


def extract_segments(signal, peaks, window_before, window_after):
    segments = []
    indices = []
    for peak in peaks:
        start = peak - window_before
        end = peak + window_after
        if start >= 0 and end < len(signal):
            segments.append(signal[start:end])
            indices.append(peak)
    return np.array(segments), np.array(indices)


def slice_and_cluster(record_name):
    # 1. Rekord betöltése
    record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
    record = wfdb.rdrecord(record_path + record_name)
    signal = record.p_signal[:, 0]
    fs = record.fs

    # 2. High-PPV csúcsdetektálás és szegmensek kivágása
    r_peaks_ppv = detect_peaks(signal, fs, mode="high_ppv")
    wb = int(0.2 * fs)
    wa = int(0.4 * fs)
    segments_ppv, idx_ppv = extract_segments(signal, r_peaks_ppv, wb, wa)

    # 3. Szegmensek középértékének levonása (centerelés)
    mean_seg = np.mean(segments_ppv, axis=0)
    centered = segments_ppv - mean_seg

    # 4. SVD dekompozíció
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)

    display = True
    # 5. Első három komponens megjelenítése
    if display:
        plt.figure(figsize=(10, 6))
        time = np.arange(segments_ppv.shape[1]) / fs * 1000  # ms
        for i in range(6):
            component = Vt[i] * S[i]
            plt.plot(time, component, label=f'PC{i+1} (singular value {S[i]:.1f})')
        plt.title("Első 3 főkomponens (PC) hullámformái")
        plt.xlabel("Idő (ms)")
        plt.ylabel("Amplitúdó")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        #plt.show()

    # 6. High-sensitivity detektálás és szegmensek
    r_peaks_sens = detect_peaks(signal, fs, mode="high_sensitivity")
    segments_sens, idx_sens = extract_segments(signal, r_peaks_sens, wb, wa)

    # 7. Korreláció PC1-gyel mint template
    template = Vt[0] * S[0]
    final_peaks = []
    for seg, idx in zip(segments_sens, idx_sens):
        corr = np.corrcoef(seg - mean_seg, template)[0, 1]
        if corr > 0.90:
            final_peaks.append(idx)
    final_peaks = np.array(sorted(set(final_peaks)))

    # 8. Eredeti jel és végső csúcsok megjelenítése
    if display:
        plt.figure(figsize=(20, 5))
        plt.plot(signal, label="ECG")
        plt.scatter(final_peaks, signal[final_peaks], color='red', s=10, label="Final R peaks")
        plt.title(f"Végső R-csúcsok (SVD alapú)")
        plt.xlabel("Mintavételi pont")
        plt.ylabel("Feszültség (mV)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return final_peaks

if __name__ == "__main__":
    slice_and_cluster("200")
