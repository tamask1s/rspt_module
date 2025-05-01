import wfdb
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import neurokit2 as nk


def detect_peaks(signal, fs, mode="high_sensitivity"):
    cleaned = nk.ecg_clean(signal, sampling_rate=fs)
    _, info = nk.ecg_peaks(cleaned, sampling_rate=fs)
    return info["ECG_R_Peaks"]


def extract_segments(signal, peaks, window_before, window_after):
    segments = []
    indices = []
    for peak in peaks:
        start = peak - window_before
        end = peak + window_after
        if start >= 0 and end < len(signal):
            segment = signal[start:end]
            segments.append(segment)
            indices.append(peak)
    return np.array(segments), np.array(indices)


def slice_and_cluster(record_name, nr_max_clusters=3):
    # --- 1. Rekord beolvasása ---
    record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
    record = wfdb.rdrecord(record_path + record_name)
    signal = record.p_signal[:, 0]
    fs = record.fs

    # --- 2. High-PPV R-csúcsok detektálása ---
    r_peaks_ppv = detect_peaks(signal, fs, mode="high_ppv")
    window_before = int(0.2 * fs)
    window_after = int(0.4 * fs)
    segments_ppv, indices_ppv = extract_segments(signal, r_peaks_ppv, window_before, window_after)

    # --- 3. Normalizálás és klaszterezés ---
    scaler = StandardScaler()
    segments_scaled = scaler.fit_transform(segments_ppv)

    best_score = -1
    best_n = 1
    best_labels = None
    best_model = None
    for n_clusters in range(1, nr_max_clusters + 1):
        model = KMeans(n_clusters=n_clusters, random_state=42).fit(segments_scaled)
        if n_clusters > 1:
            score = silhouette_score(segments_scaled, model.labels_)
            if score > best_score:
                best_score = score
                best_n = n_clusters
                best_labels = model.labels_
                best_model = model
        else:
            best_labels = model.labels_
            best_model = model

    print(f"Legjobb klaszterszám: {best_n}, Silhouette score: {best_score:.4f}")

    # --- 4. Klaszterek ábrázolása ---
    for cluster_id in range(best_n):
        cluster_segments = segments_ppv[best_labels == cluster_id]
        plt.figure(figsize=(10, 4))
        for seg in cluster_segments:
            plt.plot(seg, alpha=0.3, color='gray')
        mean_waveform = np.mean(cluster_segments, axis=0)
        plt.plot(mean_waveform, color='red', linewidth=2, label=f'Cluster {cluster_id} átlag')
        plt.title(f"Cluster {cluster_id} ({len(cluster_segments)} szegmens)")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()

    # --- 5. High-sensitivity detektálás ---
    r_peaks_sens = detect_peaks(signal, fs, mode="high_sensitivity")
    segments_sens, indices_sens = extract_segments(signal, r_peaks_sens, window_before, window_after)

    # --- 6. Klaszterátlagokkal korreláció számítása ---
    final_peaks = []
    cluster_means = [np.mean(segments_ppv[best_labels == i], axis=0) for i in range(best_n)]

    for i, seg in enumerate(segments_sens):
        for mean in cluster_means:
            corr = np.corrcoef(seg, mean)[0, 1]
            if corr > 0.90:  # hasonlít valamelyik cluster átlagra
                final_peaks.append(indices_sens[i])
                break

    final_peaks = np.array(sorted(set(final_peaks)))

    # --- 7. Eredmény megjelenítése ---
    plt.figure(figsize=(20, 5))
    plt.plot(signal, label="ECG")
    plt.scatter(final_peaks, signal[final_peaks], color='red', s=10, label="Végső R-csúcsok")
    plt.title("Végső R-csúcsok detektálása")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return final_peaks

if __name__ == "__main__":
    slice_and_cluster("101", nr_max_clusters=3)

# Futtatás:
# slice_and_cluster("100", nr_max_clusters=3)
