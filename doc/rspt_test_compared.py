import numpy as np
import wfdb
import rspt_module
import matplotlib.pyplot as plt

print("opening")
record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/200"
#record_path = "/media/sf_SharedFolder/QT/lobachevsky-university-electrocardiography-database-1.0.1/data/117"
#record_path = "/media/sf_SharedFolder/QT/st-petersburg-incart-12-lead-arrhythmia-database-1.0.0/files/I54"
record = wfdb.rdrecord(record_path)
annotation = wfdb.rdann(record_path, 'atr')
#record = wfdb.rdrecord('101', pn_dir='mitdb')
#annotation = wfdb.rdann('101', 'atr', pn_dir='mitdb')
print("opened")

ecg_signal = record.p_signal[:, :]
sampling_rate = record.fs

print("Detecting. sampling_rate:", sampling_rate)
peak_indexes = np.array(rspt_module.detect_multichannel(ecg_signal, sampling_rate, 'default')) #'high_ppv'
print(f"{len(peak_indexes)} peaks detected.")
print("First 10 peak indexes:", peak_indexes[:10])

# Annotációból az R csúcsok
#annot_r_peaks = annotation.sample[np.isin(annotation.symbol, ['N', 'L', 'R', 'V', 'A', 'F', 'j', 'E', 'e', 'a', 'J', 'S'])]
annot_r_peaks = annotation.sample[np.isin(annotation.symbol, ['N', 'L', 'R', 'B', 'A', 'a', 'J', 'S', 'V', 'r', 'F', 'e', 'j', 'n', 'E', '/', 'f', 'Q', '?'])]

# TP / FP / FN meghatározása
tolerance_samples = int(0.05 * sampling_rate)

TP = []
FP = []
FN = []

matched_annot = np.full(len(annot_r_peaks), False)

for peak in peak_indexes:
    diffs = np.abs(annot_r_peaks - peak)
    min_diff = np.min(diffs)
    match_idx = np.argmin(diffs)
    if min_diff <= tolerance_samples and not matched_annot[match_idx]:
        TP.append(peak)
        matched_annot[match_idx] = True
    else:
        FP.append(peak)

for i, ann in enumerate(annot_r_peaks):
    if not matched_annot[i]:
        FN.append(ann)

print(f"TP: {len(TP)}, FP: {len(FP)}, FN: {len(FN)}")
TP_ = len(TP)
FP_ = len(FP)
FN_ = len(FN)
sensitivity = TP_ / (TP_ + FN_) if (TP_ + FN_) > 0 else 0
ppv = TP_ / (TP_ + FP_) if (TP_ + FP_) > 0 else 0
print(f"Sensitivity (100% means no FN): {sensitivity * 100:.3f}%")
print(f"PPV (100% means no FP): {ppv * 100:.3f}%")

# --- Közös ablak, 2 subplot ---
fig, axs = plt.subplots(2, 1, figsize=(18, 8), sharex=True)

# 1. subplot: Detektált TP, FP, FN
axs[0].plot(ecg_signal, label='ECG')
axs[0].plot(TP, ecg_signal[TP], 'go', label='True Positive')
axs[0].plot(FP, ecg_signal[FP], 'ro', label='False Positive')
axs[0].plot(FN, ecg_signal[FN], 'bo', label='False Negative')
axs[0].set_title('Detected Peaks (TP, FP, FN)')
axs[0].set_ylabel('Amplitude (mV)')
axs[0].legend()
axs[0].grid(True)

# 2. subplot: Annotációs R-csúcsok
axs[1].plot(ecg_signal, label='ECG')
axs[1].plot(annot_r_peaks, ecg_signal[annot_r_peaks], 'mo', label='Annotated R peaks')
axs[1].set_title('Annotated R Peaks')
axs[1].set_xlabel('Sample index')
axs[1].set_ylabel('Amplitude (mV)')
axs[1].legend()
axs[1].grid(True)

plt.tight_layout()
plt.show()
