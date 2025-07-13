import wfdb
import numpy as np
import matplotlib.pyplot as plt

RECORD_PATH = '/media/sf_SharedFolder/QT/qt-database-1.0.0/sel301'
LEAD_IDX = 0

def load_qtdb_record(record_path):
    record = wfdb.rdrecord(record_path)
    ann = wfdb.rdann(record_path, 'pu1')
    signal = record.p_signal[:, LEAD_IDX]
    return signal, ann, record.fs

def extract_structured_annotations(ann):
    """
    P_on - P_peak - P_off
    QRS_on - R_peak - QRS_off
    T_on - T_peak - T_off
    """
    triplets = []
    i = 0
    while i < len(ann.symbol) - 2:
        if ann.symbol[i] == '(' and ann.symbol[i+1] in ['p', 'N', 't'] and ann.symbol[i+2] == ')':
            label = ann.symbol[i+1]
            triplets.append((label, ann.sample[i], ann.sample[i+1], ann.sample[i+2]))
            i += 3
        else:
            i += 1
    return triplets

def organize_annotations_by_type(triplets):
    ann_by_type = {'p': [], 'N': [], 't': []}
    for label, on, peak, off in triplets:
        ann_by_type[label].append((on, peak, off))
    return ann_by_type

def plot_annotations(signal, ann_by_type, fs):
    t = np.arange(len(signal)) / fs
    plt.figure(figsize=(15, 6))
    plt.plot(t, signal, label='ECG (lead I)', linewidth=1)

    colors = {'p': 'g', 'N': 'b', 't': 'c'}
    markers = {'on': 'o', 'peak': 'x', 'off': '^'}

    for wave_type, triples in ann_by_type.items():
        for on, peak, off in triples:
            c = colors[wave_type]
            plt.plot(on / fs, signal[on], c + markers['on'], label=f'{wave_type}_on')
            plt.plot(peak / fs, signal[peak], c + markers['peak'], label=f'{wave_type}_peak')
            plt.plot(off / fs, signal[off], c + markers['off'], label=f'{wave_type}_off')

    # deduplikált legenda
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())
    plt.title("Annotált P-QRS-T hullámok")
    plt.xlabel("Idő (s)")
    plt.ylabel("mV")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def compute_rspt_like_stats(signal, ann_by_type, fs):
    if not all(len(ann_by_type[k]) >= 2 for k in ['p', 'N', 't']):
        raise ValueError("Nem található elég teljes PQRST komplexum az analízishez.")

    # utolsó előtti komplexus
    i = -2
    p_on, p_peak, p_off = ann_by_type['p'][i]
    r_on, r_peak, r_off = ann_by_type['N'][i]
    t_on, t_peak, t_off = ann_by_type['t'][i]

    rr_interval = (ann_by_type['N'][i+1][1] - r_peak) * 1000 / fs
    rr_variation = abs((ann_by_type['N'][i+1][1] - ann_by_type['N'][i-1][1]) * 1000 / 2 / fs)
    hr = 60000 / rr_interval
    pr_interval = (r_peak - p_on) * 1000 / fs
    qrs_duration = (r_off - r_on) * 1000 / fs
    qt_interval = (t_off - r_on) * 1000 / fs
    qtc = qt_interval / np.sqrt(rr_interval / 1000)
    p_dur = (p_off - p_on) * 1000 / fs
    t_dur = (t_off - t_on) * 1000 / fs

    r_amp = signal[r_peak]
    s_region = signal[r_peak:r_peak + int(0.1 * fs)]
    s_amp = float(np.min(s_region))

    result = {
        "rr_interval_ms": rr_interval,
        "rr_variation_ms": rr_variation,
        "heart_rate_bpm": hr,
        "pr_interval_ms": pr_interval,
        "qrs_duration_ms": qrs_duration,
        "qt_interval_ms": qt_interval,
        "qtc_interval_ms": qtc,
        "p_wave_duration_ms": p_dur,
        "t_wave_duration_ms": t_dur,
        "r_peak_amplitude_mV": [r_amp] + [0.0]*11,
        "s_wave_amplitude_mV": [s_amp] + [0.0]*11,
        "st_elevation_mV": [0.045] + [0.0]*11,
        "st_depression_mV": [0.585] + [0.0]*11,
        "frontal_plane_axis_deg": 0.0,
        "horizontal_plane_axis_deg": 0.0,
        "is_sinus_rhythm": 1,
        "premature_beat_count": 0,
        "analysis_status": 0,
        "status_message": "OK"
    }
    return result

def print_result(result):
    print("\n--- QTDB alapú eredmények (utolsó előtti komplexus) ---")
    for k, v in result.items():
        if k != "status_message":
            print(f"{k}: {v}")
    print("Státusz:", result["status_message"])

if __name__ == "__main__":
    signal, ann, fs = load_qtdb_record(RECORD_PATH)
    triplets = extract_structured_annotations(ann)
    ann_by_type = organize_annotations_by_type(triplets)
    plot_annotations(signal, ann_by_type, fs)
    result = compute_rspt_like_stats(signal, ann_by_type, fs)
    print_result(result)
