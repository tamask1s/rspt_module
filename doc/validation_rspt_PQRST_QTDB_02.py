import wfdb
import numpy as np
import matplotlib.pyplot as plt
import rspt_module

RECORD_PATH = '/media/sf_SharedFolder/QT/qt-database-1.0.0/sel301'
RECORD_PATH = '/media/sf_SharedFolder/QT/qt-database-1.0.0/sel306'
DISPLAY_LEAD_IDX = 1

def load_qtdb_record(record_path):
    record = wfdb.rdrecord(record_path)
    ann = wfdb.rdann(record_path, 'pu1')
    #signal = record.p_signal[:, DISPLAY_LEAD_IDX]
    signal = record.p_signal
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


def plot_annotations(signal, ann_by_type, rspt_annotations, fs):
    import matplotlib.pyplot as plt
    import numpy as np

    # --- 1. Keressük meg az utolsó előtti komplexust (minden típusból kell legalább 2) ---
    if not all(len(ann_by_type[k]) >= 2 for k in ['p', 'N', 't']):
        raise ValueError("Nem található elég teljes PQRST komplexum a megjelenítéshez.")

    idx = -1  # utolsó előtti
    p_on, p_peak, p_off = ann_by_type['p'][idx]
    r_on, r_peak, r_off = ann_by_type['N'][idx]
    t_on, t_peak, t_off = ann_by_type['t'][idx]

    # --- 2. Kivágási ablak meghatározása ±200 ms pufferral ---
    pad = int(0.2 * fs)
    start = max(0, min(p_on, r_on, t_on) - pad)
    end = min(len(signal), max(p_off, r_off, t_off) + pad)

    t = np.arange(start, end) / fs
    sig_segment = signal[start:end]

    plt.figure(figsize=(15, 8))

    # --- 1. subplot: QTDB annotációk ---
    ax1 = plt.subplot(2, 1, 1)
    ax1.plot(t, sig_segment, label='ECG (lead I)', linewidth=1)

    colors = {'p': 'g', 'N': 'b', 't': 'c'}
    markers = {'on': 'o', 'peak': 'x', 'off': '^'}

    for wave_type, triples in ann_by_type.items():
        for on, peak, off in triples:
            if start <= on <= end:
                c = colors[wave_type]
                ax1.plot(on / fs, signal[on], c + markers['on'], label=f'{wave_type}_on')
                ax1.plot(peak / fs, signal[peak], c + markers['peak'], label=f'{wave_type}_peak')
                ax1.plot(off / fs, signal[off], c + markers['off'], label=f'{wave_type}_off')

    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax1.legend(by_label.values(), by_label.keys())
    ax1.set_title("QTDB annotációk – utolsó előtti komplexus")
    ax1.set_ylabel("mV")
    ax1.grid(True)

    # --- 2. subplot: RSPT annotációk ---
    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    ax2.plot(t, sig_segment, label='ECG (lead I)', linewidth=1)

    colors = {'p': 'g', 'r': 'b', 't': 'c'}
    markers = {0: 'o', 1: 'x', 2: '^'}  # on, peak, off

    for ann in rspt_annotations:
        for wave_type in ['p', 'r', 't']:
            if wave_type in ann:
                for pos_idx, sample_idx in enumerate(ann[wave_type]):
                    if sample_idx >= 0 and start <= sample_idx <= end:
                        c = colors[wave_type]
                        m = markers[pos_idx]
                        ax2.plot(sample_idx / fs, signal[sample_idx], c + m, label=f'{wave_type}_{["on","peak","off"][pos_idx]}')

    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax2.legend(by_label.values(), by_label.keys())
    ax2.set_title("RSPT annotációk – utolsó előtti komplexus")
    ax2.set_xlabel("Idő (s)")
    ax2.set_ylabel("mV")
    ax2.grid(True)

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
    pr_interval = (r_on - p_off) * 1000 / fs
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
    result = compute_rspt_like_stats(signal[:, DISPLAY_LEAD_IDX], ann_by_type, fs)
    print_result(result)

    # RSPT analízis meghívása
    #rspt_result = rspt_module.analyse_ecg(np.stack([signal, signal], axis=1), fs)
    rspt_result = rspt_module.analyse_ecg(signal, fs)
    #print('annots1: ', ann_by_type)
    print('annots2: ', rspt_result['annotations'])
    plot_annotations(signal[:, DISPLAY_LEAD_IDX], ann_by_type, rspt_result['annotations'], fs)
    #rspt_result = rspt.analyze_ecg(ecg_signal=signal.tolist(), sampling_rate=fs, leads=[LEAD_IDX])
    print_result(rspt_result)

    # Eredmények összehasonlítása
    print("\n--- Összehasonlítás: QTDB vs RSPT ---")
    keys_to_compare = [
        "rr_interval_ms", "rr_variation_ms", "heart_rate_bpm", "pr_interval_ms",
        "qrs_duration_ms", "qt_interval_ms", "qtc_interval_ms",
        "p_wave_duration_ms", "t_wave_duration_ms"
    ]

#    for key in keys_to_compare:
#        own_val = result.get(key)
#        rspt_val = rspt_result.get(key)
#        diff = abs(own_val - rspt_val) if own_val is not None and rspt_val is not None else None
#        print(f"{key:<25} QTDB: {own_val:<10.2f} | RSPT: {rspt_val:<10.2f} | Δ: {diff:.2f}" if diff is not None else f"{key:<25} Érték hiányzik.")

for key in keys_to_compare:
    own_val = result.get(key)
    rspt_val = rspt_result.get(key)

    if own_val is not None and rspt_val is not None:
        diff = abs(own_val - rspt_val)
        pct_diff = (diff / rspt_val * 100) if rspt_val != 0 else float('nan')
        print(f"{key:<25} QTDB: {own_val:<10.2f} | RSPT: {rspt_val:<10.2f} | Δ: {diff:<8.2f} | Δ%: {pct_diff:.2f}%")
    else:
        print(f"{key:<25} Érték hiányzik.")
