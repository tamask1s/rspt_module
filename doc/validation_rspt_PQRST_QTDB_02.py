import wfdb
import numpy as np
import matplotlib.pyplot as plt
import rspt_module

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
def plot_annotations(signal, ann_by_type, rspt_annotations, fs):
    import matplotlib.pyplot as plt
    t = np.arange(len(signal)) / fs

    plt.figure(figsize=(15, 8))

    # --- 1. subplot: QTDB annotációk ---
    ax1 = plt.subplot(2, 1, 1)
    ax1.plot(t, signal, label='ECG (lead I)', linewidth=1)

    colors = {'p': 'g', 'N': 'b', 't': 'c'}
    markers = {'on': 'o', 'peak': 'x', 'off': '^'}

    for wave_type, triples in ann_by_type.items():
        for on, peak, off in triples:
            c = colors[wave_type]
            ax1.plot(on / fs, signal[on], c + markers['on'], label=f'{wave_type}_on')
            ax1.plot(peak / fs, signal[peak], c + markers['peak'], label=f'{wave_type}_peak')
            ax1.plot(off / fs, signal[off], c + markers['off'], label=f'{wave_type}_off')

    handles, labels = ax1.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax1.legend(by_label.values(), by_label.keys())
    ax1.set_title("QTDB annotációk")
    ax1.set_ylabel("mV")
    ax1.grid(True)

    # --- 2. subplot: RSPT annotációk ---
    ax2 = plt.subplot(2, 1, 2, sharex=ax1)
    ax2.plot(t, signal, label='ECG (lead I)', linewidth=1)

    rspt_map = {}
    for item in rspt_annotations:
        idx_str, symbol = item.split(':')
        idx = int(idx_str)
        if idx not in rspt_map:
            rspt_map[idx] = []
        rspt_map[idx].append(symbol)

    symbol_color_marker = {
        '(': ('g', 'o'),
        'p': ('g', 'x'),
        ')': ('g', '^'),
        '[': ('b', 'o'),
        'N': ('b', 'x'),
        ']': ('b', '^'),
        '{': ('c', 'o'),
        't': ('c', 'x'),
        '}': ('c', '^'),
    }

    for idx, symbols in rspt_map.items():
        for sym in symbols:
            if sym in symbol_color_marker:
                color, marker = symbol_color_marker[sym]
                ax2.plot(idx / fs, signal[idx], color + marker, label=sym)

    handles, labels = ax2.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax2.legend(by_label.values(), by_label.keys())
    ax2.set_title("RSPT annotációk")
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
    result = compute_rspt_like_stats(signal, ann_by_type, fs)
    print_result(result)

    # RSPT analízis meghívása
    rspt_result = rspt_module.analyse_ecg(np.stack([signal, signal], axis=1), fs)
    #print('annots1: ', ann_by_type)
    print('annots2: ', rspt_result['annotations'])
    plot_annotations(signal, ann_by_type, rspt_result['annotations'], fs)
    #rspt_result = rspt.analyze_ecg(ecg_signal=signal.tolist(), sampling_rate=fs, leads=[LEAD_IDX])
    print_result(rspt_result)

    # Eredmények összehasonlítása
    print("\n--- Összehasonlítás: QTDB vs RSPT ---")
    keys_to_compare = [
        "rr_interval_ms", "rr_variation_ms", "heart_rate_bpm", "pr_interval_ms",
        "qrs_duration_ms", "qt_interval_ms", "qtc_interval_ms",
        "p_wave_duration_ms", "t_wave_duration_ms"
    ]

    for key in keys_to_compare:
        own_val = result.get(key)
        rspt_val = rspt_result.get(key)
        diff = abs(own_val - rspt_val) if own_val is not None and rspt_val is not None else None
        print(f"{key:<25} QTDB: {own_val:<10.2f} | RSPT: {rspt_val:<10.2f} | Δ: {diff:.2f}" if diff is not None else f"{key:<25} Érték hiányzik.")
