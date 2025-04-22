import wfdb
import numpy as np
from rspt_module import detect_peaks

def benchmark_record(record_name, tolerance=0.05):
    print(f"\nBenchmarking record: {record_name}")

    # Rekord és annotáció beolvasása PhysioNet-ről
    record = wfdb.rdrecord(record_name, pn_dir='mitdb')
    annotation = wfdb.rdann(record_name, 'atr', pn_dir='mitdb')

    fs = record.fs
    signal = record.p_signal[:, 0]  # csak az első csatorna

    # Detektálás
    detected_peaks = detect_peaks(signal, fs)
    detected_peaks = np.array(detected_peaks)

    # Annotált R csúcsok kigyűjtése
    annot_r_peaks = annotation.sample[np.isin(annotation.symbol, ['N', 'L', 'R', 'A', 'V'])]

    # Tolerancia ablak (mintákban)
    tolerance_samples = int(tolerance * fs)

    TP, FP, FN = 0, 0, 0
    matched = set()

    for peak in detected_peaks:
        diffs = np.abs(annot_r_peaks - peak)
        if np.any(diffs <= tolerance_samples):
            match_idx = np.argmin(diffs)
            if match_idx not in matched:
                TP += 1
                matched.add(match_idx)
            else:
                FP += 1  # már volt hozzá találat
        else:
            FP += 1

    FN = len(annot_r_peaks) - len(matched)

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    ppv = TP / (TP + FP) if (TP + FP) > 0 else 0

    print(f"  True Positives (TP): {TP}")
    print(f"  False Positives (FP): {FP}")
    print(f"  False Negatives (FN): {FN}")
    print(f"  Sensitivity (Recall): {sensitivity:.3f}")
    print(f"  Positive Predictive Value (PPV): {ppv:.3f}")

    return {
        'record': record_name,
        'TP': TP,
        'FP': FP,
        'FN': FN,
        'Sensitivity': sensitivity,
        'PPV': ppv
    }

if __name__ == "__main__":
    records = ['100', '101', '103', '105', '111', '112', '113', '115', '117', '118', '119', '121', '122', '200', '202', '209']    # ide még bővíthetsz

    results = []
    for rec in records:
        res = benchmark_record(rec)
        results.append(res)

    print("\n=== Summary ===")
    for r in results:
        print(f"{r['record']}: Sensitivity={r['Sensitivity']:.3f}, PPV={r['PPV']:.3f}")
