import wfdb
import numpy as np
from rspt_module import detect_peaks

def benchmark_record(record_name, tolerance=0.05):
    print(f"\nBenchmarking record: {record_name}")

    # Rekord és annotáció beolvasása PhysioNet-ről
    record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"

    record = wfdb.rdrecord(record_path + record_name)
    annotation = wfdb.rdann(record_path + record_name, 'atr')
#    record = wfdb.rdrecord(record_name, pn_dir='mitdb')
#    annotation = wfdb.rdann(record_name, 'atr', pn_dir='mitdb')

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

    records = [
    '100', '101', '102', '103', '104', '105', '106', '107',
    '108', '109', '111', '112', '113', '114', '115', '116',
    '117', '118', '119', '121', '122', '123', '124', '200',
    '201', '202', '203', '205', '207', '208', '209', '210',
    '212', '213', '214', '215', '217', '219', '220', '221',
    '222', '223', '228', '230', '231', '232', '233', '234']
    #records = ['100', '101', '103', '105', '111', '112', '113', '115', '117', '118', '119', '121', '122', '200', '202', '209']    # ide még bővíthetsz
    sensitivities = []
    ppvs = []
    results = []
    for rec in records:
        res = benchmark_record(rec)
        results.append(res)

    print("\n=== Summary ===")
    for r in results:
        print(f"{r['record']}: Sensitivity={r['Sensitivity']:.3f}, PPV={r['PPV']:.3f}")
        sensitivities.append(r['Sensitivity'])
        ppvs.append(r['PPV'])

    # Átlagos metrikák
    avg_sensitivity = np.mean(sensitivities)
    avg_ppv = np.mean(ppvs)

    print(f"Average Sensitivity (100% means no FN): {avg_sensitivity * 100:.3f}%")
    print(f"Average (100% means no FP) PPV: {avg_ppv * 100:.3f}%")
    print(f"Average SPPV: {(avg_ppv + avg_sensitivity) * 50:.3f}%")
    