import wfdb
import time
import numpy as np
from slice_and_cluster06 import slice_and_cluster

start_time2 = time.perf_counter()
end_time2 = time.perf_counter()
elapsed_ms2 = (end_time2 - start_time2) * 1000

def benchmark_record(record_name, tolerance=0.05):
    global elapsed_ms2

    record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
    record = wfdb.rdrecord(record_path + record_name)
    annotation = wfdb.rdann(record_path + record_name, 'atr')

    fs = record.fs
    signal = record.p_signal[:, 0]

    start_time2 = time.perf_counter()
    detected_peaks = slice_and_cluster(record_name);
    end_time2 = time.perf_counter()
    elapsed_ms2 = elapsed_ms2 + (end_time2 - start_time2) * 1000

    annot_r_peaks = annotation.sample[np.isin(annotation.symbol, ['N', 'L', 'R', 'V', 'A', 'F', 'j', 'E', 'e', 'a', 'J', 'S'])]

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
    #TP = []
    #FP = []
    #FN = []
    #matched_annot = np.full(len(annot_r_peaks), False)
    #for peak in detected_peaks:
    #    diffs = np.abs(annot_r_peaks - peak)
    #    min_diff = np.min(diffs)
    #    match_idx = np.argmin(diffs)
    #    if min_diff <= tolerance_samples and not matched_annot[match_idx]:
    #        TP.append(peak)
    #        matched_annot[match_idx] = True
    #    else:
    #        FP.append(peak)
    #
    #for i, ann in enumerate(annot_r_peaks):
    #    if not matched_annot[i]:
    #        FN.append(ann)
    #TP = len(TP)
    #FP = len(FP)
    #FN = len(FN)

    sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0
    ppv = TP / (TP + FP) if (TP + FP) > 0 else 0

    #print(f"  True Positives (TP): {TP}")
    #print(f"  False Positives (FP): {FP}")
    #print(f"  False Negatives (FN): {FN}")
    #print(f"  Sensitivity (Recall): {sensitivity * 100.0:.3f}")
    #print(f"  Positive Predictive Value (PPV): {ppv * 100.0:.3f}")

    return {
        'record': record_name,
        'TP': TP,
        'FP': FP,
        'FN': FN,
        'Sensitivity': sensitivity,
        'PPV': ppv
    }

if __name__ == "__main__":

    records = ['100', '101', '102', '103', '104', '105', '106', '107','108', '109', '111', '112', '113', '114', '115', '116', '117', '118', '119', '121', '122', '123', '124', '200', '201', '202', '203', '205', '207', '208', '209', '210', '212', '213', '214', '215', '217', '219', '220', '221', '222', '223', '228', '230', '231', '232', '233', '234']
    records = ['100', '101',        '103',        '105', '106',       '108', '109', '111', '112', '113', '114', '115', '116', '117', '118', '119', '121', '122', '123', '124', '200', '201', '202',        '205',               '209', '210', '212', '213', '214', '215',        '219', '220', '221', '222', '223',        '230', '231', '232', '233', '234']
    records = ['108', '114', '222', '233', '117', '214', '113', '115']
    #records = ['100', '101', '103', '105', '111', '112', '113', '115', '117', '118', '119', '121', '122', '200', '202', '209']    # ide még bővíthetsz
    sensitivities = []
    ppvs = []
    results = []
    start_time = time.perf_counter()
    for rec in records:
        res = benchmark_record(rec)
        results.append(res)
    end_time = time.perf_counter()
    elapsed_ms = (end_time - start_time) * 1000
    for r in results:
        print(f"{r['record']}: Sensitivity={r['Sensitivity']*100.0:.3f}, PPV={r['PPV']*100.0:.3f}")
        sensitivities.append((r['record'], r['Sensitivity']))
        ppvs.append((r['record'], r['PPV']))

    # Legalacsonyabb 5 Sensitivity
    sensitivities.sort(key=lambda x: x[1])  # Növekvő sorrend
    print("\n5 rekord a legalacsonyabb Sensitivity-vel:")
    for rec, se in sensitivities[:5]:
        print(f"{rec}: Sensitivity={se*100.0:.3f}%")

    # Legalacsonyabb 5 PPV
    ppvs.sort(key=lambda x: x[1])  # Növekvő sorrend
    print("\n5 rekord a legalacsonyabb PPV-vel:")
    for rec, ppv in ppvs[:5]:
        print(f"{rec}: PPV={ppv*100.0:.3f}%")

    # Átlagos metrikák
    avg_sensitivity = np.mean([se for _, se in sensitivities])
    avg_ppv = np.mean([ppv for _, ppv in ppvs])

    print(f"Average Sensitivity (100% means no FN): {avg_sensitivity * 100:.3f}%")
    print(f"Average (100% means no FP) PPV: {avg_ppv * 100:.3f}%")
    print(f"Average SPPV: {(avg_ppv + avg_sensitivity) * 50:.3f}%")
    
    print(f"Futási idő: {elapsed_ms:.3f} ms")
    print(f"Futási idő: {elapsed_ms2:.3f} ms")
