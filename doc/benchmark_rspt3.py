import wfdb
from wfdb import processing
import time
import numpy as np
from rspt_module import detect_peaks
from rspt_module import detect_multichannel
from slice_and_cluster07 import do_ica
import neurokit2 as nk

elapsed_ms2 = 0

detection_metod = 1
ICA_is_used = False
multichannel_data = True

def benchmark_record(record_name, tolerance=0.05):
    global elapsed_ms2
    global detection_metod
    global ICA_is_used
    global multichannel_data

    record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
    record_path = "/media/sf_SharedFolder/QT/st-petersburg-incart-12-lead-arrhythmia-database-1.0.0/files/"
    record = wfdb.rdrecord(record_path + record_name)
    annotation = wfdb.rdann(record_path + record_name, 'atr')

    fs = record.fs
    signal = []
    if ICA_is_used:
        signal = do_ica(record.p_signal)
        signal = signal[:, 0]
    else:
        if multichannel_data:
            signal = record.p_signal[:, :]
        else:
            signal = record.p_signal[:, 0]

    detected_peaks = []
    start_time2 = time.perf_counter()

    if detection_metod == 1 or detection_metod == 'RSPT':
        detection_metod = 'RSPT'
        detected_peaks = detect_multichannel(signal, fs, 'default')
        #detected_peaks = detect_peaks(signal, fs, 'default')
    if detection_metod == 2 or detection_metod == 'neurokit ecg_peaks':
        detection_metod = 'neurokit ecg_peaks'
        signals, info = nk.ecg_peaks(signal, sampling_rate=fs, method="neurokit")
        detected_peaks = info["ECG_R_Peaks"]
    if detection_metod == 3 or detection_metod == 'neurokit ecg_process':
        detection_metod = 'neurokit ecg_process'
        signals, info = nk.ecg_process(signal, sampling_rate=fs)
        detected_peaks = info["ECG_R_Peaks"]    
    if detection_metod == 4 or detection_metod == 'wfdb.processing':
        detection_metod = 'wfdb.processing'
        #conf = wfdb.processing.qrs.GQRS.Conf(
        #    fs=record.fs,
        #    adc_gain=record.adcgain[0],
        #    adc_zero=record.adczero[0],
        #    threshold=1.0,
        #    hr=75,
        #    RRdelta=0.2,
        #    RRmin=0.28,
        #    RRmax=2.4,
        #    QS=0.07,
        #    QT=0.35,
        #    RTmin=0.25,
        #    RTmax=0.33,
        #    QRSa=750,
        #    QRSamin=130
        #)

        #gqrs = wfdb.processing.qrs.GQRS()
        #annotations = gqrs.detect(x=signal, conf=conf, adc_zero=record.adczero[0])
        #detected_peaks = wfdb.processing.gqrs_detect(signal, fs=record.fs)

        #detected_peaks = np.array([a.time for a in annotations])
        detected_peaks = processing.gqrs_detect(sig=signal, fs=fs)

    end_time2 = time.perf_counter()
    elapsed_ms2 = elapsed_ms2 + (end_time2 - start_time2) * 1000

    #annot_r_peaks = annotation.sample[np.isin(annotation.symbol, ['N', 'L', 'R', 'V', 'A', 'F', 'j', 'E', 'e', 'a', 'J', 'S'])]
    annot_r_peaks = annotation.sample[np.isin(annotation.symbol, ['N', 'L', 'R', 'B', 'A', 'a', 'J', 'S', 'V', 'r', 'F', 'e', 'j', 'n', 'E', '/', 'f', 'Q', '?'])]

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
    records = ['I{:02d}'.format(i) for i in range(1, 16)] #76
    #records = ['108', '114', '222', '233', '117', '214', '113', '115']
    #records = ['100', '101', '103', '105', '111', '112', '113', '115', '117', '118', '119', '121', '122', '200', '202', '209']    # ide még bővíthetsz
    #records = ['100', '101']
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

    print('records:', records)
    print('Detection metod:', detection_metod, " ICA_is_used:", ICA_is_used, " multichannel_data:", multichannel_data)

    print(f"Average Sensitivity (100% means no FN): {avg_sensitivity * 100:.3f}%")
    print(f"Average (100% means no FP) PPV: {avg_ppv * 100:.3f}%")
    print(f"Average SPPV: {(avg_ppv + avg_sensitivity) * 50:.3f}%")
    
    print(f"Elapsed detection time (with data open): {elapsed_ms:.3f} ms")
    print(f"Elapsed detection time (processing only): {elapsed_ms2:.3f} ms")