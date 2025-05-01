import wfdb
import rspt_module
import numpy as np

# Segédfüggvények a metrikák számításához
def calculate_metrics(detected, annotations):
    # A metrikák kiszámítása
    tp = len(set(detected).intersection(set(annotations)))  # True Positives
    fp = len(detected) - tp  # False Positives
    fn = len(annotations) - tp  # False Negatives

    # Sensitivity és PPV számítás
    sensitivity = tp / (tp + fn) if tp + fn > 0 else 0
    ppv = tp / (tp + fp) if tp + fp > 0 else 0

    return sensitivity, ppv

# Benchmark teszt
def run_benchmark(records):
    sensitivities = []
    ppvs = []

    for record in records:
        record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/"
        record_path = record_path + record

        # Rekord beolvasása
        print(f"Processing record {record}")
        try:
            # Ha nincs annotáció, ugrik a rekordra
#            ann = wfdb.rdann(record, 'atr', pn_dir='mitdb')
            ann = wfdb.rdann(record_path, 'atr')
        except Exception as e:
            print(f"Error loading annotation for {record}: {e}")
            continue

        # Jelek beolvasása
#        ecg_signal = wfdb.rdrecord(record, pn_dir='mitdb').p_signal[:, 0]
        ecgsig = wfdb.rdrecord(record_path)
        ecg_signal = ecgsig.p_signal[:, 0]
#        sampling_rate = wfdb.rdrecord(record, pn_dir='mitdb').fs
        sampling_rate = ecgsig.fs

        # Detektálás
        peak_indexes = rspt_module.detect_peaks(ecg_signal, sampling_rate)

        # Kiszámoljuk a metrikákat
        sensitivity, ppv = calculate_metrics(peak_indexes, ann.sample)

        sensitivities.append(sensitivity)
        ppvs.append(ppv)

        # Kép kiírása a rekordonkénti detektálásról (opcionális)
        # (pl. matplotlib ábrázolás)

    # Átlagos metrikák
    avg_sensitivity = np.mean(sensitivities)
    avg_ppv = np.mean(ppvs)

    print(f"Average Sensitivity: {avg_sensitivity * 100:.2f}%")
    print(f"Average PPV: {avg_ppv * 100:.2f}%")

# Benchmark rekordok listája (már az előzőek alapján)
# records = ['100', '101', '103', '104', '107', '109', '113', '118', '119', '121', '122']
records = ['100', '101', '103', '105', '111', '112', '113', '115', '117', '118', '119', '121', '122', '200', '202', '209']

# Benchmark futtatása
run_benchmark(records)
