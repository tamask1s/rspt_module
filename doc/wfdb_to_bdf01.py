import os
import wfdb
import numpy as np
import pyedflib
from pathlib import Path

# QTDB alapkönyvtár
BASE_DIR = Path('/media/sf_SharedFolder/QT/qt-database-1.0.0/')
OUTPUT_DIR = BASE_DIR / 'converted_edf_bdf'
OUTPUT_DIR.mkdir(exist_ok=True)

def load_qtdb_record(record_path):
    record = wfdb.rdrecord(record_path)
    ann = wfdb.rdann(record_path, 'pu1')
    fs = record.fs

    # A valódi fiziológiai jelek kiszűrése (pl. ha utolsó csatorna "annotation")
    signal = record.p_signal
    labels = record.sig_name

    # Ha van 'ECG annotation' vagy hasonló, akkor az utolsó csatorna valószínűleg az
    #if labels[-1].lower().startswith('ecg') or 'ann' in labels[-1].lower():
    #    signal = signal[:, :-1]
    #    labels = labels[:-1]

    # Annotációs események kódolása
    ann_channel = np.zeros(signal.shape[0])

    # Annotációk bejárása
    for idx, (samp, symb) in enumerate(zip(ann.sample, ann.symbol)):
        code = 0
        if symb in ['(', ')']:
            code = 1  # Onset/Offset
        elif symb.upper() == 'P':
            code = 2
        elif symb.upper() in ['N', 'R']:
            code = 3
        elif symb.upper() == 'T':
            code = 4
        else:
            continue  # ha nem ismert, kihagyjuk

        if 0 <= samp < len(ann_channel):
            ann_channel[samp] = code

    # Jel és annotációs csatorna összefűzése
    signal_with_ann = np.column_stack((signal, ann_channel))
    labels_with_ann = labels + ['annotation']

    return signal_with_ann, ann, fs, labels_with_ann

def save_to_edf_or_bdf(filename, signal, fs, channel_labels, prefer_bdf=True):
    num_signals = signal.shape[1]
    duration = signal.shape[0]

    # Bemeneti adat átalakítás 16 bit integerre, ha szükséges
    if signal.dtype != np.int16:
        scaled_signal = (signal * 1000).astype(np.int16)
    else:
        scaled_signal = signal

    # Metadata
    signal_headers = []
    for label in channel_labels:
        signal_headers.append({
            'label': label,
            'dimension': 'uV',
            'sample_rate': fs,
            'physical_min': -32768,
            'physical_max': 32767,
            'digital_min': -32768,
            'digital_max': 32767,
            'transducer': '',
            'prefilter': ''
        })

    #file_type = pyedflib.FILETYPE_BDFPLUS if prefer_bdf else pyedflib.FILETYPE_EDFPLUS
    file_type = pyedflib.FILETYPE_BDF if prefer_bdf else pyedflib.FILETYPE_EDF

    try:
        with pyedflib.EdfWriter(str(filename), num_signals, file_type=file_type) as f:
            f.setSignalHeaders(signal_headers)
            f.writeSamples([scaled_signal[:, i] for i in range(num_signals)])
        print(f"✔️ Mentve: {filename}")
    except Exception as e:
        if prefer_bdf:
            print(f"⚠️ BDF mentés nem sikerült {filename}, próbálkozás EDF-fel... ({e})")
            save_to_edf_or_bdf(filename.with_suffix('.edf'), signal, fs, channel_labels, prefer_bdf=False)
        else:
            print(f"❌ EDF mentés is sikertelen: {filename} ({e})")

def convert_all_records(base_path):
    records = [f.stem for f in base_path.glob('sel*.dat')]
    print(f"{len(records)} rekord található.")

    for rec_name in records:
        full_path = base_path / rec_name
        try:
            signal, ann, fs, labels = load_qtdb_record(str(full_path))
            out_file = OUTPUT_DIR / f"{rec_name}.bdf"
            save_to_edf_or_bdf(out_file, signal, fs, labels)
        except Exception as e:
            print(f"❌ Hiba történt {rec_name} feldolgozásakor: {e}")

if __name__ == "__main__":
    convert_all_records(BASE_DIR)
