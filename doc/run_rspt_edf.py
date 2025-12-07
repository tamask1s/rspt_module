import pyedflib
import numpy as np
import matplotlib.pyplot as plt
import rspt_module

RECORD_PATH = '/media/sf_SharedFolder/CSE/CS50132D_CTS_Database_delivery/EDF_Output/ANE20000.edf'
RECORD_PATH = '/media/sf_SharedFolder/CSE/CS50132D_CTS_Database_delivery/EDF_Output/CAL50000.edf'

def load_edf_record(record_path):
    with pyedflib.EdfReader(record_path) as f:
        n_channels = f.signals_in_file
        fs = f.getSampleFrequency(0)
        signal_length = f.getNSamples()[0]
        
        signals = np.zeros((signal_length, n_channels))
        for i in range(n_channels):
            signals[:, i] = f.readSignal(i)
    
    return signals, fs

def plot_rspt_annotations(signal, rspt_annotations, fs):
    ann = rspt_annotations[0]
    p_on, p_peak, p_off = ann['p']
    r_on, r_peak, r_off = ann['r']
    t_on, t_peak, t_off = ann['t']
    
    pad = int(0.2 * fs)
    start = max(0, min(p_on, r_on, t_on) - pad)
    end = min(len(signal), max(p_off, r_off, t_off) + pad)
    
    t = np.arange(start, end) / fs
    sig_segment = signal[start:end]
    
    plt.figure(figsize=(12, 6))
    plt.plot(t, sig_segment, 'k-', linewidth=1, label='ECG')
    
    colors = {'p': 'g', 'r': 'b', 't': 'c'}
    markers = ['o', 'x', '^']
    
    for wave_type in ['p', 'r', 't']:
        for i, sample_idx in enumerate(ann[wave_type]):
            if start <= sample_idx <= end:
                plt.plot(sample_idx / fs, signal[sample_idx], 
                        colors[wave_type] + markers[i], 
                        label=f'{wave_type}_{["on","peak","off"][i]}')
    
    plt.legend()
    plt.title("RSPT annotációk")
    plt.xlabel("Idő (s)")
    plt.ylabel("Amplitúdó")
    plt.grid(True)
    plt.show()

def print_result(result):
    for k, v in result.items():
        if k != "status_message":
            print(f"{k}: {v}")
    print("Státusz:", result["status_message"])

if __name__ == "__main__":
    signal, fs = load_edf_record(RECORD_PATH)
    
    rspt_result = rspt_module.analyse_ecg(signal, fs)
    
    print("--- RSPT eredmények ---")
    print_result(rspt_result)
    
    plot_rspt_annotations(signal[:, 0], rspt_result['annotations'], fs)
