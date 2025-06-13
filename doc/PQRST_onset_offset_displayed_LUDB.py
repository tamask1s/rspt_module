import wfdb
import numpy as np
import matplotlib.pyplot as plt

# Beállítások
data_dir    = '/media/sf_SharedFolder/QT/lobachevsky-university-electrocardiography-database-1.0.1/data'
record_name = '1'
ann_ext     = 'i'
ch_idx      = 1

# Jel és annotáció beolvasása
rec_path = f"{data_dir}/{record_name}"
signals, fields = wfdb.rdsamp(rec_path)
fs       = fields['fs']
sig_name = fields['sig_name'][ch_idx]
ann      = wfdb.rdann(rec_path, ann_ext)

# Időtengely
t = np.arange(signals.shape[0]) / fs

# Szimbólum lista és színek (fentről lefelé)
symbols = ['(', 'p', ')', '[', 'N', ']', '{', 't', '}']
colors  = ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen',
           'lightblue', 'darkblue', 'purple', 'brown']
labels  = [
    'P-hullám on', 'P-hullám csúcs', 'P-hullám off',
    'QRS on', 'QRS csúcs', 'QRS off',
    'T-hullám on', 'T-hullám csúcs', 'T-hullám off'
]

# Ábra
plt.figure(figsize=(12,4))
plt.plot(t, signals[:,ch_idx], color='black', label=sig_name)

# Scatter csak pöttyökkel, szín szerint
for sym, col, lab in zip(symbols, colors, labels):
    idxs = np.where(np.array(ann.symbol) == sym)[0]
    if idxs.size > 0:
        plt.scatter(
            ann.sample[idxs] / fs,
            signals[ann.sample[idxs], ch_idx],
            c=col,
            marker='o',
            s=50,
            label=lab
        )

plt.xlabel('Idő [s]')
plt.ylabel('Amplitúdó [mV]')
plt.title(f"LUDB rekord {record_name} – csatorna {sig_name} – annotáció ({ann_ext})")
plt.legend(loc='upper right', ncol=2)
plt.tight_layout()
plt.show()
