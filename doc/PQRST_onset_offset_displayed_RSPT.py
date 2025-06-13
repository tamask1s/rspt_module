import wfdb
import numpy as np
import matplotlib.pyplot as plt
from rspt_module import analyse_ecg

# Beállítások
data_dir    = '/media/sf_SharedFolder/QT/lobachevsky-university-electrocardiography-database-1.0.1/data'
             #'/media/sf_SharedFolder/QT/qt-database-1.0.0/sel301.hea'
record_name = '1'
ch_idx      = 1  # A megjelenítendő csatorna indexe (pl. I. elvezetés)

# Jel beolvasása
rec_path = f"{data_dir}/{record_name}"
signals, fields = wfdb.rdsamp(rec_path)
fs       = fields['fs']
sig_name = fields['sig_name'][ch_idx]

# Időtengely
t = np.arange(signals.shape[0]) / fs

# 2D tömb Python-nak: (minták x csatornák)
ecg_array = signals   # shape = (N, 12)

# analyse_ecg meghívása (multichannel bemenet, részletes annotációkkal)
result = analyse_ecg(ecg_array, fs, mode="default")

# Lekérjük az annotációkat: formátumuk "index:symbol"
annots = result["annotations"]

print(result)
print("-------------------------")
print("-------------------------")
print("-------------------------")
print(annots)

# Szimbólum lista és színek (Lobachevsky-stílus):
symbols = ['(', 'p', ')', '[', 'N', ']', '{', 't', '}']
colors  = ['red', 'orange', 'yellow', 'lightgreen', 'darkgreen',
           'lightblue', 'darkblue', 'purple', 'brown']
labels  = [
    'P-hullám eleje', 'P-hullám csúcs', 'P-hullám vége',
    'QRS eleje', 'QRS (R) csúcs', 'QRS vége',
    'T-hullám eleje', 'T-hullám csúcs', 'T-hullám vége'
]

# Előkészítjük a szimbólum→(szín, címke) map-et
symbol_to_color = {sym: col for sym, col in zip(symbols, colors)}
symbol_to_label = {sym: lab for sym, lab in zip(symbols, labels)}

# Ábra
plt.figure(figsize=(12, 4))
plt.plot(t, signals[:, ch_idx], color='black', label=sig_name)

# Scatter az analyse_ecg által adott annotációkkal
# Parse: minden annotáció "idx:symbol"
for annot in annots:
    try:
        idx_str, sym = annot.split(':', 1)
        idx = int(idx_str)
    except Exception:
        continue  # ha rossz formátum, kihagyjuk

    # Csak akkor rajzoljuk, ha a szimbólum listában szerepel
    if sym in symbol_to_color:
        time = idx / fs
        amp  = signals[idx, ch_idx]
        plt.scatter(
            time,
            amp,
            c=symbol_to_color[sym],
            marker='o',
            s=50,
            label=symbol_to_label[sym]
        )

# Elkerüljük a duplikált legend-bejegyzéseket
handles, lbls = plt.gca().get_legend_handles_labels()
by_label = dict(zip(lbls, handles))
plt.legend(by_label.values(), by_label.keys(), loc='upper right', ncol=2)

plt.xlabel('Idő [s]')
plt.ylabel('Amplitúdó [mV]')
plt.title(f"LUDB rekord {record_name} – csatorna {sig_name} – analyse_ecg annotációk")
plt.tight_layout()
plt.show()
