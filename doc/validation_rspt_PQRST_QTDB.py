import wfdb
import numpy as np
import matplotlib.pyplot as plt
import rspt_module  # saját modulod, legyen elérhető PYTHONPATH-ban vagy a könyvtárban

# --- Beállítások ---
record_path = '/media/sf_SharedFolder/QT/qt-database-1.0.0/sel301'
sampling_rate = 250.0  # QT adatbázis fix mintavételezése

# --- Adat betöltése ---
record = wfdb.rdrecord(record_path)
annotation = wfdb.rdann(record_path, 'pu1')

print("📄 Jel felvételi információk:")
print(f"  Mintavételezési frekvencia: {record.fs} Hz")
print(f"  Csatornák: {record.sig_name}")
print(f"  Jel alakja: {record.p_signal.shape}")
print(f"  Első 10 annotáció: {annotation.aux_note[:10]}")

# --- Jel és annotáció megjelenítése ---
fig, ax = plt.subplots(2, 1, figsize=(15, 6), sharex=True)
time_axis = np.arange(record.p_signal.shape[0]) / sampling_rate

for i in range(min(2, record.p_signal.shape[1])):
    ax[i].plot(time_axis, record.p_signal[:, i], label=f'Elvezetés: {record.sig_name[i]}')
    ax[i].scatter(np.array(annotation.sample) / sampling_rate,
                  record.p_signal[annotation.sample, i],
                  color='red', marker='x', label='Annotációk (pu1)')
    ax[i].legend()
    ax[i].set_title(f'SEL301 - {record.sig_name[i]}')
plt.tight_layout()

# --- Algoritmus futtatása ---
ecg_data = np.array(record.p_signal, dtype=np.float64)
result = rspt_module.analyse_ecg(ecg_data, sampling_rate)

print("\n🧠 Algoritmus eredményei:")
for key, value in result.items():
    print(f"{key}: {value}")

# --- Második ábra: Detektált R-csúcsok vizualizációja ---
fig2, ax2 = plt.subplots(figsize=(15, 4))
ax2.plot(time_axis, ecg_data[:, 0], label='ECG - 1. elvezetés')
r_peak_idxs = [i for i, ann in enumerate(result['annotations']) if 'R' in ann.upper()]
ax2.scatter(np.array(r_peak_idxs) / sampling_rate,
            ecg_data[r_peak_idxs, 0],
            color='green', marker='o', label='Detektált R-csúcsok')
ax2.legend()
ax2.set_title("Detektált R-csúcsok az algoritmus alapján")
plt.tight_layout()

# --- Összehasonlítás a referenciával ---
print("\n📏 Validáció:")
manual_r = [s for s, a in zip(annotation.sample, annotation.aux_note) if '(N' in a or '(AFIB' in a]
if manual_r:
    auto_r = np.array(r_peak_idxs)
    manual_r = np.array(manual_r)
    diff_matrix = np.abs(np.subtract.outer(auto_r, manual_r))
    min_diffs = np.min(diff_matrix, axis=1)
    mean_diff_ms = np.mean(min_diffs) * 1000.0 / sampling_rate
    print(f"Átlagos eltérés annotációkhoz képest: {mean_diff_ms:.2f} ms")
    if mean_diff_ms < 20.0:
        print("✅ Megfelel az elvárásnak (<20 ms)")
    else:
        print("❌ Nem felel meg az elvárásnak (>20 ms)")
else:
    print("⚠️ Nem található összehasonlítható annotáció (pl. R hullám).")

plt.show()