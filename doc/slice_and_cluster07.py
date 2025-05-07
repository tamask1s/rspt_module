# Újraimportálás a kernel reset után
import numpy as np
import matplotlib.pyplot as plt
import wfdb
from sklearn.decomposition import FastICA

# 1. Rekord betöltése
def do_ica(signal):
    #fs = record.fs

    # Ellenőrzés, hogy legalább két csatorna van-e
    if signal.shape[1] < 2:
        raise ValueError("Legalább két csatorna szükséges az ICA-hoz.")

    # 2. Csak az első két csatornát használjuk
    X = signal[:, :]

    # 3. ICA alkalmazása
    ica = FastICA(n_components=2, random_state=42)
    S_ = ica.fit_transform(X)  # Rekonstruált forrásjelek (komponensek)
    return S_
    
if __name__ == "__main__":
    record_path = "/media/sf_SharedFolder/QT/mit-bih-arrhythmia-database-1.0.0/100"
    record = wfdb.rdrecord(record_path)
    signal = record.p_signal
    S_ = do_ica(signal)
    # 4. Komponensek megjelenítése
    plt.figure(figsize=(15, 6))
    for i in range(2):
        plt.plot(S_[:, i], label=f"ICA komponens {i+1}", alpha=0.7)
    plt.title("ICA által szétválasztott komponensek")
    plt.xlabel("Mintavételezési index")
    plt.ylabel("Intenzitás (relatív)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
