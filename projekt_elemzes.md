# RSPT Module - Projekt Elemzés

## Projekt Áttekintés

Az **rspt_module** egy Python könyvtár, amely ECG (elektrokardiogram) jelek csúcsdetektálására specializálódott. A projekt egy hibrid megközelítést alkalmaz: a számításigényes algoritmusokat C++ nyelven implementálja, majd Python binding-on keresztül teszi elérhetővé.

## Projekt Jellege

**Típus:** Tudományos/orvosi szoftver könyvtár  
**Fő cél:** ECG jelek R-csúcsainak automatikus detektálása  
**Architektúra:** C++ backend + Python frontend  
**Licenc:** Apache License 2.0  

## Funkcionalitás

### Fő Funkciók

1. **ECG Csúcsdetektálás (`detect_peaks`)**
   - Egycsatornás ECG jelek feldolgozása
   - Adaptív küszöbérték alapú detektálás
   - Három működési mód:
     - `default`: Alapértelmezett érzékenység
     - `high_sensitivity`: Magas érzékenység
     - `high_ppv`: Magas pozitív prediktív érték

2. **Többcsatornás Detektálás (`detect_multichannel`)**
   - Több ECG csatorna egyidejű feldolgozása
   - Csatornák közötti átlagolás támogatása

3. **ECG Analízis**
   - Teljes ECG elemzés PQRST hullámokkal
   - Szívritmus paraméterek kiszámítása
   - Patológiás állapotok detektálása

### Algoritmus Jellemzők

- **Sávszűrés:** IIR szűrők alkalmazása
- **Integráció:** Rövid távú fluktuációk simítása
- **Adaptív küszöbérték:** Dinamikus érzékenység beállítás
- **Késleltetés:** 220 ms detektálási késleltetés
- **Valós idejű feldolgozás:** Streaming adatok támogatása

## Technikai Implementáció

### Fájlstruktúra

```
rspt_module/
├── rspt_module/           # Fő modul könyvtár
│   ├── rspt_module.cpp    # Python binding (pybind11)
│   ├── ecg_analysis.cpp   # ECG elemzési algoritmusok
│   ├── ecg_analysis.h     # ECG struktúrák és interfészek
│   ├── peak_detector.h    # Csúcsdetektáló osztály
│   ├── filter.h           # Szűrő implementációk
│   └── lib_filter/        # Szűrő könyvtár
├── doc/                   # Dokumentáció és tesztek
├── test/                  # Unit tesztek
└── setup.py              # Python telepítő
```

### Függőségek

- **Python:** NumPy, matplotlib, wfdb
- **C++:** pybind11 (Python binding)
- **Build:** setuptools, C++ compiler

### Adatstruktúrák

```cpp
struct ecg_analysis_result {
    double heart_rate_bpm;
    double rr_interval_ms;
    double pr_interval_ms;
    double qrs_duration_ms;
    double qt_interval_ms;
    // ... további paraméterek
};
```

## Használat

### Telepítés

```bash
pip install git+https://github.com/tamask1s/rspt_module.git
```

### Alapvető Használat

```python
import numpy as np
import wfdb
import rspt_module
import matplotlib.pyplot as plt

# ECG rekord betöltése
record = wfdb.rdrecord('200', pn_dir='mitdb')
ecg_signal = record.p_signal[:, 0]
sampling_rate = record.fs

# Csúcsdetektálás
peak_indexes = rspt_module.detect_peaks(ecg_signal, sampling_rate)

# Vizualizáció
plt.figure(figsize=(15, 5))
plt.plot(ecg_signal, label='ECG')
plt.plot(peak_indexes, ecg_signal[peak_indexes], 'ro', label='Peaks')
plt.show()
```

### Speciális Módok

```python
# Magas érzékenység
peaks_sensitive = rspt_module.detect_peaks(ecg_signal, fs, "high_sensitivity")

# Magas pontosság
peaks_precise = rspt_module.detect_peaks(ecg_signal, fs, "high_ppv")

# Többcsatornás detektálás
peaks_multi = rspt_module.detect_multichannel(multichannel_ecg, fs)
```

## Tesztelés és Validáció

### Benchmark Adatbázisok

- **MIT-BIH Arrhythmia Database:** Szívritmuszavar adatok
- **CHF Database:** Szívelégtelenség adatok  
- **European ST-T Database:** ST-T változások
- **QT Database:** QT intervallum adatok

### Teljesítmény Metrikák

- **Érzékenység (Sensitivity):** Valós csúcsok detektálási aránya
- **Pozitív Prediktív Érték (PPV):** Detektált csúcsok pontossága
- **F1-score:** Harmonikus átlag Se és PPV között

### Validációs Folyamat

```python
def benchmark_record(record_name, tolerance=0.05):
    # Rekord és annotáció betöltése
    record = wfdb.rdrecord(record_name, pn_dir='mitdb')
    annotation = wfdb.rdann(record_name, 'atr', pn_dir='mitdb')
    
    # Detektálás és összehasonlítás
    detected_peaks = detect_peaks(signal, fs)
    # ... teljesítmény számítás
```

## Fejlesztési Környezet

### Build Rendszer

- **Compiler:** GCC/Clang C++11 támogatással
- **Python Build:** setuptools + pybind11
- **IDE támogatás:** Code::Blocks projekt fájlok

### Debug és Profiling

- Részletes logging rendszer
- Teljesítmény benchmark szkriptek
- Memória optimalizált implementáció

## Alkalmazási Területek

1. **Orvosi Diagnosztika**
   - Szívritmus zavarok detektálása
   - ECG automatikus elemzése
   - Telemedicina alkalmazások

2. **Kutatás és Fejlesztés**
   - Algoritmus validáció
   - Új detektálási módszerek tesztelése
   - Adatbázis annotálás

3. **Valós Idejű Monitoring**
   - Betegmegfigyelő rendszerek
   - Wearable eszközök
   - Távoli monitorozás

## Projekt Státusz

- **Aktív fejlesztés:** Folyamatos commit aktivitás
- **Stabil API:** Jól definiált interfészek
- **Dokumentáció:** Példakódok és benchmark eredmények
- **Tesztelés:** Validált orvosi adatbázisokon

## Összefoglalás

Az rspt_module egy professzionális szintű ECG feldolgozó könyvtár, amely ötvözi a C++ teljesítményét a Python egyszerűségével. A projekt kiváló példája a tudományos szoftver fejlesztésnek, ahol a pontosság és a teljesítmény egyaránt kritikus fontosságú.
