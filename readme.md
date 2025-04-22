# rspt_module

Python module for rspt's peak detection algorithm, tested with a few ECG signals.

![ECG peak detection result](doc/result_mitdb_200.png)

## Install

```bash
pip install git+https://github.com/tamask1s/rspt_module.git
```

[doc/rspt_test.py](doc/rspt_test.py)

```python
import numpy as np
import wfdb
import rspt_module
import matplotlib.pyplot as plt

print("opening")
#record = wfdb.rdrecord('100', pn_dir='mitdb')
#record = wfdb.rdrecord('chf07', pn_dir='chfdb', channels=[0], sampto=650000)
#record = wfdb.rdrecord('e0103', pn_dir='edb', channels=[0], sampto=650000)
record = wfdb.rdrecord('200', pn_dir='mitdb')
print("opened")
ecg_signal = record.p_signal[:, 0]
sampling_rate = record.fs

print("Detecting. sampling_rate:", sampling_rate)
peak_indexes = rspt_module.detect_peaks(ecg_signal, sampling_rate)
print(f"{len(peak_indexes)} peaks detected.")
print("First 10 peak indexes:", peak_indexes[:10])

plt.figure(figsize=(15, 5))
plt.plot(ecg_signal, label='ECG')
plt.plot(peak_indexes, ecg_signal[peak_indexes], 'ro', label='Peaks')
plt.xlabel('Sample index')
plt.ylabel('Amplitude (mV)')
plt.title('ECG signal and detected peaks')
plt.legend()
plt.grid(True)
plt.show()

 ```
