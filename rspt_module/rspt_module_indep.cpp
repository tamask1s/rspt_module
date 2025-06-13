#include <vector>
#include <cstdio>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <map>
#include <vector>
#include <locale>
#include <chrono>
#include <deque>
#include <cmath>
#include <algorithm>

using namespace std;

#include "filter.h"
#include "iir_filter_opt.h"
#include "peak_detector.h"

std::vector<unsigned int> detect_multichannel(const double** ecg_signal_np, double sampling_rate, std::string mode = "default")
{
    size_t len         = 10000;   // minta darabszám
    size_t nr_channels = 12;   // csatornák száma

    // Adattároló C++-ban
    std::vector<std::vector<double>> data(nr_channels, std::vector<double>(len));
    std::vector<const double*> data_ptrs(nr_channels);

    char* base_ptr = (char*)(&ecg_signal_np[0][0]);
    size_t stride0 = sizeof(double); //buf.strides[0];  // lépés egy minta között (byte)
    size_t stride1 = sizeof(double); //buf.strides[1];  // lépés egy csatorna között (byte)

    // Feltöltjük data[ch][i] = array[i, ch]
    for (size_t ch = 0; ch < nr_channels; ++ch) {
        for (size_t i = 0; i < len; ++i) {
            double* val_ptr = reinterpret_cast<double*>(base_ptr + i * stride0 + ch * stride1);
            data[ch][i] = *val_ptr;
        }
        data_ptrs[ch] = data[ch].data();
    }

    //cout << "nr_channels: " << nr_channels << " len: " << len << endl;
    std::vector<double> peak_signal(len), filt_signal(len), threshold_signal(len);
    std::vector<unsigned int> peak_indexes;

    peak_detector_offline detector(sampling_rate);
    if (mode == "high_sensitivity")
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == "high_ppv")
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);

    detector.detect_multichannel(data_ptrs.data(), nr_channels, len, peak_signal.data(), filt_signal.data(), threshold_signal.data(), &peak_indexes);

    return peak_indexes;
}

int main()
{
    return 0;
}
