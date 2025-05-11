#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>  // hogy tudjunk std::vector-t és list-eket kezelni

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

namespace py = pybind11;

// Python 2D lista inplace négyzetre emelés
void square_inplace_list(py::list array) {
    ssize_t rows = py::len(array);
    for (ssize_t i = 0; i < rows; ++i) {
        py::list row = array[i].cast<py::list>();
        ssize_t cols = py::len(row);
        for (ssize_t j = 0; j < cols; ++j) {
            auto val = py::cast<double>(row[j]);
            row[j] = val * val;
        }
        array[i] = row;
    }
}

// numpy 2D array inplace négyzetre emelés
void square_inplace_numpy(py::array_t<double> array) {
    auto buf = array.mutable_unchecked<2>();  // 2D, gyors eléréshez
    for (ssize_t i = 0; i < buf.shape(0); ++i)
        for (ssize_t j = 0; j < buf.shape(1); ++j)
            buf(i, j) = buf(i, j) * buf(i, j);
}

// Python 2D lista négyzetre emelés, visszaadott új listával
py::list square_list(py::list array) {
    py::list result;
    ssize_t rows = py::len(array);
    for (ssize_t i = 0; i < rows; ++i) {
        py::list row = array[i].cast<py::list>();
        py::list new_row;
        ssize_t cols = py::len(row);
        for (ssize_t j = 0; j < cols; ++j) {
            auto val = py::cast<double>(row[j]);
            new_row.append(val * val);
        }
        result.append(new_row);
    }
    return result;
}

// numpy 2D array négyzetre emelés, visszaadott új tömbbel
py::array_t<double> square_numpy(py::array_t<double> array) {
    auto buf = array.unchecked<2>();
    auto result = py::array_t<double>({buf.shape(0), buf.shape(1)});
    auto res_buf = result.mutable_unchecked<2>();
    for (ssize_t i = 0; i < buf.shape(0); ++i)
        for (ssize_t j = 0; j < buf.shape(1); ++j)
            res_buf(i, j) = buf(i, j) * buf(i, j);
    return result;
}

std::vector<unsigned int> detect_peaks(py::array_t<double> ecg_signal_np, double sampling_rate, std::string mode = "default") {
    py::buffer_info buf = ecg_signal_np.request();
    unsigned int len = buf.shape[0];
    std::vector<double> ecg_signal(len);
    char* base_ptr = static_cast<char*>(buf.ptr);

    if (buf.ndim == 1)
        for (unsigned int i = 0; i < len; ++i)
            ecg_signal[i] = *reinterpret_cast<double*>(base_ptr + i * buf.strides[0]);
    else if (buf.ndim == 2)
    {
        size_t nr_channels = buf.shape[1];
        for (size_t i = 0; i < len; ++i)
        {
            ecg_signal[i] = *(reinterpret_cast<double*>(base_ptr + i * buf.strides[0]));
            for (size_t ch = 1; ch < nr_channels; ++ch)
                ecg_signal[i] += *(reinterpret_cast<double*>(base_ptr + i * buf.strides[0] + ch * buf.strides[1]));
            ecg_signal[i] /= (double)nr_channels;
        }
    }
    else
        cout << "Input ECG array must be 1 or 2-dimensional (samples x channels)" << endl;

    std::vector<double> peak_signal(len), filt_signal(len), threshold_signal(len);
    std::vector<unsigned int> peak_indexes;

    peak_detector_offline detector(sampling_rate);
    if (mode == "high_sensitivity")
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == "high_ppv")
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);
    detector.detect(ecg_signal.data(), len, peak_signal.data(), filt_signal.data(), threshold_signal.data(), &peak_indexes);

    return peak_indexes;
}

std::vector<unsigned int> detect_multichannel(py::array_t<double> ecg_signal_np, double sampling_rate, std::string mode = "default") {

    // Buffer lekérése
    py::buffer_info buf = ecg_signal_np.request();
    if (buf.ndim != 2)
        throw std::runtime_error("Input ECG array must be 2-dimensional (samples x channels)");

    size_t len         = buf.shape[0];   // minta darabszám
    size_t nr_channels = buf.shape[1];   // csatornák száma

    // Adattároló C++-ban
    std::vector<std::vector<double>> data(nr_channels, std::vector<double>(len));
    std::vector<double*> data_ptrs(nr_channels);

    char* base_ptr = static_cast<char*>(buf.ptr);
    size_t stride0 = buf.strides[0];  // lépés egy minta között (byte)
    size_t stride1 = buf.strides[1];  // lépés egy csatorna között (byte)

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

PYBIND11_MODULE(rspt_module, m) {
    m.doc() = "2D array squaring module (Python lists & numpy arrays)";
    m.def("square_inplace_list", &square_inplace_list, "In-place squaring of a Python 2D list");
    m.def("square_inplace_numpy", &square_inplace_numpy, "In-place squaring of a numpy 2D array");
    m.def("square_list", &square_list, "Squaring of a Python 2D list, returning a new list");
    m.def("square_numpy", &square_numpy, "Squaring of a numpy 2D array, returning a new array");
    m.def("detect_peaks", &detect_peaks, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
    m.def("detect_multichannel", &detect_multichannel, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
}
