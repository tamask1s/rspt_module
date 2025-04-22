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

std::vector<unsigned int> detect_peaks(py::array_t<double> ecg_signal_np, double sampling_rate) {
    py::buffer_info buf = ecg_signal_np.request();
    double* ecg_signal = static_cast<double*>(buf.ptr);
    unsigned int len = buf.shape[0];

    // work bufferek
    std::vector<double> peak_signal(len);
    std::vector<double> filt_signal(len);
    std::vector<double> threshold_signal(len);
    std::vector<unsigned int> peak_indexes;

    peak_detector_offline detector(sampling_rate);
    detector.detect(ecg_signal, len, peak_signal.data(), filt_signal.data(), threshold_signal.data(), &peak_indexes);

    for (ssize_t i = 0; i < peak_indexes.size(); ++i)
        peak_indexes[i] /= 2;
    return peak_indexes;
}


// modul exportálás
PYBIND11_MODULE(rspt_module, m) {
    m.doc() = "2D array squaring module (Python lists & numpy arrays)";
    m.def("square_inplace_list", &square_inplace_list, "In-place squaring of a Python 2D list");
    m.def("square_inplace_numpy", &square_inplace_numpy, "In-place squaring of a numpy 2D array");
    m.def("square_list", &square_list, "Squaring of a Python 2D list, returning a new list");
    m.def("square_numpy", &square_numpy, "Squaring of a numpy 2D array, returning a new array");
    m.def("detect_peaks", &detect_peaks, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"));
}
