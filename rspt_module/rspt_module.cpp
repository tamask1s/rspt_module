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
#include "ecg_analysis.h"

namespace py = pybind11;

std::vector<unsigned int> detect_peaks(py::array_t<double> ecg_signal_np, double sampling_rate, std::string mode = "default")
{
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

std::vector<unsigned int> detect_multichannel(py::array_t<double> ecg_signal_np, double sampling_rate, std::string mode = "default")
{
    // Buffer lekérése
    py::buffer_info buf = ecg_signal_np.request();
    if (buf.ndim != 2)
        throw std::runtime_error("Input ECG array must be 2-dimensional (samples x channels)");

    size_t len         = buf.shape[0];   // minta darabszám
    size_t nr_channels = buf.shape[1];   // csatornák száma

    // Adattároló C++-ban
    std::vector<std::vector<double>> data(nr_channels, std::vector<double>(len));
    std::vector<const double*> data_ptrs(nr_channels);

    char* base_ptr = static_cast<char*>(buf.ptr);
    size_t stride0 = buf.strides[0];  // lépés egy minta között (byte)
    size_t stride1 = buf.strides[1];  // lépés egy csatorna között (byte)

    // Feltöltjük data[ch][i] = array[i, ch]
    for (size_t ch = 0; ch < nr_channels; ++ch)
    {
        for (size_t i = 0; i < len; ++i)
        {
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

// Python binding az analyse_ecg függvényhez
py::dict analyse_ecg(py::array_t<double, py::array::c_style | py::array::forcecast> ecg_signal_np, double sampling_rate, std::string mode = "default")
{
    // 1) Buffer ellenőrzése és adatmozgatás
    py::buffer_info buf = ecg_signal_np.request();
    if (buf.ndim != 2)
        throw std::runtime_error("Input ECG array must be 2-dimensional (samples x channels)");

    size_t len         = buf.shape[0];
    size_t nr_channels = buf.shape[1];

    // 2D C++ tömb létrehozása
    std::vector<std::vector<double>> data(nr_channels, std::vector<double>(len));
    std::vector<const double*> data_ptrs(nr_channels);
    char* base_ptr = static_cast<char*>(buf.ptr);
    size_t stride0 = buf.strides[0];
    size_t stride1 = buf.strides[1];

    for (size_t ch = 0; ch < nr_channels; ++ch)
    {
        for (size_t i = 0; i < len; ++i)
        {
            double* val_ptr = reinterpret_cast<double*>(base_ptr + i * stride0 + ch * stride1);
            data[ch][i] = *val_ptr;
        }
        data_ptrs[ch] = data[ch].data();
    }

    // 2) R-csúcsok detektálása többsávos bemenet alapján
    std::vector<double> peak_signal(len), filt_signal(len), threshold_signal(len);
    std::vector<unsigned int> peak_indexes;
    peak_detector_offline detector(sampling_rate);
    if (mode == "high_sensitivity")
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == "high_ppv")
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);

    detector.detect_multichannel(data_ptrs.data(), nr_channels, (unsigned int)len,
                                 peak_signal.data(), filt_signal.data(), threshold_signal.data(),
                                 &peak_indexes);

    // 3) Analyse eredményének előkészítése
    ecg_analysis_result result;
    std::vector<std::string> annotations;
    analyse_ecg_multichannel(data_ptrs.data(), (unsigned int)nr_channels, (unsigned int)len, sampling_rate, peak_indexes, annotations, result);

    // 4) Eredmények Python dict-be töltése
    py::dict output;
    output["annotations"] = annotations;
    output["rr_interval_ms"]       = result.rr_interval_ms;
    output["rr_variation_ms"]      = result.rr_variation_ms;
    output["heart_rate_bpm"]       = result.heart_rate_bpm;
    output["pr_interval_ms"]       = result.pr_interval_ms;
    output["qrs_duration_ms"]      = result.qrs_duration_ms;
    output["qt_interval_ms"]       = result.qt_interval_ms;
    output["qtc_interval_ms"]      = result.qtc_interval_ms;
    output["p_wave_duration_ms"]   = result.p_wave_duration_ms;
    output["t_wave_duration_ms"]   = result.t_wave_duration_ms;

    // Amplitúdók, 12 elvezetés sorrendjében
    py::list r_peaks, s_waves, st_elev, st_depr;
    for (int i = 0; i < 12; ++i)
    {
        r_peaks.append(result.r_peak_amplitude_mV[i]);
        s_waves.append(result.s_wave_amplitude_mV[i]);
        st_elev.append(result.st_elevation_mV[i]);
        st_depr.append(result.st_depression_mV[i]);
    }
    output["r_peak_amplitude_mV"]   = r_peaks;
    output["s_wave_amplitude_mV"]   = s_waves;
    output["st_elevation_mV"]       = st_elev;
    output["st_depression_mV"]      = st_depr;

    output["frontal_plane_axis_deg"]    = result.frontal_plane_axis_deg;
    output["horizontal_plane_axis_deg"] = result.horizontal_plane_axis_deg;
    output["is_sinus_rhythm"]           = result.is_sinus_rhythm;
    output["premature_beat_count"]      = result.premature_beat_count;
    output["analysis_status"]           = result.analysis_status;
    output["status_message"]            = std::string(result.status_message);

    return output;
}

PYBIND11_MODULE(rspt_module, m)
{
    m.doc() = "ECG analysis module";
    m.def("detect_peaks", &detect_peaks, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
    m.def("detect_multichannel", &detect_multichannel, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
    m.def("analyse_ecg", &analyse_ecg, "Perform ECG analysis (multichannel). Returns annotations and parameter dictionary", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
}
