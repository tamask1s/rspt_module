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

static void copy_ecg_signal_to_channels(
    py::array_t<double, py::array::c_style | py::array::forcecast> ecg_signal_np,
    size_t& len,
    size_t& nr_channels,
    std::vector<std::vector<double>>& data,
    std::vector<const double*>& data_ptrs)
{
    py::buffer_info buf = ecg_signal_np.request();
    if (buf.ndim != 1 && buf.ndim != 2)
        throw std::runtime_error("Input ECG array must be 1 or 2-dimensional (samples x channels)");

    len = buf.shape[0];
    nr_channels = (buf.ndim == 1) ? 1 : buf.shape[1];
    data.assign(nr_channels, std::vector<double>(len));
    data_ptrs.resize(nr_channels);

    char* base_ptr = static_cast<char*>(buf.ptr);
    size_t stride0 = buf.strides[0];

    if (buf.ndim == 1)
    {
        for (size_t i = 0; i < len; ++i)
            data[0][i] = *reinterpret_cast<double*>(base_ptr + i * stride0);
        data_ptrs[0] = data[0].data();
        return;
    }

    size_t stride1 = buf.strides[1];
    for (size_t ch = 0; ch < nr_channels; ++ch)
    {
        for (size_t i = 0; i < len; ++i)
            data[ch][i] = *reinterpret_cast<double*>(base_ptr + i * stride0 + ch * stride1);
        data_ptrs[ch] = data[ch].data();
    }
}

static void set_detector_mode(peak_detector_offline& detector, const std::string& mode)
{
    if (mode == "high_sensitivity")
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == "high_ppv")
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);
}

static py::dict annotation_to_dict(const pqrst_indxes& ann)
{
    py::dict complex;
    complex["p"] = py::make_tuple(ann.p[0], ann.p[1], ann.p[2]);
    complex["r"] = py::make_tuple(ann.r[0], ann.r[1], ann.r[2]);
    complex["t"] = py::make_tuple(ann.t[0], ann.t[1], ann.t[2]);
    complex["p_full"] = py::make_tuple(ann.p[0], ann.p[1], ann.p[2], ann.p[3], ann.p[4], ann.p[5]);
    complex["r_full"] = py::make_tuple(ann.r[0], ann.r[1], ann.r[2], ann.r[3], ann.r[4]);
    complex["t_full"] = py::make_tuple(ann.t[0], ann.t[1], ann.t[2], ann.t[3]);
    return complex;
}

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
    set_detector_mode(detector, mode);
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
    set_detector_mode(detector, mode);

    detector.detect_multichannel(data_ptrs.data(), nr_channels, len, peak_signal.data(), filt_signal.data(), threshold_signal.data(), &peak_indexes);

    return peak_indexes;
}

// Python binding az analyse_ecg függvényhez
py::dict analyse_ecg(py::array_t<double, py::array::c_style | py::array::forcecast> ecg_signal_np, double sampling_rate, std::string mode = "default", int analysis_ch_indx = -1, int analysis_peak_indx = 0)
{
    size_t len = 0;
    size_t nr_channels = 0;
    std::vector<std::vector<double>> data;
    std::vector<const double*> data_ptrs;
    copy_ecg_signal_to_channels(ecg_signal_np, len, nr_channels, data, data_ptrs);

    // 3) Analyse eredményének előkészítése
    std::vector<pqrst_indxes> annotations;
    ecg_analysis_result result = analyse_ecg_detect_peaks(data_ptrs.data(), nr_channels, len, sampling_rate, annotations, 0, mode, analysis_ch_indx, analysis_peak_indx);

    py::list py_annotations;
    for (const auto& ann : annotations)
        py_annotations.append(annotation_to_dict(ann));

    // 4) Eredmények Python dict-be töltése
    py::dict output;
    output["annotations"] = py_annotations;
    output["rr_interval_ms"]       = result.rr_interval_ms;
    output["rr_variation_ms"]      = result.rr_variation_ms;
    output["heart_rate_bpm"]       = result.heart_rate_bpm;
    output["pr_interval_ms"]       = result.pr_interval_ms;
    output["pr_segment_ms"]        = result.pr_segment_ms;
    output["qrs_duration_ms"]      = result.qrs_duration_ms;
    output["qt_interval_ms"]       = result.qt_interval_ms;
    output["st_segment_ms"]        = result.st_segment_ms;
    output["p_wave_duration_ms"]   = result.p_wave_duration_ms;
    output["t_wave_duration_ms"]   = result.t_wave_duration_ms;
    output["is_sinus_rhythm"]           = result.is_sinus_rhythm;
    output["premature_beat_count"]      = result.premature_beat_count;
    output["analysis_status"]           = result.analysis_status;
    output["status_message"]            = std::string(result.status_message);

    output["p1_wave_duration_ms"] = result.p1_wave_duration_ms;
    output["p1_amplitude_input_units"] = result.p1_amplitude_input_units;
    output["p2_wave_duration_ms"] = result.p2_wave_duration_ms;
    output["p2_amplitude_input_units"] = result.p2_amplitude_input_units;
    output["q_duration_ms"] = result.q_duration_ms;
    output["q_amplitude_input_units"] = result.q_amplitude_input_units;
    output["r_duration_ms"] = result.r_duration_ms;
    output["r_amplitude_input_units"] = result.r_amplitude_input_units;
    output["s_duration_ms"] = result.s_duration_ms;
    output["s_amplitude_input_units"] = result.s_amplitude_input_units;
    output["j_point_amplitude_input_units"] = result.j_point_amplitude_input_units;
    output["st20_amplitude_input_units"] = result.st20_amplitude_input_units;
    output["st40_amplitude_input_units"] = result.st40_amplitude_input_units;
    output["st60_amplitude_input_units"] = result.st60_amplitude_input_units;
    output["st80_amplitude_input_units"] = result.st80_amplitude_input_units;
    output["t_amplitude_input_units"] = result.t_amplitude_input_units;

    output["P1_DURATION"] = result.p1_wave_duration_ms;
    output["P1_AMPLITUDE"] = result.p1_amplitude_input_units;
    output["P2_DURATION"] = result.p2_wave_duration_ms;
    output["P2_AMPLITUDE"] = result.p2_amplitude_input_units;

    output["PQ_INTERVAL"] = result.pr_interval_ms;

    output["Q_DURATION"] = result.q_duration_ms;
    output["Q_AMPLITUDE"] = result.q_amplitude_input_units;
    output["R_DURATION"] = result.r_duration_ms;
    output["R_AMPLITUDE"] = result.r_amplitude_input_units;
    output["S_DURATION"] = result.s_duration_ms;
    output["S_AMPLITUDE"] = result.s_amplitude_input_units;

    output["QRS_DURATION"] = result.qrs_duration_ms;
    output["QT_INTERVAL"] = result.qt_interval_ms;

    output["J_AMPLITUDE"] = result.j_point_amplitude_input_units;
    output["ST_20_AMPLITUDE"] = result.st20_amplitude_input_units;
    output["ST_40_AMPLITUDE"] = result.st40_amplitude_input_units;
    output["ST_60_AMPLITUDE"] = result.st60_amplitude_input_units;
    output["ST_80_AMPLITUDE"] = result.st80_amplitude_input_units;
    output["T_AMPLITUDE"] = result.t_amplitude_input_units;

    return output;
}

py::dict analyse_ecg_beats(py::array_t<double, py::array::c_style | py::array::forcecast> ecg_signal_np, double sampling_rate, std::string mode = "default", int analysis_ch_indx = -1)
{
    size_t len = 0;
    size_t nr_channels = 0;
    std::vector<std::vector<double>> data;
    std::vector<const double*> data_ptrs;
    copy_ecg_signal_to_channels(ecg_signal_np, len, nr_channels, data, data_ptrs);

    if (analysis_ch_indx < 0)
        analysis_ch_indx = 0;
    if (analysis_ch_indx >= static_cast<int>(nr_channels))
        analysis_ch_indx = static_cast<int>(nr_channels) - 1;

    std::vector<double> peak_signal(len), filt_signal(len), threshold_signal(len);
    std::vector<unsigned int> peak_indexes;
    peak_detector_offline detector(sampling_rate);
    set_detector_mode(detector, mode);
    detector.detect(data[analysis_ch_indx].data(), len, peak_signal.data(), filt_signal.data(), threshold_signal.data(), &peak_indexes);

    std::vector<pqrst_indxes> annotations;
    std::vector<ecg_analysis_result> results;
    analyse_ecg_all_beats(data_ptrs.data(), nr_channels, len, sampling_rate, peak_indexes, annotations, results, analysis_ch_indx);

    py::list py_peak_indexes;
    for (const auto peak_index : peak_indexes)
        py_peak_indexes.append(peak_index);

    py::list py_beats;
    for (size_t i = 0; i < results.size(); ++i)
    {
        const auto& result = results[i];
        py::dict beat;
        beat["beat_index"] = i;
        if (i < peak_indexes.size())
            beat["r_peak_sample"] = peak_indexes[i];
        else
            beat["r_peak_sample"] = py::none();

        beat["analysis_channel_index"] = analysis_ch_indx;
        beat["analysis_status"] = result.analysis_status;
        beat["status_message"] = std::string(result.status_message);

        if (i < annotations.size())
            beat["annotation"] = annotation_to_dict(annotations[i]);
        else
            beat["annotation"] = py::none();

        if (i > 0 && i < peak_indexes.size() && peak_indexes[i] > peak_indexes[i - 1])
        {
            double rr_ms = (peak_indexes[i] - peak_indexes[i - 1]) / sampling_rate * 1000.0;
            beat["rr_interval_ms"] = rr_ms;
            beat["heart_rate_bpm"] = (rr_ms > 0.0) ? 60000.0 / rr_ms : 0.0;
        }
        else
        {
            beat["rr_interval_ms"] = py::none();
            beat["heart_rate_bpm"] = py::none();
        }

        beat["p_wave_duration_ms"] = result.p_wave_duration_ms;
        beat["p1_wave_duration_ms"] = result.p1_wave_duration_ms;
        beat["p1_amplitude_input_units"] = result.p1_amplitude_input_units;
        beat["p2_wave_duration_ms"] = result.p2_wave_duration_ms;
        beat["p2_amplitude_input_units"] = result.p2_amplitude_input_units;
        beat["pr_interval_ms"] = result.pr_interval_ms;
        beat["pr_segment_ms"] = result.pr_segment_ms;
        beat["q_duration_ms"] = result.q_duration_ms;
        beat["q_amplitude_input_units"] = result.q_amplitude_input_units;
        beat["r_duration_ms"] = result.r_duration_ms;
        beat["r_amplitude_input_units"] = result.r_amplitude_input_units;
        beat["s_duration_ms"] = result.s_duration_ms;
        beat["s_amplitude_input_units"] = result.s_amplitude_input_units;
        beat["qrs_duration_ms"] = result.qrs_duration_ms;
        beat["qt_interval_ms"] = result.qt_interval_ms;
        beat["st_segment_ms"] = result.st_segment_ms;
        beat["t_wave_duration_ms"] = result.t_wave_duration_ms;
        beat["j_point_amplitude_input_units"] = result.j_point_amplitude_input_units;
        beat["st20_amplitude_input_units"] = result.st20_amplitude_input_units;
        beat["st40_amplitude_input_units"] = result.st40_amplitude_input_units;
        beat["st60_amplitude_input_units"] = result.st60_amplitude_input_units;
        beat["st80_amplitude_input_units"] = result.st80_amplitude_input_units;
        beat["t_amplitude_input_units"] = result.t_amplitude_input_units;

        beat["P1_DURATION"] = result.p1_wave_duration_ms;
        beat["P1_AMPLITUDE"] = result.p1_amplitude_input_units;
        beat["P2_DURATION"] = result.p2_wave_duration_ms;
        beat["P2_AMPLITUDE"] = result.p2_amplitude_input_units;
        beat["PQ_INTERVAL"] = result.pr_interval_ms;
        beat["Q_DURATION"] = result.q_duration_ms;
        beat["Q_AMPLITUDE"] = result.q_amplitude_input_units;
        beat["R_DURATION"] = result.r_duration_ms;
        beat["R_AMPLITUDE"] = result.r_amplitude_input_units;
        beat["S_DURATION"] = result.s_duration_ms;
        beat["S_AMPLITUDE"] = result.s_amplitude_input_units;
        beat["QRS_DURATION"] = result.qrs_duration_ms;
        beat["QT_INTERVAL"] = result.qt_interval_ms;
        beat["J_AMPLITUDE"] = result.j_point_amplitude_input_units;
        beat["ST_20_AMPLITUDE"] = result.st20_amplitude_input_units;
        beat["ST_40_AMPLITUDE"] = result.st40_amplitude_input_units;
        beat["ST_60_AMPLITUDE"] = result.st60_amplitude_input_units;
        beat["ST_80_AMPLITUDE"] = result.st80_amplitude_input_units;
        beat["T_AMPLITUDE"] = result.t_amplitude_input_units;

        py_beats.append(beat);
    }

    py::dict output;
    output["analysis_channel_index"] = analysis_ch_indx;
    output["r_peak_indexes"] = py_peak_indexes;
    output["beats"] = py_beats;
    return output;
}

PYBIND11_MODULE(rspt_module, m)
{
    m.doc() = "ECG analysis module";
    m.def("detect_peaks", &detect_peaks, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
    m.def("detect_multichannel", &detect_multichannel, "Detect ECG peaks", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default");
    m.def("analyse_ecg", &analyse_ecg, "Perform ECG analysis (multichannel). Returns annotations and parameter dictionary", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default", py::arg("analysis_ch_indx") = -1, py::arg("analysis_peak_indx") = 0);
    m.def("analyse_ecg_beats", &analyse_ecg_beats, "Perform ECG analysis for every detected beat", py::arg("ecg_signal"), py::arg("sampling_rate"), py::arg("mode") = "default", py::arg("analysis_ch_indx") = -1);
}
