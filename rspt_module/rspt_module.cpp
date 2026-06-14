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

    output["P1_DURATION"] = result.result.P1_DURATION;
    output["P1_AMPLITUDE"] = result.result.P1_AMPLITUDE;
    output["P2_DURATION"] = result.result.P2_DURATION;
    output["P2_AMPLITUDE"] = result.result.P2_AMPLITUDE;

    output["PQ_INTERVAL"] = result.result.PQ_INTERVAL;

    output["Q_DURATION"] = result.result.Q_DURATION;
    output["Q_AMPLITUDE"] = result.result.Q_AMPLITUDE;
    output["R_DURATION"] = result.result.R_DURATION;
    output["R_AMPLITUDE"] = result.result.R_AMPLITUDE;
    output["S_DURATION"] = result.result.S_DURATION;
    output["S_AMPLITUDE"] = result.result.S_AMPLITUDE;

    output["QRS_DURATION"] = result.result.QRS_DURATION;
    output["QT_INTERVAL"] = result.result.QT_INTERVAL;

    output["J_AMPLITUDE"] = result.result.J_AMPLITUDE;
    output["ST_20_AMPLITUDE"] = result.result.ST_20_AMPLITUDE;
    output["ST_40_AMPLITUDE"] = result.result.ST_40_AMPLITUDE;
    output["ST_60_AMPLITUDE"] = result.result.ST_60_AMPLITUDE;
    output["ST_80_AMPLITUDE"] = result.result.ST_80_AMPLITUDE;
    output["T_AMPLITUDE"] = result.result.T_AMPLITUDE;

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
            beat["previous_rr_interval_ms"] = rr_ms;
            beat["rr_interval_ms"] = rr_ms;
            beat["heart_rate_bpm"] = (rr_ms > 0.0) ? 60000.0 / rr_ms : 0.0;
            if (rr_ms > 0.0 && result.result.QT_INTERVAL > 0)
            {
                double qtc_ms = result.result.QT_INTERVAL / std::sqrt(rr_ms / 1000.0);
                beat["qtc_bazett_ms"] = qtc_ms;
                beat["qtc_interval_ms"] = qtc_ms;
            }
            else
            {
                beat["qtc_bazett_ms"] = py::none();
                beat["qtc_interval_ms"] = py::none();
            }
        }
        else
        {
            beat["previous_rr_interval_ms"] = py::none();
            beat["rr_interval_ms"] = py::none();
            beat["heart_rate_bpm"] = py::none();
            beat["qtc_bazett_ms"] = py::none();
            beat["qtc_interval_ms"] = py::none();
        }

        if (i + 1 < peak_indexes.size() && peak_indexes[i + 1] > peak_indexes[i])
            beat["next_rr_interval_ms"] = (peak_indexes[i + 1] - peak_indexes[i]) / sampling_rate * 1000.0;
        else
            beat["next_rr_interval_ms"] = py::none();

        if (i > 0 && i < annotations.size() && i - 1 < annotations.size() && annotations[i].p[1] >= 0 && annotations[i - 1].p[1] >= 0)
            beat["pp_interval_ms"] = (annotations[i].p[1] - annotations[i - 1].p[1]) / sampling_rate * 1000.0;
        else
            beat["pp_interval_ms"] = py::none();

        beat["p_wave_duration_ms"] = result.p_wave_duration_ms;
        beat["pr_interval_ms"] = result.pr_interval_ms;
        beat["pr_segment_ms"] = result.pr_segment_ms;
        beat["qrs_duration_ms"] = result.qrs_duration_ms;
        beat["qt_interval_ms"] = result.qt_interval_ms;
        beat["st_segment_ms"] = result.st_segment_ms;
        beat["t_wave_duration_ms"] = result.t_wave_duration_ms;

        beat["P1_DURATION"] = result.result.P1_DURATION;
        beat["P1_AMPLITUDE"] = result.result.P1_AMPLITUDE;
        beat["P2_DURATION"] = result.result.P2_DURATION;
        beat["P2_AMPLITUDE"] = result.result.P2_AMPLITUDE;
        beat["PQ_INTERVAL"] = result.result.PQ_INTERVAL;
        beat["Q_DURATION"] = result.result.Q_DURATION;
        beat["Q_AMPLITUDE"] = result.result.Q_AMPLITUDE;
        beat["R_DURATION"] = result.result.R_DURATION;
        beat["R_AMPLITUDE"] = result.result.R_AMPLITUDE;
        beat["S_DURATION"] = result.result.S_DURATION;
        beat["S_AMPLITUDE"] = result.result.S_AMPLITUDE;
        beat["QRS_DURATION"] = result.result.QRS_DURATION;
        beat["QT_INTERVAL"] = result.result.QT_INTERVAL;
        beat["J_AMPLITUDE"] = result.result.J_AMPLITUDE;
        beat["ST_20_AMPLITUDE"] = result.result.ST_20_AMPLITUDE;
        beat["ST_40_AMPLITUDE"] = result.result.ST_40_AMPLITUDE;
        beat["ST_60_AMPLITUDE"] = result.result.ST_60_AMPLITUDE;
        beat["ST_80_AMPLITUDE"] = result.result.ST_80_AMPLITUDE;
        beat["T_AMPLITUDE"] = result.result.T_AMPLITUDE;

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
