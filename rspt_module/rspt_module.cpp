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
#include <limits>

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

static int mode_to_c_api_value(const std::string& mode)
{
    if (mode == "high_sensitivity")
        return RSPT_MODE_HIGH_SENSITIVITY;
    if (mode == "high_ppv")
        return RSPT_MODE_HIGH_PPV;
    return RSPT_MODE_DEFAULT;
}

static double numeric_or_zero(double value)
{
    return std::isfinite(value) ? value : 0.0;
}

static py::object numeric_or_none(double value)
{
    if (std::isfinite(value))
        return py::float_(value);
    return py::none();
}

static py::dict annotation_to_dict(const rspt_pqrst_annotation& ann)
{
    py::dict complex;
    complex["p"] = py::make_tuple(ann.p1_onset_sample, ann.p1_peak_sample, ann.p1_offset_sample);
    complex["r"] = py::make_tuple(ann.qrs_onset_sample, ann.r_peak_sample, ann.qrs_offset_sample);
    complex["t"] = py::make_tuple(ann.t_onset_sample, ann.t_peak_sample, ann.t_offset_sample);
    complex["p_full"] = py::make_tuple(
        ann.p1_onset_sample,
        ann.p1_peak_sample,
        ann.p1_offset_sample,
        ann.p2_onset_sample,
        ann.p2_peak_sample,
        ann.p2_offset_sample);
    complex["r_full"] = py::make_tuple(
        ann.qrs_onset_sample,
        ann.r_peak_sample,
        ann.qrs_offset_sample,
        ann.q_peak_sample,
        ann.s_peak_sample);
    complex["t_full"] = py::make_tuple(
        ann.t_onset_sample,
        ann.t_peak_sample,
        ann.t_offset_sample,
        ann.j_point_sample);
    return complex;
}

static bool has_any_pqrst_annotation(const rspt_pqrst_annotation& ann)
{
    return ann.p1_onset_sample >= 0 ||
           ann.p1_peak_sample >= 0 ||
           ann.p1_offset_sample >= 0 ||
           ann.p2_onset_sample >= 0 ||
           ann.p2_peak_sample >= 0 ||
           ann.p2_offset_sample >= 0 ||
           ann.qrs_onset_sample >= 0 ||
           ann.r_peak_sample >= 0 ||
           ann.qrs_offset_sample >= 0 ||
           ann.q_peak_sample >= 0 ||
           ann.s_peak_sample >= 0 ||
           ann.t_onset_sample >= 0 ||
           ann.t_peak_sample >= 0 ||
           ann.t_offset_sample >= 0 ||
           ann.j_point_sample >= 0;
}

static void add_standard_metric_fields(py::dict& output, const rspt_ecg_beat_result& result, const rspt_ecg_summary_result* summary)
{
    double qtc_bazett_ms = result.qtc_bazett_ms;
    if (summary && std::isfinite(summary->rr_interval_ms.mean) && summary->rr_interval_ms.mean > 0.0 && result.qt_interval_ms > 0.0)
        qtc_bazett_ms = result.qt_interval_ms / std::sqrt(summary->rr_interval_ms.mean / 1000.0);

    output["rr_interval_ms"]       = summary ? numeric_or_zero(summary->rr_interval_ms.mean) : numeric_or_zero(result.rr_interval_ms);
    output["rr_variation_ms"]      = summary ? numeric_or_zero(summary->rr_variation_ms) : 0.0;
    output["heart_rate_bpm"]       = summary ? numeric_or_zero(summary->heart_rate_bpm.mean) : numeric_or_zero(result.heart_rate_bpm);
    output["pr_interval_ms"]       = numeric_or_zero(result.pr_interval_ms);
    output["pr_segment_ms"]        = numeric_or_zero(result.pr_segment_ms);
    output["pp_interval_ms"]       = summary ? numeric_or_zero(summary->pp_interval_ms.mean) : numeric_or_zero(result.pp_interval_ms);
    output["p_peak_to_end_interval_ms"] = summary ? numeric_or_zero(summary->p_peak_to_end_interval_ms.mean) : numeric_or_zero(result.p_peak_to_end_interval_ms);
    output["qrs_duration_ms"]      = numeric_or_zero(result.qrs_duration_ms);
    output["qt_interval_ms"]       = numeric_or_zero(result.qt_interval_ms);
    output["qtc_bazett_ms"]        = numeric_or_zero(qtc_bazett_ms);
    output["qtc_interval_ms"]      = numeric_or_zero(qtc_bazett_ms);
    output["qt_dispersion_ms"]     = summary ? numeric_or_zero(summary->qt_dispersion_ms.mean) : numeric_or_zero(result.qt_dispersion_ms);
    output["st_segment_ms"]        = numeric_or_zero(result.st_segment_ms);
    output["p_wave_duration_ms"]   = numeric_or_zero(result.p_wave_duration_ms);
    output["t_wave_duration_ms"]   = numeric_or_zero(result.t_wave_duration_ms);
    output["t_peak_to_end_interval_ms"] = summary ? numeric_or_zero(summary->t_peak_to_end_interval_ms.mean) : numeric_or_zero(result.t_peak_to_end_interval_ms);
}

static void add_peak_interval_fields(py::dict& output, const std::vector<uint32_t>& peak_indexes, double sampling_rate, const rspt_ecg_beat_result& result)
{
    output["is_sinus_rhythm"] = 0;
    output["premature_beat_count"] = 0;

    if (peak_indexes.size() < 2 || !std::isfinite(sampling_rate) || sampling_rate <= 0.0)
        return;

    double total_ms = 0.0;
    double min_rr = std::numeric_limits<double>::max();
    double max_rr = 0.0;
    std::vector<double> rr_intervals;
    rr_intervals.reserve(peak_indexes.size() - 1);

    for (size_t i = 1; i < peak_indexes.size(); ++i)
    {
        if (peak_indexes[i] <= peak_indexes[i - 1])
            continue;

        double rr_ms = (peak_indexes[i] - peak_indexes[i - 1]) / sampling_rate * 1000.0;
        if (!std::isfinite(rr_ms) || rr_ms <= 0.0)
            continue;

        rr_intervals.push_back(rr_ms);
        total_ms += rr_ms;
        min_rr = std::min(min_rr, rr_ms);
        max_rr = std::max(max_rr, rr_ms);
    }

    if (rr_intervals.empty())
        return;

    double mean_rr = total_ms / static_cast<double>(rr_intervals.size());
    double rr_variation = max_rr - min_rr;
    int premature_beat_count = 0;
    for (double rr_ms : rr_intervals)
    {
        if (mean_rr > 0.0 && rr_ms < 0.80 * mean_rr)
            ++premature_beat_count;
    }

    output["rr_interval_ms"] = mean_rr;
    output["rr_variation_ms"] = rr_variation;
    output["heart_rate_bpm"] = mean_rr > 0.0 ? 60000.0 / mean_rr : 0.0;
    output["premature_beat_count"] = premature_beat_count;
    output["is_sinus_rhythm"] = (premature_beat_count || (mean_rr > 0.0 && rr_variation / mean_rr > 0.10)) ? 0 : 1;

    if (std::isfinite(result.qt_interval_ms) && result.qt_interval_ms > 0.0 && mean_rr > 0.0)
    {
        double qtc_bazett_ms = result.qt_interval_ms / std::sqrt(mean_rr / 1000.0);
        output["qtc_bazett_ms"] = qtc_bazett_ms;
        output["qtc_interval_ms"] = qtc_bazett_ms;
    }
}

static void add_validation_metric_fields(py::dict& output, const rspt_ecg_beat_result& result)
{
    output["p1_wave_duration_ms"] = numeric_or_zero(result.p1_wave_duration_ms);
    output["p1_amplitude_input_units"] = numeric_or_zero(result.p1_amplitude_input_units);
    output["p2_wave_duration_ms"] = numeric_or_zero(result.p2_wave_duration_ms);
    output["p2_amplitude_input_units"] = numeric_or_zero(result.p2_amplitude_input_units);
    output["q_duration_ms"] = numeric_or_zero(result.q_duration_ms);
    output["q_amplitude_input_units"] = numeric_or_zero(result.q_amplitude_input_units);
    output["r_duration_ms"] = numeric_or_zero(result.r_duration_ms);
    output["r_amplitude_input_units"] = numeric_or_zero(result.r_amplitude_input_units);
    output["s_duration_ms"] = numeric_or_zero(result.s_duration_ms);
    output["s_amplitude_input_units"] = numeric_or_zero(result.s_amplitude_input_units);
    output["j_point_amplitude_input_units"] = numeric_or_zero(result.j_point_amplitude_input_units);
    output["st20_amplitude_input_units"] = numeric_or_zero(result.st20_amplitude_input_units);
    output["st40_amplitude_input_units"] = numeric_or_zero(result.st40_amplitude_input_units);
    output["st60_amplitude_input_units"] = numeric_or_zero(result.st60_amplitude_input_units);
    output["st80_amplitude_input_units"] = numeric_or_zero(result.st80_amplitude_input_units);
    output["t_amplitude_input_units"] = numeric_or_zero(result.t_amplitude_input_units);

    output["P1_DURATION"] = numeric_or_zero(result.p1_wave_duration_ms);
    output["P1_AMPLITUDE"] = numeric_or_zero(result.p1_amplitude_input_units);
    output["P2_DURATION"] = numeric_or_zero(result.p2_wave_duration_ms);
    output["P2_AMPLITUDE"] = numeric_or_zero(result.p2_amplitude_input_units);
    output["PQ_INTERVAL"] = numeric_or_zero(result.pr_interval_ms);
    output["Q_DURATION"] = numeric_or_zero(result.q_duration_ms);
    output["Q_AMPLITUDE"] = numeric_or_zero(result.q_amplitude_input_units);
    output["R_DURATION"] = numeric_or_zero(result.r_duration_ms);
    output["R_AMPLITUDE"] = numeric_or_zero(result.r_amplitude_input_units);
    output["S_DURATION"] = numeric_or_zero(result.s_duration_ms);
    output["S_AMPLITUDE"] = numeric_or_zero(result.s_amplitude_input_units);
    output["QRS_DURATION"] = numeric_or_zero(result.qrs_duration_ms);
    output["QT_INTERVAL"] = numeric_or_zero(result.qt_interval_ms);
    output["J_AMPLITUDE"] = numeric_or_zero(result.j_point_amplitude_input_units);
    output["ST_20_AMPLITUDE"] = numeric_or_zero(result.st20_amplitude_input_units);
    output["ST_40_AMPLITUDE"] = numeric_or_zero(result.st40_amplitude_input_units);
    output["ST_60_AMPLITUDE"] = numeric_or_zero(result.st60_amplitude_input_units);
    output["ST_80_AMPLITUDE"] = numeric_or_zero(result.st80_amplitude_input_units);
    output["T_AMPLITUDE"] = numeric_or_zero(result.t_amplitude_input_units);
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

    ///cout << "Input ECG array must be 1 or 2-dimensional (samples x channels)" << endl;
    std::vector<uint32_t> peak_indexes;
    rspt::detect_peaks(ecg_signal.data(), len, sampling_rate, static_cast<rspt_detection_mode>(mode_to_c_api_value(mode)), peak_indexes);

    return std::vector<unsigned int>(peak_indexes.begin(), peak_indexes.end());
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

    std::vector<rspt_ecg_beat_result> results;
    rspt_ecg_summary_result summary;
    std::vector<uint32_t> peak_indexes;
    std::vector<uint32_t> beat_indexes_to_analyze;
    beat_indexes_to_analyze.push_back(static_cast<uint32_t>(std::max(0, analysis_peak_indx)));

    int32_t status = rspt::analyze_ecg(
        data_ptrs.data(),
        nr_channels,
        len,
        sampling_rate,
        results,
        summary,
        analysis_ch_indx,
        nullptr,
        &peak_indexes,
        static_cast<rspt_detection_mode>(mode_to_c_api_value(mode)),
        &beat_indexes_to_analyze);

    const rspt_ecg_beat_result* result = results.empty() ? nullptr : &results[0];

    py::list py_annotations;
    if (result && has_any_pqrst_annotation(result->annotation))
        py_annotations.append(annotation_to_dict(result->annotation));

    py::dict output;
    output["annotations"] = py_annotations;
    if (status == RSPT_STATUS_OK && result)
    {
        add_standard_metric_fields(output, *result, nullptr);
        add_peak_interval_fields(output, peak_indexes, sampling_rate, *result);
        output["analysis_status"] = result->status;
        output["status_message"] = std::string(rspt::status_message(result->status));
        add_validation_metric_fields(output, *result);
    }
    else
    {
        output["analysis_status"] = status;
        output["status_message"] = std::string(rspt::status_message(status));
    }

    return output;
}

py::dict analyse_ecg_beats(py::array_t<double, py::array::c_style | py::array::forcecast> ecg_signal_np, double sampling_rate, std::string mode = "default", int analysis_ch_indx = -1)
{
    size_t len = 0;
    size_t nr_channels = 0;
    std::vector<std::vector<double>> data;
    std::vector<const double*> data_ptrs;
    copy_ecg_signal_to_channels(ecg_signal_np, len, nr_channels, data, data_ptrs);

    std::vector<rspt_ecg_beat_result> results;
    rspt_ecg_summary_result summary;
    std::vector<uint32_t> peak_indexes;
    int32_t status = rspt::analyze_ecg(
        data_ptrs.data(),
        nr_channels,
        len,
        sampling_rate,
        results,
        summary,
        analysis_ch_indx,
        nullptr,
        &peak_indexes,
        static_cast<rspt_detection_mode>(mode_to_c_api_value(mode)),
        nullptr);

    py::list py_peak_indexes;
    for (const auto peak_index : peak_indexes)
        py_peak_indexes.append(peak_index);

    py::list py_beats;
    for (size_t i = 0; i < results.size(); ++i)
    {
        const auto& result = results[i];
        py::dict beat;
        beat["beat_index"] = i;
        if (result.r_peak_sample >= 0)
            beat["r_peak_sample"] = result.r_peak_sample;
        else
            beat["r_peak_sample"] = py::none();

        beat["analysis_channel_index"] = result.analysis_channel_index;
        beat["analysis_status"] = result.status;
        beat["status_message"] = std::string(rspt::status_message(result.status));
        beat["annotation"] = annotation_to_dict(result.annotation);

        add_standard_metric_fields(beat, result, nullptr);
        add_validation_metric_fields(beat, result);
        beat["rr_interval_ms"] = numeric_or_none(result.rr_interval_ms);
        beat["heart_rate_bpm"] = numeric_or_none(result.heart_rate_bpm);

        py_beats.append(beat);
    }

    py::dict output;
    output["analysis_channel_index"] = results.empty() ? analysis_ch_indx : static_cast<int>(results[0].analysis_channel_index);
    output["analysis_status"] = status;
    output["status_message"] = std::string(rspt::status_message(status));
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
