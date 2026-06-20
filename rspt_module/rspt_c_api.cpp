#include "rspt_c_api.h"

#include "ecg_analysis.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

using namespace std;

#include "filter.h"
#include "iir_filter_opt.h"
#include "peak_detector.h"

namespace
{
constexpr size_t kStatusMessageSize = 64;

double missing_value()
{
    return std::numeric_limits<double>::quiet_NaN();
}

bool is_valid_mode(int32_t mode)
{
    return mode == RSPT_MODE_DEFAULT ||
           mode == RSPT_MODE_HIGH_SENSITIVITY ||
           mode == RSPT_MODE_HIGH_PPV;
}

void configure_detector(peak_detector_offline& detector, int32_t mode)
{
    if (mode == RSPT_MODE_HIGH_SENSITIVITY)
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == RSPT_MODE_HIGH_PPV)
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);
}

void copy_message(char* dest, const char* message)
{
    if (!dest)
        return;
    std::snprintf(dest, kStatusMessageSize, "%s", message ? message : "");
}

bool positive_finite(double value)
{
    return std::isfinite(value) && value > 0.0;
}

bool non_negative_finite(double value)
{
    return std::isfinite(value) && value >= 0.0;
}

void set_if_positive(double value, double& field, uint64_t& valid_fields, uint64_t bit)
{
    if (positive_finite(value))
    {
        field = value;
        valid_fields |= bit;
    }
}

void set_if_non_negative(double value, double& field, uint64_t& valid_fields, uint64_t bit)
{
    if (non_negative_finite(value))
    {
        field = value;
        valid_fields |= bit;
    }
}

void initialise_annotation(rspt_pqrst_annotation& annotation)
{
    annotation.p1_onset_sample = -1;
    annotation.p1_peak_sample = -1;
    annotation.p1_offset_sample = -1;
    annotation.p2_onset_sample = -1;
    annotation.p2_peak_sample = -1;
    annotation.p2_offset_sample = -1;
    annotation.qrs_onset_sample = -1;
    annotation.r_peak_sample = -1;
    annotation.qrs_offset_sample = -1;
    annotation.q_peak_sample = -1;
    annotation.s_peak_sample = -1;
    annotation.t_onset_sample = -1;
    annotation.t_peak_sample = -1;
    annotation.t_offset_sample = -1;
    annotation.j_point_sample = -1;
}

void initialise_beat_result(rspt_ecg_beat_result& result)
{
    std::memset(&result, 0, sizeof(result));
    result.status = RSPT_STATUS_OK;
    copy_message(result.status_message, "OK");
    initialise_annotation(result.annotation);

    result.rr_interval_ms = missing_value();
    result.heart_rate_bpm = missing_value();
    result.p_wave_duration_ms = missing_value();
    result.p1_wave_duration_ms = missing_value();
    result.p1_amplitude_input_units = missing_value();
    result.p2_wave_duration_ms = missing_value();
    result.p2_amplitude_input_units = missing_value();
    result.pr_interval_ms = missing_value();
    result.pr_segment_ms = missing_value();
    result.q_duration_ms = missing_value();
    result.q_amplitude_input_units = missing_value();
    result.r_duration_ms = missing_value();
    result.r_amplitude_input_units = missing_value();
    result.s_duration_ms = missing_value();
    result.s_amplitude_input_units = missing_value();
    result.qrs_duration_ms = missing_value();
    result.qt_interval_ms = missing_value();
    result.qtc_bazett_ms = missing_value();
    result.st_segment_ms = missing_value();
    result.t_wave_duration_ms = missing_value();
    result.j_point_amplitude_input_units = missing_value();
    result.st20_amplitude_input_units = missing_value();
    result.st40_amplitude_input_units = missing_value();
    result.st60_amplitude_input_units = missing_value();
    result.st80_amplitude_input_units = missing_value();
    result.t_amplitude_input_units = missing_value();
}

void initialise_summary_result(rspt_ecg_summary_result& result)
{
    std::memset(&result, 0, sizeof(result));
    result.status = RSPT_STATUS_OK;
    copy_message(result.status_message, "OK");

    result.mean_rr_interval_ms = missing_value();
    result.rr_variation_ms = missing_value();
    result.heart_rate_bpm = missing_value();
    result.mean_p_wave_duration_ms = missing_value();
    result.mean_pr_interval_ms = missing_value();
    result.mean_qrs_duration_ms = missing_value();
    result.mean_qt_interval_ms = missing_value();
    result.mean_qtc_bazett_ms = missing_value();
    result.mean_st_segment_ms = missing_value();
    result.mean_t_wave_duration_ms = missing_value();
    result.is_sinus_rhythm = -1;
    result.premature_beat_count = -1;
}

int32_t detect_peak_indexes(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    std::vector<unsigned int>& peak_indexes)
{
    peak_indexes.clear();

    if (!signal || sample_count == 0 || !std::isfinite(sampling_rate) || sampling_rate <= 0.0)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (sample_count > static_cast<size_t>(std::numeric_limits<unsigned int>::max()))
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!is_valid_mode(mode))
        return RSPT_STATUS_INVALID_ARGUMENT;

    std::vector<double> peak_signal(sample_count);
    std::vector<double> filtered_signal(sample_count);
    std::vector<double> threshold_signal(sample_count);

    peak_detector_offline detector(sampling_rate);
    configure_detector(detector, mode);
    detector.detect(
        signal,
        static_cast<unsigned int>(sample_count),
        peak_signal.data(),
        filtered_signal.data(),
        threshold_signal.data(),
        &peak_indexes);

    std::sort(peak_indexes.begin(), peak_indexes.end());
    peak_indexes.erase(std::unique(peak_indexes.begin(), peak_indexes.end()), peak_indexes.end());

    return peak_indexes.empty() ? RSPT_STATUS_NO_R_PEAKS : RSPT_STATUS_OK;
}

struct AnalysisData
{
    int32_t analysis_channel_index = 0;
    std::vector<unsigned int> peak_indexes;
    std::vector<pqrst_indxes> annotations;
    std::vector<ecg_analysis_result> results;
};

int32_t run_core_analysis(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    AnalysisData& analysis)
{
    analysis = AnalysisData{};

    if (!channels || channel_count == 0 || samples_per_channel == 0 || !std::isfinite(sampling_rate) || sampling_rate <= 0.0)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (channel_count > static_cast<size_t>(std::numeric_limits<unsigned int>::max()) ||
        samples_per_channel > static_cast<size_t>(std::numeric_limits<unsigned int>::max()))
        return RSPT_STATUS_INVALID_ARGUMENT;

    for (size_t channel = 0; channel < channel_count; ++channel)
    {
        if (!channels[channel])
            return RSPT_STATUS_NO_CHANNEL_DATA;
    }

    if (analysis_channel_index < 0)
        analysis_channel_index = 0;
    if (static_cast<size_t>(analysis_channel_index) >= channel_count)
        return RSPT_STATUS_INVALID_CHANNEL_INDEX;
    analysis.analysis_channel_index = analysis_channel_index;

    if (r_peak_indexes && r_peak_count > 0)
    {
        analysis.peak_indexes.reserve(r_peak_count);
        for (size_t i = 0; i < r_peak_count; ++i)
        {
            if (r_peak_indexes[i] >= samples_per_channel)
                return RSPT_STATUS_INVALID_PEAK_INDEX;
            analysis.peak_indexes.push_back(static_cast<unsigned int>(r_peak_indexes[i]));
        }
        std::sort(analysis.peak_indexes.begin(), analysis.peak_indexes.end());
        analysis.peak_indexes.erase(std::unique(analysis.peak_indexes.begin(), analysis.peak_indexes.end()), analysis.peak_indexes.end());
    }
    else
    {
        int32_t status = detect_peak_indexes(
            channels[analysis.analysis_channel_index],
            samples_per_channel,
            sampling_rate,
            mode,
            analysis.peak_indexes);
        if (status != RSPT_STATUS_OK)
            return status;
    }

    if (analysis.peak_indexes.empty())
        return RSPT_STATUS_NO_R_PEAKS;

    std::vector<const double*> channel_ptrs(channel_count);
    for (size_t channel = 0; channel < channel_count; ++channel)
        channel_ptrs[channel] = channels[channel];

    analyse_ecg_all_beats(
        channel_ptrs.data(),
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis.peak_indexes,
        analysis.annotations,
        analysis.results,
        analysis.analysis_channel_index);

    return RSPT_STATUS_OK;
}

void copy_annotation(const pqrst_indxes& source, rspt_pqrst_annotation& destination, uint64_t& valid_fields)
{
    initialise_annotation(destination);

    destination.p1_onset_sample = source.p[0];
    destination.p1_peak_sample = source.p[1];
    destination.p1_offset_sample = source.p[2];

    if (source.p[3] >= 0 && source.p[4] >= 0 && source.p[5] >= 0)
    {
        destination.p2_onset_sample = source.p[3];
        destination.p2_peak_sample = source.p[4];
        destination.p2_offset_sample = source.p[5];
    }

    destination.qrs_onset_sample = source.r[0];
    destination.r_peak_sample = source.r[1];
    destination.qrs_offset_sample = source.r[2];
    destination.q_peak_sample = source.r[3];
    destination.s_peak_sample = source.r[4];

    destination.t_onset_sample = source.t[0];
    destination.t_peak_sample = source.t[1];
    destination.t_offset_sample = source.t[2];
    destination.j_point_sample = source.t[3];

    if (source.p[0] >= 0 && source.p[1] >= 0 && source.p[2] >= 0 &&
        source.r[0] >= 0 && source.r[1] >= 0 && source.r[2] >= 0 &&
        source.t[0] >= 0 && source.t[1] >= 0 && source.t[2] >= 0)
    {
        valid_fields |= RSPT_VALID_PQRS_T_ANNOTATION;
    }
}

bool rr_interval_ms_for_index(
    const std::vector<unsigned int>& peak_indexes,
    size_t right_peak_index,
    double sampling_rate,
    double& rr_ms)
{
    if (right_peak_index == 0 || right_peak_index >= peak_indexes.size())
        return false;
    if (peak_indexes[right_peak_index] <= peak_indexes[right_peak_index - 1])
        return false;

    rr_ms = (peak_indexes[right_peak_index] - peak_indexes[right_peak_index - 1]) / sampling_rate * 1000.0;
    return positive_finite(rr_ms);
}

void fill_beat_result(
    const ecg_analysis_result& core_result,
    const pqrst_indxes* annotation,
    const std::vector<unsigned int>& peak_indexes,
    size_t beat_index,
    double sampling_rate,
    uint32_t analysis_channel_index,
    rspt_ecg_beat_result& output)
{
    initialise_beat_result(output);
    output.status = core_result.analysis_status;
    output.analysis_channel_index = analysis_channel_index;
    output.beat_index = static_cast<uint32_t>(beat_index);

    if (beat_index < peak_indexes.size())
    {
        output.r_peak_sample = peak_indexes[beat_index];
        output.valid_fields |= RSPT_VALID_R_PEAK_SAMPLE;
    }

    copy_message(output.status_message, core_result.status_message[0] ? core_result.status_message : rspt_status_message(output.status));

    if (core_result.analysis_status != RSPT_STATUS_OK)
        return;

    if (annotation)
        copy_annotation(*annotation, output.annotation, output.valid_fields);

    set_if_positive(core_result.p_wave_duration_ms, output.p_wave_duration_ms, output.valid_fields, RSPT_VALID_P_WAVE_DURATION_MS);
    output.p1_wave_duration_ms = core_result.p1_wave_duration_ms;
    output.p1_amplitude_input_units = core_result.p1_amplitude_input_units;
    output.p2_wave_duration_ms = core_result.p2_wave_duration_ms;
    output.p2_amplitude_input_units = core_result.p2_amplitude_input_units;

    set_if_positive(core_result.pr_interval_ms, output.pr_interval_ms, output.valid_fields, RSPT_VALID_PR_INTERVAL_MS);
    set_if_non_negative(core_result.pr_segment_ms, output.pr_segment_ms, output.valid_fields, RSPT_VALID_PR_SEGMENT_MS);

    output.q_duration_ms = core_result.q_duration_ms;
    output.q_amplitude_input_units = core_result.q_amplitude_input_units;
    output.r_duration_ms = core_result.r_duration_ms;
    output.r_amplitude_input_units = core_result.r_amplitude_input_units;
    output.s_duration_ms = core_result.s_duration_ms;
    output.s_amplitude_input_units = core_result.s_amplitude_input_units;

    set_if_positive(core_result.qrs_duration_ms, output.qrs_duration_ms, output.valid_fields, RSPT_VALID_QRS_DURATION_MS);
    set_if_positive(core_result.qt_interval_ms, output.qt_interval_ms, output.valid_fields, RSPT_VALID_QT_INTERVAL_MS);
    set_if_non_negative(core_result.st_segment_ms, output.st_segment_ms, output.valid_fields, RSPT_VALID_ST_SEGMENT_MS);
    set_if_positive(core_result.t_wave_duration_ms, output.t_wave_duration_ms, output.valid_fields, RSPT_VALID_T_WAVE_DURATION_MS);

    output.j_point_amplitude_input_units = core_result.j_point_amplitude_input_units;
    output.st20_amplitude_input_units = core_result.st20_amplitude_input_units;
    output.st40_amplitude_input_units = core_result.st40_amplitude_input_units;
    output.st60_amplitude_input_units = core_result.st60_amplitude_input_units;
    output.st80_amplitude_input_units = core_result.st80_amplitude_input_units;
    output.t_amplitude_input_units = core_result.t_amplitude_input_units;

    double previous_rr_ms = 0.0;
    if (rr_interval_ms_for_index(peak_indexes, beat_index, sampling_rate, previous_rr_ms))
    {
        output.rr_interval_ms = previous_rr_ms;
        output.heart_rate_bpm = 60000.0 / previous_rr_ms;
        output.valid_fields |= RSPT_VALID_RR_INTERVAL_MS | RSPT_VALID_HEART_RATE_BPM;
        if (positive_finite(output.qt_interval_ms))
        {
            output.qtc_bazett_ms = output.qt_interval_ms / std::sqrt(previous_rr_ms / 1000.0);
            output.valid_fields |= RSPT_VALID_QTC_BAZETT_MS;
        }
    }
}

void add_to_mean(double value, double& sum, size_t& count)
{
    if (positive_finite(value))
    {
        sum += value;
        ++count;
    }
}

void set_mean(double sum, size_t count, double& destination, uint64_t& valid_fields, uint64_t bit)
{
    if (count > 0)
    {
        destination = sum / static_cast<double>(count);
        valid_fields |= bit;
    }
}
}

extern "C"
{

const char* rspt_status_message(int32_t status)
{
    switch (status)
    {
        case RSPT_STATUS_OK:
            return "OK";
        case RSPT_STATUS_NO_R_PEAKS:
            return "No R peaks found";
        case RSPT_STATUS_NO_CHANNEL_DATA:
            return "No channel data given";
        case RSPT_STATUS_INVALID_ARGUMENT:
            return "Invalid argument";
        case RSPT_STATUS_INVALID_CHANNEL_INDEX:
            return "Invalid channel index";
        case RSPT_STATUS_INVALID_PEAK_INDEX:
            return "Invalid peak index";
        case RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL:
            return "Output buffer too small";
        default:
            return "Unknown status";
    }
}

uint32_t rspt_c_api_version(void)
{
    return 3;
}

int32_t rspt_detect_peaks_double(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count)
{
    if (!out_r_peak_count)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!out_r_peak_indexes && r_peak_capacity > 0)
        return RSPT_STATUS_INVALID_ARGUMENT;

    std::vector<unsigned int> peak_indexes;
    int32_t status = detect_peak_indexes(signal, sample_count, sampling_rate, mode, peak_indexes);
    *out_r_peak_count = peak_indexes.size();

    if (status != RSPT_STATUS_OK)
        return status;
    if (out_r_peak_indexes && r_peak_capacity < peak_indexes.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;

    if (out_r_peak_indexes)
    {
        for (size_t i = 0; i < peak_indexes.size(); ++i)
            out_r_peak_indexes[i] = peak_indexes[i];
    }

    return RSPT_STATUS_OK;
}

int32_t rspt_analyse_ecg_beats_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_beat_result* out_beats,
    size_t beat_capacity,
    size_t* out_beat_count,
    uint32_t* out_detected_r_peak_indexes,
    size_t detected_r_peak_capacity,
    size_t* out_detected_r_peak_count)
{
    if (!out_beat_count)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!out_beats && beat_capacity > 0)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!out_detected_r_peak_indexes && detected_r_peak_capacity > 0)
        return RSPT_STATUS_INVALID_ARGUMENT;

    AnalysisData analysis;
    int32_t status = run_core_analysis(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis_channel_index,
        r_peak_indexes,
        r_peak_count,
        mode,
        analysis);

    *out_beat_count = analysis.results.size();
    if (out_detected_r_peak_count)
        *out_detected_r_peak_count = analysis.peak_indexes.size();

    if (status != RSPT_STATUS_OK)
        return status;

    if (out_beats && beat_capacity < analysis.results.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;
    if (out_detected_r_peak_indexes && detected_r_peak_capacity < analysis.peak_indexes.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;

    if (out_detected_r_peak_indexes)
    {
        for (size_t i = 0; i < analysis.peak_indexes.size(); ++i)
            out_detected_r_peak_indexes[i] = analysis.peak_indexes[i];
    }

    if (out_beats)
    {
        for (size_t i = 0; i < analysis.results.size(); ++i)
        {
            const pqrst_indxes* annotation = (i < analysis.annotations.size()) ? &analysis.annotations[i] : nullptr;
            fill_beat_result(
                analysis.results[i],
                annotation,
                analysis.peak_indexes,
                i,
                sampling_rate,
                static_cast<uint32_t>(analysis.analysis_channel_index),
                out_beats[i]);
        }
    }

    return RSPT_STATUS_OK;
}

int32_t rspt_analyse_ecg_summary_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_summary_result* out_summary)
{
    if (!out_summary)
        return RSPT_STATUS_INVALID_ARGUMENT;

    initialise_summary_result(*out_summary);

    AnalysisData analysis;
    int32_t status = run_core_analysis(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis_channel_index,
        r_peak_indexes,
        r_peak_count,
        mode,
        analysis);

    out_summary->status = status;
    out_summary->analysis_channel_index = static_cast<uint32_t>(analysis.analysis_channel_index);
    out_summary->r_peak_count = analysis.peak_indexes.size();
    copy_message(out_summary->status_message, rspt_status_message(status));

    if (status != RSPT_STATUS_OK)
        return status;

    if (analysis.peak_indexes.size() >= 2)
    {
        double rr_sum = 0.0;
        double rr_min = std::numeric_limits<double>::max();
        double rr_max = 0.0;
        std::vector<double> rr_intervals;
        rr_intervals.reserve(analysis.peak_indexes.size() - 1);

        for (size_t i = 1; i < analysis.peak_indexes.size(); ++i)
        {
            double rr_ms = 0.0;
            if (!rr_interval_ms_for_index(analysis.peak_indexes, i, sampling_rate, rr_ms))
                continue;
            rr_intervals.push_back(rr_ms);
            rr_sum += rr_ms;
            rr_min = std::min(rr_min, rr_ms);
            rr_max = std::max(rr_max, rr_ms);
        }

        if (!rr_intervals.empty())
        {
            out_summary->mean_rr_interval_ms = rr_sum / static_cast<double>(rr_intervals.size());
            out_summary->rr_variation_ms = rr_max - rr_min;
            out_summary->heart_rate_bpm = 60000.0 / out_summary->mean_rr_interval_ms;
            out_summary->valid_fields |= RSPT_VALID_RR_INTERVAL_MS | RSPT_VALID_HEART_RATE_BPM;

            int32_t premature_count = 0;
            for (double rr_ms : rr_intervals)
            {
                if (rr_ms < 0.80 * out_summary->mean_rr_interval_ms)
                    ++premature_count;
            }
            out_summary->premature_beat_count = premature_count;
            out_summary->is_sinus_rhythm =
                (premature_count == 0 &&
                 out_summary->mean_rr_interval_ms > 0.0 &&
                 out_summary->rr_variation_ms / out_summary->mean_rr_interval_ms <= 0.10)
                ? 1
                : 0;
        }
    }

    double p_sum = 0.0, pr_sum = 0.0, qrs_sum = 0.0, qt_sum = 0.0, qtc_sum = 0.0, st_sum = 0.0, t_sum = 0.0;
    size_t p_count = 0, pr_count = 0, qrs_count = 0, qt_count = 0, qtc_count = 0, st_count = 0, t_count = 0;
    int32_t first_failure_status = RSPT_STATUS_OK;

    for (size_t i = 0; i < analysis.results.size(); ++i)
    {
        if (analysis.results[i].analysis_status != RSPT_STATUS_OK)
        {
            if (first_failure_status == RSPT_STATUS_OK)
                first_failure_status = analysis.results[i].analysis_status;
            continue;
        }

        ++out_summary->analysed_beat_count;

        rspt_ecg_beat_result beat;
        const pqrst_indxes* annotation = (i < analysis.annotations.size()) ? &analysis.annotations[i] : nullptr;
        fill_beat_result(
            analysis.results[i],
            annotation,
            analysis.peak_indexes,
            i,
            sampling_rate,
            static_cast<uint32_t>(analysis.analysis_channel_index),
            beat);

        add_to_mean(beat.p_wave_duration_ms, p_sum, p_count);
        add_to_mean(beat.pr_interval_ms, pr_sum, pr_count);
        add_to_mean(beat.qrs_duration_ms, qrs_sum, qrs_count);
        add_to_mean(beat.qt_interval_ms, qt_sum, qt_count);
        add_to_mean(beat.qtc_bazett_ms, qtc_sum, qtc_count);
        add_to_mean(beat.st_segment_ms, st_sum, st_count);
        add_to_mean(beat.t_wave_duration_ms, t_sum, t_count);
    }

    if (out_summary->analysed_beat_count == 0)
    {
        out_summary->status = first_failure_status == RSPT_STATUS_OK ? RSPT_STATUS_INVALID_PEAK_INDEX : first_failure_status;
        copy_message(out_summary->status_message, rspt_status_message(out_summary->status));
        return out_summary->status;
    }

    set_mean(p_sum, p_count, out_summary->mean_p_wave_duration_ms, out_summary->valid_fields, RSPT_VALID_P_WAVE_DURATION_MS);
    set_mean(pr_sum, pr_count, out_summary->mean_pr_interval_ms, out_summary->valid_fields, RSPT_VALID_PR_INTERVAL_MS);
    set_mean(qrs_sum, qrs_count, out_summary->mean_qrs_duration_ms, out_summary->valid_fields, RSPT_VALID_QRS_DURATION_MS);
    set_mean(qt_sum, qt_count, out_summary->mean_qt_interval_ms, out_summary->valid_fields, RSPT_VALID_QT_INTERVAL_MS);
    set_mean(qtc_sum, qtc_count, out_summary->mean_qtc_bazett_ms, out_summary->valid_fields, RSPT_VALID_QTC_BAZETT_MS);
    set_mean(st_sum, st_count, out_summary->mean_st_segment_ms, out_summary->valid_fields, RSPT_VALID_ST_SEGMENT_MS);
    set_mean(t_sum, t_count, out_summary->mean_t_wave_duration_ms, out_summary->valid_fields, RSPT_VALID_T_WAVE_DURATION_MS);

    out_summary->status = RSPT_STATUS_OK;
    copy_message(out_summary->status_message, "OK");
    return RSPT_STATUS_OK;
}

}
