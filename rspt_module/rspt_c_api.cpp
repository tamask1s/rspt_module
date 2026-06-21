#include "rspt_c_api.h"

#include "ecg_analysis.h"

#include <vector>

extern "C"
{

const char* rspt_status_message(int32_t status)
{
    return rspt::status_message(status);
}

uint32_t rspt_api_version(void)
{
    return rspt::api_version();
}

int32_t rspt_detect_peaks(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    rspt_detection_mode peak_detection_mode,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count)
{
    if (!out_r_peak_count)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!out_r_peak_indexes && r_peak_capacity > 0)
        return RSPT_STATUS_INVALID_ARGUMENT;

    std::vector<uint32_t> peak_indexes;
    int32_t status = rspt::detect_peaks(signal, sample_count, sampling_rate, peak_detection_mode, peak_indexes);
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

int32_t rspt_analyze_ecg(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    rspt_ecg_beat_result* out_beat_results,
    size_t beat_result_capacity,
    size_t* out_beat_result_count,
    rspt_ecg_summary_result* out_summary,
    int32_t analysis_channel_index,
    const uint32_t* input_r_peak_indexes,
    size_t input_r_peak_count,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count,
    rspt_detection_mode peak_detection_mode,
    const uint32_t* beat_indexes_to_analyze,
    size_t beat_index_count)
{
    if (!out_beat_result_count || !out_summary)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!out_beat_results && beat_result_capacity > 0)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!out_r_peak_indexes && r_peak_capacity > 0)
        return RSPT_STATUS_INVALID_ARGUMENT;

    std::vector<uint32_t> input_peaks;
    if (input_r_peak_indexes && input_r_peak_count > 0)
        input_peaks.assign(input_r_peak_indexes, input_r_peak_indexes + input_r_peak_count);

    std::vector<uint32_t> selected_beats;
    if (beat_indexes_to_analyze && beat_index_count > 0)
        selected_beats.assign(beat_indexes_to_analyze, beat_indexes_to_analyze + beat_index_count);

    std::vector<rspt_ecg_beat_result> beats;
    std::vector<uint32_t> effective_r_peak_indexes;
    int32_t status = rspt::analyze_ecg(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        beats,
        *out_summary,
        analysis_channel_index,
        input_peaks.empty() ? nullptr : &input_peaks,
        &effective_r_peak_indexes,
        peak_detection_mode,
        selected_beats.empty() ? nullptr : &selected_beats);

    *out_beat_result_count = beats.size();
    if (out_r_peak_count)
        *out_r_peak_count = effective_r_peak_indexes.size();

    if (status != RSPT_STATUS_OK)
        return status;
    if (out_beat_results && beat_result_capacity < beats.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;
    if (out_r_peak_indexes && r_peak_capacity < effective_r_peak_indexes.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;

    if (out_beat_results)
    {
        for (size_t i = 0; i < beats.size(); ++i)
            out_beat_results[i] = beats[i];
    }

    if (out_r_peak_indexes)
    {
        for (size_t i = 0; i < effective_r_peak_indexes.size(); ++i)
            out_r_peak_indexes[i] = effective_r_peak_indexes[i];
    }

    return RSPT_STATUS_OK;
}

}
