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

    std::vector<uint32_t> peak_indexes;
    int32_t status = rspt::detect_peaks_double(signal, sample_count, sampling_rate, mode, peak_indexes);
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

int32_t rspt_analyze_ecg_beats_double(
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

    std::vector<rspt_ecg_beat_result> beats;
    std::vector<uint32_t> detected_r_peak_indexes;
    int32_t status = rspt::analyze_ecg_beats_double(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis_channel_index,
        r_peak_indexes,
        r_peak_count,
        mode,
        beats,
        &detected_r_peak_indexes);

    *out_beat_count = beats.size();
    if (out_detected_r_peak_count)
        *out_detected_r_peak_count = detected_r_peak_indexes.size();

    if (status != RSPT_STATUS_OK)
        return status;
    if (out_beats && beat_capacity < beats.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;
    if (out_detected_r_peak_indexes && detected_r_peak_capacity < detected_r_peak_indexes.size())
        return RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL;

    if (out_beats)
    {
        for (size_t i = 0; i < beats.size(); ++i)
            out_beats[i] = beats[i];
    }

    if (out_detected_r_peak_indexes)
    {
        for (size_t i = 0; i < detected_r_peak_indexes.size(); ++i)
            out_detected_r_peak_indexes[i] = detected_r_peak_indexes[i];
    }

    return RSPT_STATUS_OK;
}

int32_t rspt_analyze_ecg_summary_double(
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

    return rspt::analyze_ecg_summary_double(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis_channel_index,
        r_peak_indexes,
        r_peak_count,
        mode,
        *out_summary);
}

}
