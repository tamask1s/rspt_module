#pragma once

#include "rspt_types.h"

#if defined(_WIN32)
#  if defined(RSPT_BUILD_SHARED)
#    define RSPT_API __declspec(dllexport)
#  else
#    define RSPT_API
#  endif
#else
#  define RSPT_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** Return a stable text description for an rspt_status value. */
RSPT_API const char* rspt_status_message(int32_t status);

/** Return the public API version of this library. */
RSPT_API uint32_t rspt_api_version(void);

/**
 * Detect R peaks in one ECG channel.
 *
 * The caller owns all buffers. To query the required output size, pass
 * out_r_peak_indexes as NULL with r_peak_capacity set to 0; out_r_peak_count
 * is still written. If out_r_peak_indexes is non-NULL and the capacity is too
 * small, RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL is returned after writing the
 * required count.
 *
 * @param signal ECG samples for one channel, stored as double precision values.
 * @param sample_count Number of samples in signal.
 * @param sampling_rate Sampling rate in Hz.
 * @param peak_detection_mode RSPT_MODE_DEFAULT is used for CTS/CSE validation; RSPT_MODE_HIGH_SENSITIVITY
 *        and RSPT_MODE_HIGH_PPV are optional detector modes.
 * @param out_r_peak_indexes Optional output buffer for detected R-peak sample indexes.
 * @param r_peak_capacity Number of uint32_t entries available in out_r_peak_indexes.
 * @param out_r_peak_count Required output. Receives the number of detected R peaks.
 * @return rspt_status value.
 */
RSPT_API int32_t rspt_detect_peaks(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    rspt_detection_mode peak_detection_mode,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count);

/**
 * Analyze ECG beat morphology, interval metrics, and summary statistics.
 *
 * This is the C wrapper for rspt::analyze_ecg. It is intended to cover PDF,
 * multi-beat Flutter analysis, and CTS/CSE validation. If input_r_peak_indexes
 * is NULL or input_r_peak_count is 0, R peaks are detected first. If
 * beat_indexes_to_analyze is NULL or beat_index_count is 0, every available
 * beat is analyzed; otherwise only the selected beat indexes are analyzed.
 *
 * Beat morphology is measured on analysis_channel_index. QT dispersion is the
 * per-beat multichannel QT spread max(QT) - min(QT) across available channels,
 * so it requires at least two analyzable channels. The current implementation
 * computes QT dispersion in all-beat mode and for explicit selected-beat
 * requests with at least two analyzed beats; otherwise per-beat QT dispersion
 * is left invalid.
 * Result structs use NaN for unavailable double metrics and -1 for unavailable
 * sample indexes. Summary metric statistics use valid_count to report how many
 * valid beat values contributed to mean/std.
 *
 * The caller owns all buffers. To query output sizes, pass out_beat_results
 * and/or out_r_peak_indexes as NULL with the corresponding capacity set to 0.
 * out_beat_result_count is required and receives the number of analyzed beats.
 * out_r_peak_count is optional and receives the number of effective R peaks.
 * If a non-NULL output buffer is too small, RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL
 * is returned after writing the required count fields.
 *
 * @param channels ECG channel array. channels[ch][sample] must be valid for every channel and sample.
 * @param channel_count Number of ECG channels in channels.
 * @param samples_per_channel Number of samples in each channel.
 * @param sampling_rate Sampling rate in Hz.
 * @param out_beat_results Optional output buffer, one rspt_ecg_beat_result for each analyzed beat.
 * @param beat_result_capacity Number of entries available in out_beat_results.
 * @param out_beat_result_count Required output. Receives the number of analyzed beat results.
 * @param out_summary Required output aggregate metrics. RR/heart-rate fields are derived from the effective
 *        R-peak list; morphology and per-beat QT dispersion fields are aggregated over out_beat_results.
 * @param analysis_channel_index -1 lets the library choose the analysis channel; otherwise this is a zero-based channel index.
 * @param input_r_peak_indexes Optional input R-peak sample indexes.
 * @param input_r_peak_count Number of entries in input_r_peak_indexes.
 * @param out_r_peak_indexes Optional output buffer for the effective R-peak sample indexes.
 * @param r_peak_capacity Number of uint32_t entries available in out_r_peak_indexes.
 * @param out_r_peak_count Optional output. Receives the number of effective R peaks.
 * @param peak_detection_mode RSPT_MODE_DEFAULT is used for CTS/CSE validation; used only when no input R peaks are provided.
 * @param beat_indexes_to_analyze Optional indexes into the effective R-peak list.
 * @param beat_index_count Number of entries in beat_indexes_to_analyze.
 * @return rspt_status value.
 */
RSPT_API int32_t rspt_analyze_ecg(
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
    size_t beat_index_count);

#ifdef __cplusplus
}
#endif
