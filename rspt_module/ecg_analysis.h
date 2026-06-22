#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "rspt_types.h"

namespace rspt
{

/** Return a stable text description for an rspt_status value. */
const char* status_message(int32_t status);

/** Return the public API version of this library. */
uint32_t api_version();

/**
 * Detect R peaks in one ECG channel.
 *
 * @param signal ECG samples for one channel, stored as double precision values.
 * @param sample_count Number of samples in signal.
 * @param sampling_rate Sampling rate in Hz.
 * @param peak_detection_mode RSPT_MODE_DEFAULT is used for CTS/CSE validation; RSPT_MODE_HIGH_SENSITIVITY
 *        and RSPT_MODE_HIGH_PPV are optional detector modes.
 * @param out_r_peak_indexes Output sample indexes of detected R peaks.
 * @return rspt_status value.
 */
int32_t detect_peaks(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    rspt_detection_mode peak_detection_mode,
    std::vector<uint32_t>& out_r_peak_indexes);

/**
 * Analyze ECG beat morphology, interval metrics, and summary statistics.
 *
 * This single analysis entry point is intended to cover PDF, multi-beat Flutter analysis, and CTS/CSE validation.
 * If input_r_peak_indexes is not provided, R peaks are detected first. If beat_indexes_to_analyze is not provided
 * or is empty, every available beat is analyzed; otherwise only the selected beat indexes are analyzed.
 * Beat morphology is measured on analysis_channel_index. QT dispersion is the per-beat multichannel QT spread
 * max(QT) - min(QT) across available channels, so it requires at least two analyzable channels. The current
 * implementation computes QT dispersion in all-beat mode and for explicit selected-beat requests with at least
 * two analyzed beats; otherwise per-beat QT dispersion is left invalid.
 * Result structs use NaN for unavailable double metrics and -1 for unavailable sample indexes.
 * Summary metric statistics use valid_count to report how many valid beat values contributed to mean/std.
 *
 * @param channels ECG channel array. channels[ch][sample] must be valid for every channel and sample.
 * @param channel_count Number of ECG channels in channels.
 * @param samples_per_channel Number of samples in each channel.
 * @param sampling_rate Sampling rate in Hz.
 * @param out_beat_results Output beat results, one entry for each analyzed beat.
 * @param out_summary Output aggregate metrics. RR/heart-rate fields are derived from the effective R-peak list;
 *        morphology and per-beat QT dispersion fields are aggregated over out_beat_results.
 * @param analysis_channel_index -1 lets the library choose the analysis channel; otherwise this is a zero-based channel index.
 * @param input_r_peak_indexes Optional input R-peak sample indexes. If null or empty, peaks are detected internally.
 * @param out_r_peak_indexes Optional output of the effective R-peak sample indexes, whether detected internally or supplied by input_r_peak_indexes.
 * @param peak_detection_mode RSPT_MODE_DEFAULT is used for CTS/CSE validation; used only when input_r_peak_indexes is null or empty.
 * @param beat_indexes_to_analyze Optional indexes into the effective R-peak list. If null or empty, all beats are analyzed.
 * @return rspt_status value.
 */
int32_t analyze_ecg(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    std::vector<rspt_ecg_beat_result>& out_beat_results,
    rspt_ecg_summary_result& out_summary,
    int32_t analysis_channel_index = -1,
    const std::vector<uint32_t>* input_r_peak_indexes = nullptr,
    std::vector<uint32_t>* out_r_peak_indexes = nullptr,
    rspt_detection_mode peak_detection_mode = RSPT_MODE_DEFAULT,
    const std::vector<uint32_t>* beat_indexes_to_analyze = nullptr);

}
