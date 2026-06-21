#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "rspt_types.h"

namespace rspt
{

const char* status_message(int32_t status);
uint32_t api_version();

int32_t detect_peaks_double(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    std::vector<uint32_t>& out_r_peak_indexes);

int32_t analyze_ecg_beats_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    std::vector<rspt_ecg_beat_result>& out_beats,
    std::vector<uint32_t>* out_detected_r_peak_indexes = nullptr);

int32_t analyze_ecg_beat_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    int32_t analysis_peak_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_beat_result& out_beat,
    std::vector<uint32_t>* out_detected_r_peak_indexes = nullptr);

int32_t analyze_ecg_summary_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_summary_result& out_summary);

}
