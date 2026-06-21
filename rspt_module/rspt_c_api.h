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

RSPT_API const char* rspt_status_message(int32_t status);
RSPT_API uint32_t rspt_api_version(void);

RSPT_API int32_t rspt_detect_peaks_double(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count);

RSPT_API int32_t rspt_analyze_ecg_beats_double(
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
    size_t* out_detected_r_peak_count);

RSPT_API int32_t rspt_analyze_ecg_beat_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    int32_t analysis_peak_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_beat_result* out_beat,
    uint32_t* out_detected_r_peak_indexes,
    size_t detected_r_peak_capacity,
    size_t* out_detected_r_peak_count);

RSPT_API int32_t rspt_analyze_ecg_summary_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_summary_result* out_summary);

#ifdef __cplusplus
}
#endif
