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

RSPT_API int32_t rspt_detect_peaks(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    rspt_detection_mode peak_detection_mode,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count);

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
