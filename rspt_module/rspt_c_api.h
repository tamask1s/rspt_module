#pragma once

#include <stddef.h>
#include <stdint.h>

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

typedef enum rspt_status {
    RSPT_STATUS_OK = 0,
    RSPT_STATUS_NO_R_PEAKS = 1,
    RSPT_STATUS_NO_CHANNEL_DATA = 2,
    RSPT_STATUS_INVALID_ARGUMENT = 3,
    RSPT_STATUS_INVALID_CHANNEL_INDEX = 4,
    RSPT_STATUS_INVALID_PEAK_INDEX = 5,
    RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL = 6
} rspt_status;

typedef enum rspt_detection_mode {
    RSPT_MODE_DEFAULT = 0,
    RSPT_MODE_HIGH_SENSITIVITY = 1,
    RSPT_MODE_HIGH_PPV = 2
} rspt_detection_mode;

#define RSPT_VALID_R_PEAK_SAMPLE             (1ULL << 0)
#define RSPT_VALID_RR_INTERVAL_MS            (1ULL << 1)
#define RSPT_VALID_HEART_RATE_BPM            (1ULL << 2)
#define RSPT_VALID_P_WAVE_DURATION_MS        (1ULL << 3)
#define RSPT_VALID_PR_INTERVAL_MS            (1ULL << 4)
#define RSPT_VALID_PR_SEGMENT_MS             (1ULL << 5)
#define RSPT_VALID_QRS_DURATION_MS           (1ULL << 6)
#define RSPT_VALID_QT_INTERVAL_MS            (1ULL << 7)
#define RSPT_VALID_QTC_BAZETT_MS             (1ULL << 8)
#define RSPT_VALID_ST_SEGMENT_MS             (1ULL << 9)
#define RSPT_VALID_T_WAVE_DURATION_MS        (1ULL << 10)
#define RSPT_VALID_PQRS_T_ANNOTATION         (1ULL << 11)

typedef struct rspt_pqrst_annotation {
    int32_t p1_onset_sample;
    int32_t p1_peak_sample;
    int32_t p1_offset_sample;
    int32_t p2_onset_sample;
    int32_t p2_peak_sample;
    int32_t p2_offset_sample;

    int32_t qrs_onset_sample;
    int32_t r_peak_sample;
    int32_t qrs_offset_sample;
    int32_t q_peak_sample;
    int32_t s_peak_sample;

    int32_t t_onset_sample;
    int32_t t_peak_sample;
    int32_t t_offset_sample;
    int32_t j_point_sample;
} rspt_pqrst_annotation;

typedef struct rspt_ecg_beat_result {
    int32_t status;
    uint64_t valid_fields;
    uint32_t analysis_channel_index;
    uint32_t beat_index;
    uint32_t r_peak_sample;
    char status_message[64];

    rspt_pqrst_annotation annotation;

    double rr_interval_ms;
    double heart_rate_bpm;

    double p_wave_duration_ms;
    double p1_wave_duration_ms;
    double p1_amplitude_input_units;
    double p2_wave_duration_ms;
    double p2_amplitude_input_units;
    double pr_interval_ms;
    double pr_segment_ms;

    double q_duration_ms;
    double q_amplitude_input_units;
    double r_duration_ms;
    double r_amplitude_input_units;
    double s_duration_ms;
    double s_amplitude_input_units;
    double qrs_duration_ms;
    double qt_interval_ms;
    double qtc_bazett_ms;
    double st_segment_ms;
    double t_wave_duration_ms;

    double j_point_amplitude_input_units;
    double st20_amplitude_input_units;
    double st40_amplitude_input_units;
    double st60_amplitude_input_units;
    double st80_amplitude_input_units;
    double t_amplitude_input_units;
} rspt_ecg_beat_result;

typedef struct rspt_ecg_summary_result {
    int32_t status;
    uint64_t valid_fields;
    uint32_t analysis_channel_index;
    size_t r_peak_count;
    size_t analysed_beat_count;
    char status_message[64];

    double mean_rr_interval_ms;
    double rr_variation_ms;
    double heart_rate_bpm;

    double mean_p_wave_duration_ms;
    double mean_pr_interval_ms;
    double mean_qrs_duration_ms;
    double mean_qt_interval_ms;
    double mean_qtc_bazett_ms;
    double mean_st_segment_ms;
    double mean_t_wave_duration_ms;

    int32_t is_sinus_rhythm;
    int32_t premature_beat_count;
} rspt_ecg_summary_result;

RSPT_API const char* rspt_status_message(int32_t status);
RSPT_API uint32_t rspt_c_api_version(void);

RSPT_API int32_t rspt_detect_peaks_double(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    uint32_t* out_r_peak_indexes,
    size_t r_peak_capacity,
    size_t* out_r_peak_count);

RSPT_API int32_t rspt_analyse_ecg_beats_double(
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

RSPT_API int32_t rspt_analyse_ecg_summary_double(
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
