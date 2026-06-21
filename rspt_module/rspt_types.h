#pragma once

#include <stddef.h>
#include <stdint.h>

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
#define RSPT_VALID_PP_INTERVAL_MS            (1ULL << 12)
#define RSPT_VALID_P_PEAK_TO_END_INTERVAL_MS (1ULL << 13)
#define RSPT_VALID_T_PEAK_TO_END_INTERVAL_MS (1ULL << 14)
#define RSPT_VALID_QT_DISPERSION_MS          (1ULL << 15)

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
    double pp_interval_ms;
    double p_peak_to_end_interval_ms;

    double q_duration_ms;
    double q_amplitude_input_units;
    double r_duration_ms;
    double r_amplitude_input_units;
    double s_duration_ms;
    double s_amplitude_input_units;
    double qrs_duration_ms;
    double qt_interval_ms;
    double qtc_bazett_ms;
    double qt_dispersion_ms;
    double st_segment_ms;
    double t_wave_duration_ms;
    double t_peak_to_end_interval_ms;

    double j_point_amplitude_input_units;
    double st20_amplitude_input_units;
    double st40_amplitude_input_units;
    double st60_amplitude_input_units;
    double st80_amplitude_input_units;
    double t_amplitude_input_units;
} rspt_ecg_beat_result;

typedef struct rspt_metric_statistics {
    size_t count;
    double mean;
    double standard_deviation;
} rspt_metric_statistics;

typedef struct rspt_ecg_summary_result {
    int32_t status;
    uint64_t valid_fields;
    uint32_t analysis_channel_index;
    size_t r_peak_count;
    size_t analysed_beat_count;
    char status_message[64];

    rspt_metric_statistics rr_interval_ms;
    double rr_variation_ms;
    rspt_metric_statistics heart_rate_bpm;
    rspt_metric_statistics p_wave_duration_ms;
    rspt_metric_statistics pr_interval_ms;
    rspt_metric_statistics pp_interval_ms;
    rspt_metric_statistics p_peak_to_end_interval_ms;
    rspt_metric_statistics qrs_duration_ms;
    rspt_metric_statistics qt_interval_ms;
    rspt_metric_statistics qtc_bazett_ms;
    rspt_metric_statistics qt_dispersion_ms;
    rspt_metric_statistics st_segment_ms;
    rspt_metric_statistics t_wave_duration_ms;
    rspt_metric_statistics t_peak_to_end_interval_ms;

    int32_t is_sinus_rhythm;
    int32_t premature_beat_count;
} rspt_ecg_summary_result;
