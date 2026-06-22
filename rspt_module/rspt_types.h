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

/**
 * PQRST annotation sample indexes. Unavailable sample indexes are reported as -1.
 */
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

/**
 * Per-beat ECG analysis result.
 *
 * Double metrics use NaN when the value is unavailable. Sample indexes use -1
 * when unavailable. status describes whether this beat was analyzed
 * successfully; use rspt_status_message/status_message to convert it to text.
 */
typedef struct rspt_ecg_beat_result {
    int32_t status;
    int32_t analysis_channel_index;
    uint32_t beat_index;
    int32_t r_peak_sample;

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

/**
 * Summary statistics for a metric. If count is 0, mean and standard_deviation
 * are NaN.
 */
typedef struct rspt_metric_statistics {
    size_t count;
    double mean;
    double standard_deviation;
} rspt_metric_statistics;

/**
 * Aggregate ECG analysis result.
 *
 * Double metrics use NaN when unavailable. status describes whether the
 * aggregate analysis completed successfully; use rspt_status_message/
 * status_message to convert it to text.
 */
typedef struct rspt_ecg_summary_result {
    int32_t status;
    int32_t analysis_channel_index;
    size_t r_peak_count;
    size_t analysed_beat_count;

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
