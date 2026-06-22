#pragma once

#include <stddef.h>
#include <stdint.h>

typedef enum rspt_status {
    /** Operation completed successfully. */
    RSPT_STATUS_OK = 0,

    /** No R peaks were available or detected. */
    RSPT_STATUS_NO_R_PEAKS = 1,

    /** Input channel array or selected channel data is missing. */
    RSPT_STATUS_NO_CHANNEL_DATA = 2,

    /** One or more scalar arguments are invalid. */
    RSPT_STATUS_INVALID_ARGUMENT = 3,

    /** Requested analysis channel index is outside the input channel range. */
    RSPT_STATUS_INVALID_CHANNEL_INDEX = 4,

    /** Supplied R-peak or beat index is outside the valid sample/peak range. */
    RSPT_STATUS_INVALID_PEAK_INDEX = 5,

    /** Caller-provided output buffer is smaller than the required output size. */
    RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL = 6
} rspt_status;

/**
 * R-peak detector operating mode.
 */
typedef enum rspt_detection_mode {
    /** Default detector mode used for CTS/CSE validation. */
    RSPT_MODE_DEFAULT = 0,

    /** More sensitive detector mode; may produce more candidate peaks. */
    RSPT_MODE_HIGH_SENSITIVITY = 1,

    /** Higher positive predictive value mode; may reject borderline peaks. */
    RSPT_MODE_HIGH_PPV = 2
} rspt_detection_mode;

/**
 * PQRST annotation sample indexes.
 *
 * All values are zero-based sample indexes in the input ECG signal.
 * Unavailable or not detected sample indexes are reported as -1.
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
 * when unavailable. Durations and intervals are milliseconds. Amplitudes are
 * expressed in the input signal units. status describes whether this beat was
 * analyzed successfully; use rspt_status_message/status_message to convert it
 * to text.
 */
typedef struct rspt_ecg_beat_result {
    /** rspt_status value for this beat. */
    int32_t status;

    /** Zero-based channel used for morphology analysis, or -1 if unavailable. */
    int32_t analysis_channel_index;

    /** Zero-based index of this beat in the effective R-peak list. */
    uint32_t beat_index;

    /** R-peak sample index for this beat, or -1 if unavailable. */
    int32_t r_peak_sample;

    /** PQRST sample indexes for this beat. */
    rspt_pqrst_annotation annotation;

    /** RR interval ending at this beat, in milliseconds. NaN if unavailable. */
    double rr_interval_ms;

    /** Instantaneous heart rate derived from rr_interval_ms, in bpm. NaN if unavailable. */
    double heart_rate_bpm;

    /** P-wave duration and morphology metrics. */
    double p_wave_duration_ms;
    double p1_wave_duration_ms;
    double p1_amplitude_input_units;
    double p2_wave_duration_ms;
    double p2_amplitude_input_units;
    double pr_interval_ms;
    double pr_segment_ms;
    double pp_interval_ms;
    double p_peak_to_end_interval_ms;

    /** QRS, QT, corrected QT, QT dispersion, ST, and T-wave metrics. */
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

    /** J-point, ST level, and T-wave amplitude metrics. */
    double j_point_amplitude_input_units;
    double st20_amplitude_input_units;
    double st40_amplitude_input_units;
    double st60_amplitude_input_units;
    double st80_amplitude_input_units;
    double t_amplitude_input_units;
} rspt_ecg_beat_result;

/**
 * Summary statistics for one ECG metric.
 *
 * valid_count is the number of valid, non-NaN beat values used to calculate
 * mean and standard_deviation. If valid_count is 0, mean and standard_deviation
 * are NaN. If valid_count is 1, mean is that single value and
 * standard_deviation is 0.
 */
typedef struct rspt_metric_statistics {
    /** Number of valid values included in the aggregate calculation. */
    size_t valid_count;

    /** Arithmetic mean of the valid values, or NaN if valid_count is 0. */
    double mean;

    /** Population standard deviation of the valid values, or NaN if valid_count is 0. */
    double standard_deviation;
} rspt_metric_statistics;

/**
 * Aggregate ECG analysis result.
 *
 * Double metrics use NaN when unavailable. Durations and intervals are
 * milliseconds. Amplitudes are expressed in the input signal units. status
 * describes whether the aggregate analysis completed successfully; use
 * rspt_status_message/status_message to convert it to text.
 */
typedef struct rspt_ecg_summary_result {
    /** rspt_status value for the aggregate analysis call. */
    int32_t status;

    /** Zero-based channel used for morphology analysis, or -1 if unavailable. */
    int32_t analysis_channel_index;

    /** Number of effective R peaks used by the analysis call. */
    size_t r_peak_count;

    /** Number of beat result entries produced in the analysis call. */
    size_t analysed_beat_count;

    /** Summary statistics calculated from per-beat RR interval values. */
    rspt_metric_statistics rr_interval_ms;

    /** max(RR) - min(RR), in milliseconds. NaN if unavailable. */
    double rr_variation_ms;

    /** Summary statistics calculated from per-beat heart-rate values. */
    rspt_metric_statistics heart_rate_bpm;

    /** Summary statistics calculated from per-beat morphology and interval metrics. */
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

    /** 1 if rhythm is classified as sinus, 0 if not, -1 if unavailable. */
    int32_t is_sinus_rhythm;

    /** Number of premature beats detected in the analyzed rhythm, or -1 if unavailable. */
    int32_t premature_beat_count;
} rspt_ecg_summary_result;
