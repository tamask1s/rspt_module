#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct pqrst_indxes
{
    int32_t p[6]; // p1_on, p1_peak, p1_off, p2_on, p2_peak, p2_off; p2 values are -1 when absent.
    int32_t r[5]; // qrs_on, r_peak, qrs_off, q_peak, s_peak.
    int32_t t[4]; // t_on, t_peak, t_off, j_point.
};

struct ecg_analysis_result
{
    int analysis_status;
    char status_message[64];

    double rr_interval_ms;
    double rr_variation_ms;
    double heart_rate_bpm;
    int is_sinus_rhythm;
    int premature_beat_count;

    double p_wave_duration_ms;
    double pr_interval_ms;
    double pr_segment_ms;
    double qrs_duration_ms;
    double qt_interval_ms;
    double qtc_interval_ms;
    double st_segment_ms;
    double t_wave_duration_ms;

    double p1_wave_duration_ms;
    double p1_amplitude_input_units;
    double p2_wave_duration_ms;
    double p2_amplitude_input_units;

    double q_duration_ms;
    double q_amplitude_input_units;
    double r_duration_ms;
    double r_amplitude_input_units;
    double s_duration_ms;
    double s_amplitude_input_units;

    double j_point_amplitude_input_units;
    double st20_amplitude_input_units;
    double st40_amplitude_input_units;
    double st60_amplitude_input_units;
    double st80_amplitude_input_units;
    double t_amplitude_input_units;

    double r_peak_amplitude_input_units[12];
    double s_wave_amplitude_input_units[12];
    double st_elevation_input_units[12];
    double st_depression_input_units[12];

    double frontal_plane_axis_deg;
    double horizontal_plane_axis_deg;
};

ecg_analysis_result analyse_ecg_detect_peaks(
    const double** data,
    size_t nr_channels,
    size_t nr_samples_per_channel,
    double sampling_rate,
    std::vector<pqrst_indxes>& annotations,
    std::vector<unsigned int>* peak_indexes = nullptr,
    std::string mode = "default",
    int analysis_ch_indx = -1,
    int analysis_peak_indx = 0);

void analyse_ecg_all_beats(
    const double** data,
    size_t nr_channels,
    size_t nr_samples_per_channel,
    double sampling_rate,
    const std::vector<unsigned int>& peak_indexes,
    std::vector<pqrst_indxes>& annotations,
    std::vector<ecg_analysis_result>& results,
    int analysis_ch_indx = -1);
