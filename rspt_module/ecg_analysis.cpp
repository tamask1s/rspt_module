#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <numeric>

using namespace std;

#include "ecg_analysis.h"

template <typename Compare>
static int find_extreme_from_to(const double* signal, unsigned int from, unsigned int to, unsigned int length, Compare comp)
{
    if (from >= length)
        return -1;

    unsigned int start = from;
    unsigned int end = std::min<unsigned int>(length - 1, std::max<unsigned int>(from, to));

    unsigned int idx_extreme = start;
    double val_extreme = signal[start];
    for (unsigned int i = start + 1; i <= end; ++i)
        if (comp(signal[i], val_extreme))
        {
            val_extreme = signal[i];
            idx_extreme = i;
        }
    return idx_extreme;
}

static int find_min_from_to(const double* signal, unsigned int from, unsigned int to, unsigned int length)
{
    return find_extreme_from_to(signal, from, to, length, std::less<double>());
}

static int find_max_from_to(const double* signal, unsigned int from, unsigned int to, unsigned int length)
{
    return find_extreme_from_to(signal, from, to, length, std::greater<double>());
}

struct pqrst_positions
{
    unsigned int p_on_idx;
    unsigned int p_idx;
    unsigned int p_off_idx;
    unsigned int q_idx;
    unsigned int r_idx;
    unsigned int s_idx;
    unsigned int t_on_idx;
    unsigned int t_idx;
    unsigned int t_off_idx;
};

static pqrst_positions detect_pqrst_positions(const double* lead2, unsigned int r_idx, double sampling_rate, unsigned int nr_samples)
{
    pqrst_positions pos{};
    pos.r_idx = r_idx;

    // Q
    unsigned int window_q = 0.04 * sampling_rate;
    unsigned int start_q = (r_idx > window_q) ? (r_idx - window_q) : 0;
    pos.q_idx = find_min_from_to(lead2, start_q, r_idx, nr_samples);

    // S
    unsigned int window_s = 0.04 * sampling_rate;
    unsigned int end_s = std::min(r_idx + window_s, nr_samples - 1);
    pos.s_idx = find_min_from_to(lead2, r_idx, end_s, nr_samples);

    // T
    unsigned int t_search_start = r_idx + 0.15 * sampling_rate;
    unsigned int t_search_end = std::min(r_idx + 0.40 * sampling_rate, (double)(nr_samples - 1));
    pos.t_idx = find_max_from_to(lead2, t_search_start + 1, t_search_end, nr_samples);
    pos.t_on_idx = find_min_from_to(lead2, t_search_start, pos.t_idx, nr_samples);

    // P
    unsigned int p_search_start = (r_idx > 0.25 * sampling_rate) ? (r_idx - 0.25 * sampling_rate) : 0;
    unsigned int p_search_end = (r_idx > 0.10 * sampling_rate) ? (r_idx - 0.10 * sampling_rate) : 0;
    pos.p_idx = find_max_from_to(lead2, p_search_start + 1, p_search_end, nr_samples);
    double p_val = lead2[pos.p_idx];
    double half_p_val = p_val * 0.15;

    pos.p_on_idx = pos.p_idx;
    for (int i = pos.p_idx; i >= (int)p_search_start; --i)
        if (lead2[i] < half_p_val)
        {
            pos.p_on_idx = i;
            break;
        }

    pos.p_off_idx = pos.p_idx;
    unsigned int p_off_limit = std::min(pos.p_idx + 0.25 * sampling_rate, (double)r_idx);
    for (unsigned int i = pos.p_idx + 1; i <= p_off_limit; ++i)
        if (lead2[i] < half_p_val)
        {
            pos.p_off_idx = i;
            break;
        }

    // T off
    unsigned int t_offset_end = std::min(pos.t_idx + 0.20 * sampling_rate, (double)(nr_samples - 1));
    pos.t_off_idx = pos.t_idx;
    for (unsigned int i = pos.t_idx + 1; i <= t_offset_end; ++i)
    {
        if (std::fabs(lead2[i] - 0.0) < 0.05 * std::fabs(lead2[pos.t_idx]))
        {
            pos.t_off_idx = i;
            break;
        }
    }

    return pos;
}

static void fill_analysis_result(ecg_analysis_result& result, const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples, double sampling_rate, const pqrst_positions& pos)
{
    result.pr_interval_ms = ((pos.q_idx - pos.p_on_idx) / sampling_rate) * 1000.0;
    result.p_wave_duration_ms = ((pos.p_off_idx - pos.p_on_idx) / sampling_rate) * 1000.0;
    result.qrs_duration_ms = ((pos.s_idx - pos.q_idx) / sampling_rate) * 1000.0;
    result.qt_interval_ms = ((pos.t_idx - pos.q_idx) / sampling_rate) * 1000.0;

    double rr_sec = result.rr_interval_ms / 1000.0;
    result.qtc_interval_ms = (rr_sec > 0.0) ? result.qt_interval_ms / std::sqrt(rr_sec) : 0.0;
    result.t_wave_duration_ms = ((pos.t_off_idx - pos.t_idx) / sampling_rate) * 1000.0;

    unsigned int st_point = pos.q_idx + 0.06 * sampling_rate;
    if (st_point >= nr_samples) st_point = nr_samples - 1;

    for (unsigned int ch = 0; ch < nr_ch && ch < 12; ++ch)
    {
        const double* lead = ecg_signal[ch];
        result.r_peak_amplitude_mV[ch] = lead[pos.r_idx];

        unsigned int window_s = 0.04 * sampling_rate;
        double s_val_ch = lead[pos.r_idx];
        unsigned int s_loc = pos.r_idx;
        unsigned int end_s = std::min(pos.r_idx + window_s, nr_samples - 1);
        for (unsigned int i = pos.r_idx; i <= end_s; ++i)
            if (lead[i] < s_val_ch)
            {
                s_val_ch = lead[i];
                s_loc = i;
            }
        result.s_wave_amplitude_mV[ch] = lead[s_loc];

        double pr_baseline = lead[pos.p_on_idx + (unsigned int)(0.05 * sampling_rate)];
        double st_val = lead[st_point] - pr_baseline;
        if (st_val >= 0.0)
        {
            result.st_elevation_mV[ch] = st_val;
            result.st_depression_mV[ch] = 0.0;
        }
        else
        {
            result.st_elevation_mV[ch] = 0.0;
            result.st_depression_mV[ch] = -st_val;
        }
    }

    for (unsigned int ch = nr_ch; ch < 12; ++ch)
    {
        result.r_peak_amplitude_mV[ch] = 0.0;
        result.s_wave_amplitude_mV[ch] = 0.0;
        result.st_elevation_mV[ch] = 0.0;
        result.st_depression_mV[ch] = 0.0;
    }

    double net_I = result.r_peak_amplitude_mV[0] - std::fabs(result.s_wave_amplitude_mV[0]);
    double net_aVF = result.r_peak_amplitude_mV[5] - std::fabs(result.s_wave_amplitude_mV[5]);
    result.frontal_plane_axis_deg = std::atan2(net_aVF, net_I) * 180.0 / M_PI;
    result.horizontal_plane_axis_deg = 0.0;
}

static void fill_annotations(vector<string>& annotations, const pqrst_positions& pos)
{
    annotations.clear();
    annotations.push_back(to_string(pos.p_on_idx) + ":(");
    annotations.push_back(to_string(pos.p_idx) + ":p");
    annotations.push_back(to_string(pos.p_off_idx) + ":)");

    annotations.push_back(to_string(pos.q_idx) + ":[");
    annotations.push_back(to_string(pos.r_idx) + ":N");
    annotations.push_back(to_string(pos.s_idx) + ":]");

    annotations.push_back(to_string(pos.t_on_idx) + ":{");
    annotations.push_back(to_string(pos.t_idx) + ":t");
    annotations.push_back(to_string(pos.t_off_idx) + ":}");
}

static void check_sinus_rhythm(const vector<double>& rr_intervals, double rr_mean, double max_rr, double min_rr, double sampling_rate, ecg_analysis_result& result)
{
    result.premature_beat_count = 0;
    for (double v : rr_intervals)
        if (rr_mean > 0.0 && v < 0.80 * rr_mean)
            result.premature_beat_count++;

    result.is_sinus_rhythm = (result.premature_beat_count || (rr_mean > 0.0 && (max_rr - min_rr) / rr_mean > 0.10)) ? 0 : 1;
}

static void calculate_rr_statistics(const vector<unsigned int>& peak_indexes, double sampling_rate, ecg_analysis_result& result)
{
    vector<double> rr_intervals;
    for (size_t i = 1; i < peak_indexes.size(); ++i)
        rr_intervals.push_back((peak_indexes[i] - peak_indexes[i - 1]) / sampling_rate * 1000.0);

    double total_ms = 0.0, max_rr = 0.0, min_rr = 1e9;
    for (unsigned int rr_ms : rr_intervals)
    {
        total_ms += rr_ms;
        if (rr_ms < min_rr) min_rr = rr_ms;
        if (rr_ms > max_rr) max_rr = rr_ms;
    }

    result.rr_interval_ms = total_ms / rr_intervals.size();
    result.rr_variation_ms = max_rr - min_rr;
    //result.rr_interval_ms = ((last_idx - peak_indexes.front()) / sampling_rate / (peak_indexes.size() - 1)) * 1000.0;
    result.heart_rate_bpm = (result.rr_interval_ms > 0.0) ? 60000.0 / result.rr_interval_ms : 0.0;
    check_sinus_rhythm(rr_intervals, result.rr_interval_ms, max_rr, min_rr, sampling_rate, result);
}

void analyse_ecg_multichannel(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const vector<unsigned int>& peak_indexes, vector<string>& annotations, ecg_analysis_result& result)
{
    annotations.clear();
    result.analysis_status = 0;
    std::snprintf(result.status_message, sizeof(result.status_message), "OK");

    if (peak_indexes.size() < 2)
    {
        result.analysis_status = 1;
        std::snprintf(result.status_message, sizeof(result.status_message), "Nincs elég R-csúcs az RR számításhoz");
        return;
    }

    const unsigned int lead_for_timing = 1;
    if (lead_for_timing >= nr_ch)
    {
        result.analysis_status = 2;
        std::snprintf(result.status_message, sizeof(result.status_message), "Hiányzó II. elvezetés a feldolgozáshoz");
        return;
    }

    pqrst_positions pos = detect_pqrst_positions(ecg_signal[lead_for_timing], peak_indexes[peak_indexes.size() - 2], sampling_rate, nr_samples_per_ch);
    calculate_rr_statistics(peak_indexes, sampling_rate, result);
    fill_analysis_result(result, ecg_signal, nr_ch, nr_samples_per_ch, sampling_rate, pos);
    fill_annotations(annotations, pos);
}
