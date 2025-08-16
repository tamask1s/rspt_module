#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string.h>

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
    int p_on_idx;
    int p_idx;
    int p_off_idx;
    int q_idx;
    int r_idx;
    int s_idx;
    int t_on_idx;
    int t_idx;
    int t_off_idx;
};

unsigned int find_min_with_baseline_correction(const double* signal, unsigned int start_idx, unsigned int end_idx)
{
    if (end_idx <= start_idx) return start_idx;

    double x1 = start_idx;
    double y1 = signal[start_idx];
    double x2 = end_idx;
    double y2 = signal[end_idx];

    double m = (y2 - y1) / (x2 - x1);   // meredekség
    double b = y1 - m * x1;             // y-metszéspont

    double min_val = std::numeric_limits<double>::max();
    unsigned int min_idx = start_idx;

    for (unsigned int i = start_idx; i <= end_idx; ++i)
    {
        double baseline = m * (double)(i) + b;
        double corrected = signal[i] - baseline;
        if (corrected < min_val)
        {
            min_val = corrected;
            min_idx = i;
        }
    }
    return min_idx;
}

unsigned int find_max_with_baseline_correction(const double* signal, unsigned int start_idx, unsigned int end_idx)
{
    if (end_idx <= start_idx) return start_idx;

    double x1 = start_idx;
    double y1 = signal[start_idx];
    double x2 = end_idx;
    double y2 = signal[end_idx];

    double m = (y2 - y1) / (x2 - x1);
    double b = y1 - m * x1;

    double max_val = std::numeric_limits<double>::lowest();
    unsigned int max_idx = start_idx;

    for (unsigned int i = start_idx; i <= end_idx; ++i)
    {
        double corrected = signal[i] - (m * (double)i + b);
        if (corrected > max_val)
        {
            max_val = corrected;
            max_idx = i;
        }
    }
    return max_idx;
}

unsigned int find_isoelectric_point_before(const double* signal, unsigned int start_idx, unsigned int peak_idx, double tolerance = 0.05)
{
    if (peak_idx <= start_idx) return start_idx;

    double peak_val = signal[peak_idx];

    for (int i = static_cast<int>(peak_idx); i >= static_cast<int>(start_idx); --i)
    {
        if (std::fabs(signal[i]) <= std::fabs(peak_val) * tolerance)
            return i;
    }

    return start_idx;
}

unsigned int find_isoelectric_point_after(const double* signal, unsigned int peak_idx, unsigned int end_idx, double tolerance = 0.05)
{
    if (end_idx <= peak_idx) return end_idx;

    double peak_val = signal[peak_idx];

    for (unsigned int i = peak_idx; i <= end_idx; ++i)
    {
        if (std::fabs(signal[i]) <= std::fabs(peak_val) * tolerance)
            return i;
    }

    return end_idx;
}

double compute_rms(const double* sig, size_t sig_size, size_t start, size_t end)
{
    if (end > sig_size) end = sig_size;
    if (start >= end) return 0.0;
    double sum_sq = 0.0;
    for (size_t i = start; i < end; ++i)
        sum_sq += sig[i] * sig[i];
    return std::sqrt(sum_sq / (end - start));
}

double max_abs_peak(const double* sig, size_t sig_size, size_t start, size_t end, int& peak_indx, int& sign)
{
    if (end > sig_size) end = sig_size;
    if (start >= end) return 0.0;
    double max_val = -std::numeric_limits<double>::infinity();
    double min_val = std::numeric_limits<double>::infinity();
    int max_indx = -1;
    int min_indx = -1;
    for (size_t i = start; i < end; ++i)
    {
        if (sig[i] > max_val)
        {
            max_val = sig[i];
            max_indx = i;
        }
        if (sig[i] < min_val)
        {
            min_val = sig[i];
            min_indx = i;
        }
    }
    if (max_val > -min_val)
    {
        peak_indx = max_indx;
        sign = 1;
        return max_val;
    }
    else
    {
        peak_indx = min_indx;
        sign = -1;
        return min_val;
    }
}

int xcor_align(const double** leads, size_t nr_samples, size_t ref_ch, size_t ch, size_t t_start, size_t t_end, size_t max_lag_samples)
{
    if (t_end > nr_samples) t_end = nr_samples;
    if (t_start >= t_end) return 0;

    const double* ref = leads[ref_ch];
    const double* sig = leads[ch];

    double best_corr = -std::numeric_limits<double>::infinity();
    int best_lag = 0;

    for (int lag = -(int)max_lag_samples; lag <= (int)max_lag_samples; ++lag)
    {
        size_t start_ref = t_start;
        size_t start_sig = t_start + lag;

        if (lag < 0)
        {
            start_ref = t_start - lag;
            start_sig = t_start;
        }

        size_t len = t_end - t_start;
        if (start_ref + len > nr_samples || start_sig + len > nr_samples)
            len = std::min(nr_samples - start_ref, nr_samples - start_sig);

        if (len <= 1) continue;

        double sum_ref = 0.0, sum_sig = 0.0;
        double sum_ref2 = 0.0, sum_sig2 = 0.0;
        double sum_cross = 0.0;

        for (size_t i = 0; i < len; ++i)
        {
            double xr = ref[start_ref + i];
            double xs = sig[start_sig + i];
            sum_ref += xr;
            sum_sig += xs;
            sum_ref2 += xr * xr;
            sum_sig2 += xs * xs;
            sum_cross += xr * xs;
        }

        double num = sum_cross - (sum_ref * sum_sig) / len;
        double den = std::sqrt((sum_ref2 - (sum_ref * sum_ref) / len) * (sum_sig2 - (sum_sig * sum_sig) / len));

        if (den > 1e-12)
        {
            double corr = num / den;
            if (corr > best_corr)
            {
                best_corr = corr;
                best_lag = lag;
            }
        }
    }

    return best_lag;
}

void accumulate_aligned(double* res, const double* src, size_t nr_src_samples, size_t t_start, size_t t_end, int dt, double w)
{
    if (t_end > nr_src_samples) t_end = nr_src_samples;
    size_t len = t_end - t_start;
    if (len == 0) return;

    for (size_t i = 0; i < len; ++i)
    {
        int si = (int)t_start + (int)i + dt;
        if (si >= 0 && si < (int)nr_src_samples)
        {
            res[i] += w * src[si];
        }
    }
}

struct LeadScore
{
    size_t ch;
    double snr;
    int sign;
    int dt;
}; // sign: +1/-1, dt: mintabeli eltolás

int create_ideal_signal(double* res, const double** leads, size_t nr_channels, size_t nr_samples, double sampling_rate, size_t t_start, size_t t_end, size_t iso_start, size_t iso_end)
{
    if (t_end > nr_samples) t_end = nr_samples;
    size_t len = t_end - t_start;
    if (len == 0) return -1;

    // 1) SNR becslés és rangsor
    std::vector<LeadScore> scores;
    scores.reserve(nr_channels);
    for (size_t ch = 0; ch < nr_channels; ++ch)
    {
        double sigma = compute_rms(leads[ch], nr_samples, iso_start, iso_end);
        int peak_indx;
        int sign;
        double peak = max_abs_peak(leads[ch], nr_samples, t_start, t_end, peak_indx, sign);
        double snr = peak / (sigma + 1e-9);
        scores.push_back({ch, snr, sign, 0});
    }

    std::sort(scores.begin(), scores.end(), [](auto &a, auto &b)
    {
        return a.snr > b.snr;
    });
    int K = std::min(3, (int)scores.size());

    // 2) időigazítás a legjobb csatornához
    int ref_ch = scores[0].ch;
    for (int k = 1; k < K; ++k)
    {
        int ch = scores[k].ch;
        int max_lag_samples = 0.015 * sampling_rate;
        scores[k].dt = xcor_align(leads, nr_samples, ref_ch, ch, t_start, t_end, max_lag_samples);
    }

    // 3) súlyozott összeg
    std::fill(res, res + len, 0.0);
    double wsum = 0.0;
    for (int k = 0; k < K; ++k)
        wsum += scores[k].snr;

    for (int k = 0; k < K; ++k)
    {
        double w = scores[k].snr / (wsum + 1e-12);
        cout << "scores[k].sign: " << scores[k].sign << endl;
        accumulate_aligned(res, leads[scores[k].ch], nr_samples, t_start, t_end, scores[k].dt, scores[k].sign * w);
    }

    return 0;
}

static pqrst_positions detect_pqrst_positions(const double** leads, unsigned int nr_ch, unsigned int r_idx, double sampling_rate, unsigned int nr_samples)
{
    pqrst_positions pos{};
    pos.r_idx = r_idx;

    // Q
    unsigned int window_q = 0.04 * sampling_rate;
    unsigned int start_q = (r_idx > window_q) ? (r_idx - window_q) : 0;
    pos.q_idx = find_min_from_to(leads[0], start_q, r_idx, nr_samples);

    // S
    unsigned int window_s = 0.04 * sampling_rate;
    unsigned int end_s = std::min(r_idx + window_s, nr_samples - 1);
    pos.s_idx = find_min_from_to(leads[0], r_idx, end_s, nr_samples);

    // T
    unsigned int t_search_start = r_idx + 0.15 * sampling_rate;
    unsigned int t_search_end = std::min(r_idx + 0.40 * sampling_rate, (double)(nr_samples - 1));

    size_t search_samples = t_search_end - t_search_start + 1; ///(size_t)sampling_rate * 2;
    double sig_window[search_samples] = {};
    create_ideal_signal(sig_window, leads, nr_ch, nr_samples, sampling_rate, t_search_start, t_search_end, t_search_start, t_search_end);

    pos.t_idx = find_max_from_to(sig_window, 1, search_samples, search_samples);
    pos.t_on_idx = find_min_with_baseline_correction(sig_window, 0, pos.t_idx);
    unsigned int t_offset_end = std::min(pos.t_idx + (int)(0.20 * sampling_rate), (int)search_samples - 1);
    pos.t_off_idx = find_min_with_baseline_correction(sig_window, pos.t_idx, t_offset_end);
    pos.t_idx += t_search_start;
    pos.t_on_idx += t_search_start;
    pos.t_off_idx += t_search_start;

    // P
    unsigned int p_search_start = (r_idx > 0.25 * sampling_rate) ? (r_idx - 0.25 * sampling_rate) : 0;
    unsigned int p_search_end = pos.q_idx;

    pos.p_idx = find_max_with_baseline_correction(leads[0], p_search_start + 1, p_search_end);
    pos.p_on_idx = find_min_with_baseline_correction(leads[0], p_search_start, pos.p_idx);
    pos.p_off_idx = find_min_with_baseline_correction(leads[0], pos.p_idx, pos.q_idx);

    return pos;
}

static void fill_analysis_result(ecg_analysis_result& result, const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples, double sampling_rate, pqrst_positions& pos)
{
    result.pr_interval_ms = ((pos.q_idx - pos.p_off_idx) / sampling_rate) * 1000.0;
    result.p_wave_duration_ms = ((pos.p_off_idx - pos.p_on_idx) / sampling_rate) * 1000.0;
    unsigned int window_isoelectric_search = 0.005 * sampling_rate;
    int qrs_start = find_isoelectric_point_before(ecg_signal[1], pos.q_idx - window_isoelectric_search, pos.q_idx);
    window_isoelectric_search = 0.02 * sampling_rate;
    int qrs_end = find_isoelectric_point_after(ecg_signal[1], pos.s_idx, pos.s_idx + window_isoelectric_search);
    pos.q_idx = qrs_start;
    pos.s_idx = qrs_end;
    result.qrs_duration_ms = ((qrs_end - qrs_start) / sampling_rate) * 1000.0;
    result.qt_interval_ms = ((pos.t_off_idx - pos.q_idx) / sampling_rate) * 1000.0;

    double rr_sec = result.rr_interval_ms / 1000.0;
    result.qtc_interval_ms = (rr_sec > 0.0) ? result.qt_interval_ms / std::sqrt(rr_sec) : 0.0;
    result.t_wave_duration_ms = ((pos.t_off_idx - pos.t_on_idx) / sampling_rate) * 1000.0;

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

static void fill_annotations(std::vector<pqrst_indxes>& annotations, const pqrst_positions& pos)
{
    annotations.clear();
    annotations.push_back(
    {
        {pos.p_on_idx, pos.p_idx, pos.p_off_idx},
        {pos.q_idx, pos.r_idx, pos.s_idx},
        {pos.t_on_idx, pos.t_idx, pos.t_off_idx}
    });
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
    //result.rr_interval_ms = (peak_indexes[peak_indexes.size() - 2] - peak_indexes[peak_indexes.size() - 3]) / sampling_rate * 1000.0;
    result.rr_variation_ms = max_rr - min_rr;
    //result.rr_interval_ms = ((last_idx - peak_indexes.front()) / sampling_rate / (peak_indexes.size() - 1)) * 1000.0;
    result.heart_rate_bpm = (result.rr_interval_ms > 0.0) ? 60000.0 / result.rr_interval_ms : 0.0;
    check_sinus_rhythm(rr_intervals, result.rr_interval_ms, max_rr, min_rr, sampling_rate, result);
}

void analyse_ecg_multichannel(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, ecg_analysis_result& result)
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

    pqrst_positions pos = detect_pqrst_positions(ecg_signal, nr_ch, peak_indexes[peak_indexes.size() - 2], sampling_rate, nr_samples_per_ch);
    calculate_rr_statistics(peak_indexes, sampling_rate, result);
    fill_analysis_result(result, ecg_signal, nr_ch, nr_samples_per_ch, sampling_rate, pos);
    fill_annotations(annotations, pos);
}
