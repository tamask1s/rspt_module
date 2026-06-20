#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <string.h>
#include <cstdint>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>

using namespace std;

#include "ecg_analysis.h"
#include "filter.h"
#include "iir_filter_opt.h"
#include "peak_detector.h"

namespace
{

constexpr size_t kStatusMessageSize = 64;

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
    double qtc_bazett_ms;
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
};

}

struct ecg_analysis_config
{
    struct artifact_t
    {
        double artifact_removal_ratio = 5.0;
        double global_std_multiplier = 4.0;
        double local_std_multiplier = 6.0;
        bool enable_artifact_removal = true;
    } artifact;

    struct
    {
        double l1norm_dt1_multiplier = 2.0;
        double l1norm_dt2_multiplier = 4.0;
        double r_wave_deriv_ms = 6.0;
    } deriv;

    struct
    {
        double qrs_isoel_tolerance = 0.035;
        double p_isoel_tolerance = 0.21;
        double t_isoel_tolerance = 0.25;
        double qrs_dt_ms = 2.0;
        double p_dt_ms = 2.0;
        double t_dt_ms = 4.0;
        double qrs_perimeter_ms = 25.0;
        double p_perimeter_ms = 30.0;
        double t_perimeter_ms = 35.0;
        double qrs_max_search_ms = 110.0;
        double p_max_search_ms = 200.0;
        double t_max_search_ms = 200.0;
        int qrs_extend_mode = 0;
        int p_extend_mode = 1;
        int t_extend_mode = 1;
    } isoel;

    struct
    {
        double wave_presence_ratio = 0.01;
        double wave_bound_ratio = 0.01;
        double r_wave_bound_ratio = 0.05;
        double qrs_shrink_trshld = 0.05;
        double qrs_shrink_bound_ms = 10.0;
        double isoel_r_search_ms = 30.0;
    } wave;

    struct
    {
        double min_wave_presence_ratio = 0.01;
    } phys;

    struct
    {
        double max_s_duration_sec = 0.05;
    } s_wave;

    struct
    {
        double pq_interval_offset_ms = 4.0;
        double qt_interval_offset_ms = 8.0;
    } timing;

    struct
    {
        int nr_filter_iterations = 2;
        double low_cut_hz = 0.1;
        double high_cut_hz = 40.0;
    } filter;

    struct
    {
        double p_isoel_offset_ms = 15.0;
        double p_search_left_ms = 250.0;
        double t_search_right_ms = 300.0;
        double biphasic_p_max_gap_ms = 20.0;
        double biphasic_p_ratio = 0.3;
    } search;

    struct
    {
        double amplitude_scale = 1.0;
    } bias;
};

void write_binmx_to_file(const char* filename, const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate)
{
    FILE* f = fopen(filename, "wb");
    if (!f)
    {
        //perror("fopen");
        return;
    }

    // --- HEADER --- 0–7: sampling rate
    fwrite(&sampling_rate, sizeof(double), 1, f);

    // 8–11: nr_channels (int32)
    int32_t ch = (int32_t)nr_channels;
    fwrite(&ch, sizeof(int32_t), 1, f);

    // 12–23: padding (zeros)
    uint8_t padding[12] = {0};
    fwrite(padding, sizeof(padding), 1, f);

    // --- DATA --- interleaved: ch0_s0, ch1_s0, ..., chN_s0, ch0_s1, ...
    for (size_t s = 0; s < nr_samples_per_channel; ++s)
    {
        for (size_t c = 0; c < nr_channels; ++c)
        {
            double value = data[c][s];
            fwrite(&value, sizeof(double), 1, f);
        }
    }

    fclose(f);
}

double median(const double* data, int len)
{
    if (len <= 0) return 0.0;

    std::vector<double> v(data, data + len);
    const int mid = len / 2;

    std::nth_element(v.begin(), v.begin() + mid, v.end());
    double m = v[mid];

    if (len % 2 == 0)
    {
        std::nth_element(v.begin(), v.begin() + mid - 1, v.begin() + mid);
        m = (m + v[mid - 1]) * 0.5;
    }

    return m;
}

void remove_artifacts(double* signal, int N, double fs, const ecg_analysis_config::artifact_t& c)
{
    if (!signal || N < 3)
        return;

    /* -------- 1. Első sample speciális kezelése -------- */
    double d01 = fabs(signal[0] - signal[1]);
    double d02 = fabs(signal[0] - signal[2]);
    double d12 = fabs(signal[1] - signal[2]);

    // Ha az első minta nagyon eltér, de a következők egymáshoz közel vannak
    if (d01 > c.artifact_removal_ratio * d12 && d02 > c.artifact_removal_ratio * d12)
        signal[0] = signal[1];

    /* -------- 2. Globális statisztika -------- */
    double mean = 0.0;
    for (int i = 0; i < N; i++)
        mean += signal[i];
    mean /= N;

    double var = 0.0;
    for (int i = 0; i < N; i++)
    {
        double d = signal[i] - mean;
        var += d * d;
    }
    double std = sqrt(var / N);

    if (std < 1e-6)
        return;

    const double GLOBAL_TH = c.global_std_multiplier * std;
    const double LOCAL_TH  = c.local_std_multiplier  * std;

    /* -------- 3. Többi artifact eltávolítása -------- */
    for (int i = 1; i < N - 1; i++)
    {
        int is_artifact = 0;

        if (fabs(signal[i] - mean) > GLOBAL_TH)
            is_artifact = 1;
        else if (fabs(signal[i] - signal[i - 1]) > LOCAL_TH && fabs(signal[i] - signal[i + 1]) > LOCAL_TH)
            is_artifact = 1;

        if (!is_artifact)
            continue;

        signal[i] = 0.5 * (signal[i - 1] + signal[i + 1]);
    }

    /* -------- 4. Utolsó sample -------- */
    if (fabs(signal[N - 1] - signal[N - 2]) > GLOBAL_TH)
        signal[N - 1] = signal[N - 2];
}

inline int clamp(int i, int len)
{
    if (i < 0) return 0;
    if (i >= len) return len - 1;
    return i;
}

/** Multi-scale absolute temporal difference operator. A rövid skála a QRS meredekségét, a hosszabb skála a T-hullám energiáját emeli ki. */
void abs_L1_norm(const double* signal, double* deriv, int len, double dt, double fs, double dt1_multiplier, double dt2_multiplier)
{
    std::fill(deriv, deriv + len, 0.0);
    int step = std::max(1, int(dt1_multiplier * fs / 1000.0));
    int nr_dt_samples = dt * fs / 1000.0;
    for (int i = 0; i < len; ++i)
    {
        for (int j = 0; j < nr_dt_samples; j += step)
            deriv[i] += fabs(signal[clamp(i + j, len)] - signal[clamp(i + j + nr_dt_samples, len)]);
        for (int j = 0; j < nr_dt_samples; j += step)
            deriv[i] += fabs(signal[clamp(i - j, len)] - signal[clamp(i - j - nr_dt_samples, len)]);
    }

    int nr_dt_samples_m = nr_dt_samples * dt2_multiplier;
    for (int i = 0; i < len; ++i)
    {
        for (int j = 0; j < nr_dt_samples_m; j += step)
            deriv[i] += fabs(signal[clamp(i + j, len)] - signal[clamp(i + j + nr_dt_samples_m, len)]);
        for (int j = 0; j < nr_dt_samples_m; j += step)
            deriv[i] += fabs(signal[clamp(i - j, len)] - signal[clamp(i - j - nr_dt_samples_m, len)]);
    }
}

inline void min_max(const double* arr, int len, double& min_val, double& max_val)
{
    min_val = 1e9;
    max_val = -1e9;
    for (int i = 0; i < len; ++i)
    {
        if (arr[i] > max_val) max_val = arr[i];
        if (arr[i] < min_val) min_val = arr[i];
    }
}

inline void min_max_idx(const double* arr, int len, double& min_val, double& max_val, int& min_indx, int& max_indx)
{
    min_val = 1e9;
    max_val = -1e9;
    for (int i = 0; i < len; ++i)
    {
        if (arr[i] > max_val)
        {
            max_val = arr[i];
            max_indx = i;
        }
        if (arr[i] < min_val)
        {
            min_val = arr[i];
            min_indx = i;
        }
    }
}

void search_isoel_bounds(double* signal, int len, double fs, int peak_indx, double max_search_ms, double isoel_ms, double perimeter_ms, int& isoel_l, int& isoel_r, const ecg_analysis_config& c, double isoel_tolerance, int extend_mode = 1, int a_stop_indx_r = -1, int* peak_p2 = 0, double* p_amp1 = 0, double* p_amp2 =0)
{
    if (!signal || len <= 0)
        return;

    int nr_max_search_samples = max_search_ms * fs / 1000.0;
    int nr_perimeter_samples = perimeter_ms * fs / 1000.0;

    double peak_indx_l = peak_indx - nr_perimeter_samples;
    if (peak_indx_l >= len) peak_indx_l = len - 1;
    if (peak_indx_l < 0) peak_indx_l = 0;
    int stop_indx_l = peak_indx - nr_max_search_samples;
    if (stop_indx_l < 0) stop_indx_l = 0;
    if (stop_indx_l >= len) stop_indx_l = len - 1;

    double peak_indx_r = peak_indx + nr_perimeter_samples;
    int stop_indx_r = a_stop_indx_r;
    if (stop_indx_r == -1)
        stop_indx_r = peak_indx + nr_max_search_samples;
    if (stop_indx_r < 0) stop_indx_r = 0;
    if (stop_indx_r >= len) stop_indx_r = len - 1;

    std::vector<double> deriv(len);
    abs_L1_norm(&signal[0], &deriv[0], len, isoel_ms, fs, c.deriv.l1norm_dt1_multiplier, c.deriv.l1norm_dt2_multiplier);

    double min_val, max_val;
    min_max(&deriv[0], len, min_val, max_val);

    double max_allowed_dev_r = min_val + (max_val - min_val) * isoel_tolerance * 0.1;

    if (extend_mode == 1) /** In case of P detection, climb down from the right (from the peak of the QRS complex, which scales nonlinearly) */
    {
        int deriv_samples = c.deriv.r_wave_deriv_ms * fs / 1000.0;
        if (stop_indx_r + deriv_samples >= len)
            stop_indx_r -= 1 + stop_indx_r + deriv_samples - len;
        for (int i = stop_indx_r; i >= peak_indx_r; --i)
            if (deriv[i] > max_allowed_dev_r && deriv[i] < deriv[i + deriv_samples])
                --stop_indx_r;
            else
                break;
        stop_indx_r -= deriv_samples * 2.0;
    }

    min_max(&deriv[stop_indx_l], stop_indx_r - stop_indx_l, min_val, max_val);
    double qrs_amp = max_val - min_val;
    //double deriv_baseline = median(&deriv[0], len);
    //double signal_baseline = median(&signal[0], len);
    double max_allowed_dev = /*deriv_baseline + */ min_val + qrs_amp * isoel_tolerance;

    vector<double> baselinev(len, max_allowed_dev);
    write_binmx_to_file("c:/Tamas/test002.bin", (const double**)&baselinev, 1, len, fs);

    /** climb down from the peak to left and right, skipping the middle nr_perimeter_samples * 2 */
    double max_allowed_dev_orig = max_allowed_dev;
    isoel_l = -1;
    int safetycount = 0;

repeat_l:
    for (int i = peak_indx_l; i >= stop_indx_l; --i)
    {
        if (deriv[i] < max_allowed_dev)
        {
            isoel_l = i;
            break;
        }
    }

    if (isoel_l == -1 && safetycount++ < 1000)
    {
        max_allowed_dev *= 1.01;
        goto repeat_l;
    }

    max_allowed_dev = max_allowed_dev_orig;
    safetycount = 0;
    isoel_r = -1;

repeat_r:
    for (int i = peak_indx_r; i < stop_indx_r; ++i)
    {
        if (deriv[i] < max_allowed_dev)
        {
            isoel_r = i;
            break;
        }
    }
    if (isoel_r == -1 && safetycount++ < 1000)
    {
        max_allowed_dev *= 1.01;
        goto repeat_r;
    }

    if (isoel_l == -1)
        isoel_l = std::max(0, std::min(len - 1, (int)(peak_indx - nr_perimeter_samples * 1.0)));
    if (isoel_r == -1)
        isoel_r = std::max(0, std::min(len - 1, (int)(peak_indx + nr_perimeter_samples * 1.0)));

    write_binmx_to_file("c:/Tamas/test003.bin", (const double**)&deriv, 1, len, fs);

    if (extend_mode == 0) /** case for QRS - shrink the complex */
    {
        min_max(&signal[isoel_l], isoel_r - isoel_l, min_val, max_val);
        double max_allowed_dev_l = std::max(fabs(min_val), max_val) * c.wave.qrs_shrink_trshld;

        //vector<double> baselinev(len, max_allowed_dev_l);
        //write_binmx_to_file("c:/Tamas/test002.bin", (const double**)&baselinev, 1, len, fs);

        int shrink_samples = c.wave.qrs_shrink_bound_ms * fs / 1000.0;
        double isoel_l_start = signal[isoel_l];
        peak_indx_l += shrink_samples;
        for (int i = isoel_l; i <= peak_indx_l; ++i)
            if (fabs(signal[i] - isoel_l_start) < max_allowed_dev_l)
                isoel_l = i;
            else
                break;

        double isoel_r_start = signal[isoel_r];
        peak_indx_r -= shrink_samples;
        for (int i = isoel_r; i >= peak_indx_r; --i)
            if (fabs(signal[i] - isoel_r_start) < max_allowed_dev_l)
                isoel_r = i;
            else
                break;
    }
}

void search_p_and_t_peaks(double* signal, int len, double fs, double max_search_ms_l, double max_search_ms_r, int isoel_l, int isoel_r, /*double max_allowed_dev, */int& peak_p1, int& peak_p2, int& peak_t, const ecg_analysis_config& c)
{
    int nr_max_search_samples_l = max_search_ms_l * fs / 1000.0;
    int searc_stop_l = isoel_l - nr_max_search_samples_l;
    if (searc_stop_l < 0) searc_stop_l = 0;
    if (isoel_l >= len) isoel_l = len - 1;
    peak_p1 = -1;
    peak_p2 = -1;
    double maxval = -1e9;
    double base = signal[isoel_l];
    if (isoel_l < 0) isoel_l = 0;
    for (int i = isoel_l; i >= searc_stop_l; --i)
        if (maxval < fabs(signal[i] - base))
        {
            maxval = fabs(signal[i] - base);
            peak_p1 = i;
        }

    int nr_max_search_samples_r = max_search_ms_r * fs / 1000.0;
    int searc_stop_r = nr_max_search_samples_r + isoel_r;
    if (searc_stop_r >= len) searc_stop_r = len - 1;
    peak_t = -1;
    int isoel_r_search_samples = c.wave.isoel_r_search_ms * fs / 1000.0;
    int base_indx = std::min(len - 1, isoel_r + isoel_r_search_samples);
    maxval = -1e9;
    base = signal[base_indx];
    for (int i = isoel_r; i < searc_stop_r; ++i)
        if (maxval < fabs(signal[i] - base))
        {
            maxval = fabs(signal[i] - base);
            peak_t = i;
        }
}

template<bool wave_is_positive>
inline void find_wave_bounds(const double* signal, int peak_idx, int left_limit, int right_limit, double baseline, double amplitude, int& lbound, int& rbound, const ecg_analysis_config& c)
{
    double trshld = amplitude * c.wave.wave_bound_ratio * (wave_is_positive ? 1.0 : -1.0);

    lbound = left_limit;
    rbound = right_limit;
    if constexpr (wave_is_positive)
    {
        for (int i = peak_idx; i > left_limit; --i)
            if (signal[i] - baseline <= trshld)
            {
                lbound = i;
                break;
            }
        for (int i = peak_idx; i < right_limit; ++i)
            if (signal[i] - baseline <= trshld)
            {
                rbound = i;
                break;
            }
    }
    else
    {
        for (int i = peak_idx; i > left_limit; --i)
            if (baseline - signal[i] <= trshld)
            {
                lbound = i;
                break;
            }
        for (int i = peak_idx; i < right_limit; ++i)
            if (baseline - signal[i] <= trshld)
            {
                rbound = i;
                break;
            }
    }
}

void fill_wave_measurements(double* signal, int len, pqrst_indxes& a, double fs, const ecg_analysis_config& c, ecg_analysis_result& result)
{
    double baseline = median(signal, len);
    vector<double> baselinev(len, baseline);
    write_binmx_to_file("c:/Tamas/test004.bin", (const double**)&baselinev, 1, len, fs);

    /* ================= P WAVES ================= */
    int overall_p_end = a.p[2];
    // P1
    if (a.p[0] >= 0 && a.p[1] >= 0 && a.p[2] >= 0 && a.p[0] < a.p[1] && a.p[1] < a.p[2] && a.p[2] < len)
    {
        result.p1_wave_duration_ms = (a.p[2] - a.p[0]) * 1000.0 / fs;
        result.p_wave_duration_ms = result.p1_wave_duration_ms;
        result.p1_amplitude_input_units = signal[a.p[1]] - baseline;
    }

    // P2 (optional)
    if (a.p[3] >= 0 && a.p[4] >= 0 && a.p[5] >= 0 && a.p[3] < a.p[4] && a.p[4] < a.p[5] && a.p[5] < len)
    {
        result.p2_wave_duration_ms = (a.p[5] - a.p[3]) * 1000.0 / fs;
        result.p_wave_duration_ms = result.p1_wave_duration_ms + result.p2_wave_duration_ms;
        result.p2_amplitude_input_units = signal[a.p[4]] - baseline;
        overall_p_end = a.p[5];
    }
    if (a.p[0] >= 0 && a.p[1] >= 0 && a.r[0] > a.p[0])
        result.pr_interval_ms = (a.r[0] - a.p[0]) * 1000.0 / fs - c.timing.pq_interval_offset_ms;

    /* ================= QRS ================= */
    if (a.r[0] >= 0 && a.r[1] >= 0 && a.r[2] >= 0 && a.r[0] < a.r[1] && a.r[1] < a.r[2] && a.r[2] < len)
    {
        result.qrs_duration_ms = (a.r[2] - a.r[0]) * 1000.0 / fs;
        result.r_amplitude_input_units = signal[a.r[1]] - baseline;
        bool no_r = false;
        double q_amp_if_no_r = 0;
        if (result.r_amplitude_input_units < 0)
        {
            double min_val, max_val, min_val2, max_val2;
            int min_indx = 0, max_indx = 0, min_indx2 = 0, max_indx2 = 0;
            min_max_idx(&signal[overall_p_end], a.r[1] - overall_p_end, min_val, max_val, min_indx, max_indx);
            min_max_idx(&signal[overall_p_end], a.r[2] - overall_p_end, min_val2, max_val2, min_indx2, max_indx2);
            result.r_amplitude_input_units = max_val - signal[overall_p_end];
            double r_prime_amp = max_val2 - signal[overall_p_end];
            a.r[1] = max_indx + overall_p_end;
            if (result.r_amplitude_input_units < r_prime_amp)
                result.r_amplitude_input_units = r_prime_amp;

            if ((min_val != 0) && (max_val / fabs(min_val) < c.wave.r_wave_bound_ratio))
            {
                no_r = true;
                a.r[1] = min_indx;
                q_amp_if_no_r = min_val;
            }
        }

        /* ---- Q wave ---- */
        if (no_r)
        {
            result.q_duration_ms = result.qrs_duration_ms;
            result.q_amplitude_input_units = q_amp_if_no_r - baseline;
            a.r[3] = a.r[1];
        }
        else
        {
            int q_idx = -1;
            double q_ext = signal[a.r[1]];

            for (int i = a.r[0]; i < a.r[1]; ++i)
                if (signal[i] < q_ext)
                {
                    q_ext = signal[i];
                    q_idx = i;
                }

            if (q_idx >= 0)
            {
                result.q_duration_ms = (q_idx - a.r[0]) * 1000.0 / fs;
                result.q_amplitude_input_units = q_ext - baseline;

                if (-result.q_amplitude_input_units >= result.r_amplitude_input_units * c.wave.wave_presence_ratio)
                {
                    int lbound, rbound;
                    find_wave_bounds<false>(signal, q_idx, a.r[0], a.r[2], baseline, result.q_amplitude_input_units, lbound, rbound, c);
                    result.q_duration_ms = (rbound - lbound) * 1000.0 / fs;
                    a.r[3] = q_idx;
                }
                else
                {
                    result.q_duration_ms = 0;
                    result.q_amplitude_input_units = 0;
                    a.r[3] = -1;
                }
            }
        }

        /* ---- S wave ---- */
        int s_idx = -1;
        double s_ext = signal[a.r[1]];
        int s_end = std::min(a.r[2], a.r[1] + int(c.s_wave.max_s_duration_sec * fs));

        for (int i = a.r[1]; i < s_end; ++i)
            if (signal[i] < s_ext)
            {
                s_ext = signal[i];
                s_idx = i;
            }

        if (s_idx >= 0)
        {
            result.s_duration_ms = (s_idx - a.r[1]) * 1000.0 / fs;
            result.s_amplitude_input_units = s_ext - baseline;

            if (-result.s_amplitude_input_units >= result.r_amplitude_input_units * c.wave.wave_presence_ratio)
            {
                int lbound, rbound;
                find_wave_bounds<false>(signal, s_idx, a.r[1], a.r[2], baseline, result.s_amplitude_input_units, lbound, rbound, c);
                result.s_duration_ms = (rbound - lbound) * 1000.0 / fs;

                if (fabs(signal[rbound] - s_ext) < fabs(result.s_amplitude_input_units) * c.phys.min_wave_presence_ratio) /// there is no significant rising right edge of S -> no S
                {
                    result.s_duration_ms = 0;
                    result.s_amplitude_input_units = 0;
                    a.r[4] = -1;
                }
                else
                    a.r[4] = s_idx;
            }
            else
            {
                result.s_duration_ms = 0;
                result.s_amplitude_input_units = 0;
                a.r[4] = -1;
            }
        }

        /* ---- R wave duration ---- */
        if (no_r)
        {
            result.r_duration_ms = 0;
            result.r_amplitude_input_units = 0;
        }
        else if (result.r_amplitude_input_units != 0)
        {
            int lbound, rbound;
            find_wave_bounds<true>(signal, a.r[1], a.r[0], a.r[2], baseline, result.r_amplitude_input_units, lbound, rbound, c);
            result.r_duration_ms = (rbound - lbound) * 1000.0 / fs;
        }
        else
            result.r_duration_ms = std::max(0.0, result.qrs_duration_ms - result.q_duration_ms - result.s_duration_ms);
    }

    if (a.t[1] >= 0)
        result.qt_interval_ms = (a.t[2] - a.r[0]) * 1000.0 / fs + c.timing.qt_interval_offset_ms;

    /* ================= ST / J ================= */
    int j = a.r[2];
    if (j >= 0 && j < len)
    {
        result.j_point_amplitude_input_units = signal[j] - baseline;
        a.t[3] = j;

        int s20 = 20.0 * fs / 1000.0;
        if (j + s20 < len) result.st20_amplitude_input_units = signal[j + s20] - baseline;
        if (j + s20 * 2 < len) result.st40_amplitude_input_units = signal[j + s20 * 2] - baseline;
        if (j + s20 * 3 < len) result.st60_amplitude_input_units = signal[j + s20 * 3] - baseline;
        if (j + s20 * 4 < len) result.st80_amplitude_input_units = signal[j + s20 * 4] - baseline;
    }

    /* ================= T ================= */
    if (a.t[1] >= 0 && a.t[1] < len)
        result.t_amplitude_input_units = signal[a.t[1]] - baseline;
}

void analyse_ecg(double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, ecg_analysis_result& result, const ecg_analysis_config& c, unsigned int analysis_ch_indx = 0, int analysis_peak_indx = 0)
{
    annotations.clear();
    memset(&result, 0, sizeof(result));
    result.analysis_status = 0;
    std::snprintf(result.status_message, sizeof(result.status_message), "OK");

    if (peak_indexes.size() < 1)
    {
        result.analysis_status = 1;
        std::snprintf(result.status_message, sizeof(result.status_message), "No R peaks found");
        cout << "!!! WARNING !!!: No R peaks found" << endl;
        return;
    }

    if (!ecg_signal || nr_ch < 1)
    {
        result.analysis_status = 2;
        std::snprintf(result.status_message, sizeof(result.status_message), "No channel data given");
        return;
    }

    if (nr_ch <= analysis_ch_indx)
    {
        result.analysis_status = 4;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid channel index");
        return;
    }

    if (analysis_peak_indx < 0 || analysis_peak_indx >= static_cast<int>(peak_indexes.size()))
    {
        result.analysis_status = 5;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid analysis peak index");
        return;
    }

    double* lead_cpy = new double[nr_samples_per_ch];
    for (unsigned int i = 0; i < nr_samples_per_ch; i++)
    {
        ecg_signal[analysis_ch_indx][i] *= c.bias.amplitude_scale;
        lead_cpy[i] = ecg_signal[analysis_ch_indx][i];
    }

    if (c.filter.nr_filter_iterations)
    {
        iir_filter_4th_order bandpass_filter_;
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 2, sampling_rate, c.filter.low_cut_hz, c.filter.high_cut_hz);

        auto forward_backward_filter = [&](double* data, unsigned int length, unsigned int iterations)
        {
            for (unsigned int it = 0; it < iterations; ++it)
            {
                bandpass_filter_.init_history_values(data[0], sampling_rate);
                for (unsigned int i = 0; i < length; ++i)
                    data[i] = bandpass_filter_.filter(data[i]);
                bandpass_filter_.init_history_values(data[length - 1], sampling_rate);
                for (int i = static_cast<int>(length) - 1; i >= 0; --i)
                    data[i] = bandpass_filter_.filter(data[i]);
            }
        };

        forward_backward_filter(lead_cpy, nr_samples_per_ch, c.filter.nr_filter_iterations);
    }

    if (c.artifact.enable_artifact_removal)
        remove_artifacts(lead_cpy, nr_samples_per_ch, sampling_rate, c.artifact);
    write_binmx_to_file("c:/Tamas/test001.bin", (const double**)&lead_cpy, 1, nr_samples_per_ch, sampling_rate);

    int isoel_start = -1, isoel_r = -1;
    search_isoel_bounds(lead_cpy, nr_samples_per_ch, sampling_rate, peak_indexes[analysis_peak_indx], c.isoel.qrs_max_search_ms, c.isoel.qrs_dt_ms, c.isoel.qrs_perimeter_ms, isoel_start, isoel_r, c, c.isoel.qrs_isoel_tolerance, c.isoel.qrs_extend_mode);
    if (isoel_start < 0 || isoel_r < 0 || isoel_start >= (int)nr_samples_per_ch || isoel_r >= (int)nr_samples_per_ch || isoel_start >= isoel_r)
    {
        result.analysis_status = 5;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid QRS boundaries");
        delete[] lead_cpy;
        return;
    }

    int peak_p1, peak_p2 = -1, peak_t;
    int p_isoel_offset_samples = c.search.p_isoel_offset_ms * sampling_rate / 1000.0;
    search_p_and_t_peaks(lead_cpy, nr_samples_per_ch, sampling_rate, c.search.p_search_left_ms, c.search.t_search_right_ms, isoel_start - p_isoel_offset_samples, isoel_r, peak_p1, peak_p2, peak_t, c);
    if (peak_p1 < 0 || peak_p1 >= (int)nr_samples_per_ch || peak_t < 0 || peak_t >= (int)nr_samples_per_ch)
    {
        result.analysis_status = 5;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid P or T peak");
        delete[] lead_cpy;
        return;
    }

    int isoel_p_l = -1, isoel_p_r = -1;
    double p_amp1 = 0, p_amp2 = 0;
    search_isoel_bounds(lead_cpy, nr_samples_per_ch, sampling_rate, peak_p1, c.isoel.p_max_search_ms, c.isoel.p_dt_ms, c.isoel.p_perimeter_ms, isoel_p_l, isoel_p_r, c, c.isoel.p_isoel_tolerance, c.isoel.p_extend_mode, isoel_start, &peak_p2, &p_amp1, &p_amp2);
    if (isoel_p_l < 0 || isoel_p_r < 0 || isoel_p_l >= (int)nr_samples_per_ch || isoel_p_r >= (int)nr_samples_per_ch || isoel_p_l >= isoel_p_r)
    {
        result.analysis_status = 5;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid P wave boundaries");
        delete[] lead_cpy;
        return;
    }

    int p2_indx = -1;
    double min_val, max_val;
    int min_indx = 0, max_indx = 0;
    double baseline = (lead_cpy[isoel_p_l] + lead_cpy[isoel_p_r]) / 2;
    min_max_idx(lead_cpy + isoel_p_l, isoel_p_r - isoel_p_l, min_val, max_val, min_indx, max_indx);
    double pos_amp = max_val - baseline;
    double neg_amp = baseline - min_val;
    min_indx += isoel_p_l;
    max_indx += isoel_p_l;

    if (pos_amp > 0 && neg_amp > 0) /** Handle biphasic P */
    {
        double smaller = std::min(pos_amp, neg_amp);
        double larger  = std::max(pos_amp, neg_amp);
        int tol_samples = c.search.biphasic_p_max_gap_ms * sampling_rate / 1000.0;
        if (smaller >= c.search.biphasic_p_ratio * larger)
        {
            if (larger == pos_amp)
            {
                if (abs(peak_p1 - max_indx) < tol_samples)
                    p2_indx = min_indx;
            }
            else
            {
                if (abs(peak_p1 - min_indx) < tol_samples)
                    p2_indx = max_indx;
            }
        }
        if (p2_indx < peak_p1)
            p2_indx = -1;
    }

    int isoel_t_l = -1, isoel_t_r = -1;
    search_isoel_bounds(lead_cpy, nr_samples_per_ch, sampling_rate, peak_t, c.isoel.t_max_search_ms, c.isoel.t_dt_ms, c.isoel.t_perimeter_ms, isoel_t_l, isoel_t_r, c, c.isoel.t_isoel_tolerance, c.isoel.t_extend_mode);
    if (isoel_t_l < 0 || isoel_t_r < 0 || isoel_t_l >= (int)nr_samples_per_ch || isoel_t_r >= (int)nr_samples_per_ch || isoel_t_l >= isoel_t_r)
    {
        result.analysis_status = 5;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid T wave boundaries");
        delete[] lead_cpy;
        return;
    }

    annotations.push_back({});
    for (int i = 0; i < 6; ++i)
        annotations[0].p[i] = -1;
    for (int i = 0; i < 5; ++i)
        annotations[0].r[i] = -1;
    for (int i = 0; i < 4; ++i)
        annotations[0].t[i] = -1;
    annotations[0].p[0] = isoel_p_l;
    annotations[0].p[1] = peak_p1;
    annotations[0].p[2] = isoel_p_r;

    if (p2_indx != -1)
    {
        int middle_indx_from_l = peak_p1 + (peak_p1 - isoel_p_l);
        int middle_indx_from_r = p2_indx - (isoel_p_r - p2_indx);
        annotations[0].p[3] = (middle_indx_from_l + middle_indx_from_r) / 2;
        annotations[0].p[2] = annotations[0].p[3];
        annotations[0].p[4] = p2_indx;
        annotations[0].p[5] = isoel_p_r;
    }

    int pad_samples = c.timing.pq_interval_offset_ms * sampling_rate / 1000.0;
    annotations[0].r[0] = isoel_start - pad_samples;
    annotations[0].r[1] = peak_indexes[analysis_peak_indx];
    annotations[0].r[2] = isoel_r + pad_samples;

    annotations[0].t[0] = isoel_t_l;
    annotations[0].t[1] = peak_t;
    annotations[0].t[2] = isoel_t_r;

    fill_wave_measurements((double*)ecg_signal[analysis_ch_indx], nr_samples_per_ch, annotations[0], sampling_rate, c, result);

    delete[] lead_cpy;
}

static void fill_legacy_fields(ecg_analysis_result& result, double sampling_rate, const std::vector<unsigned int>& peak_indexes, const std::vector<pqrst_indxes>& annotations)
{
    /* ---- RR statistics ---- */
    if (peak_indexes.size() >= 2)
    {
        double total_ms = 0.0, max_rr = 0.0, min_rr = 1e9;
        std::vector<double> rr_intervals;
        for (size_t i = 1; i < peak_indexes.size(); ++i)
        {
            double rr_ms = (peak_indexes[i] - peak_indexes[i - 1]) / sampling_rate * 1000.0;
            rr_intervals.push_back(rr_ms);
            total_ms += rr_ms;
            if (rr_ms < min_rr) min_rr = rr_ms;
            if (rr_ms > max_rr) max_rr = rr_ms;
        }
        result.rr_interval_ms = total_ms / rr_intervals.size();
        result.rr_variation_ms = max_rr - min_rr;
        result.heart_rate_bpm = (result.rr_interval_ms > 0.0) ? 60000.0 / result.rr_interval_ms : 0.0;

        result.premature_beat_count = 0;
        for (double v : rr_intervals)
            if (result.rr_interval_ms > 0.0 && v < 0.80 * result.rr_interval_ms)
                result.premature_beat_count++;
        result.is_sinus_rhythm = (result.premature_beat_count || (result.rr_interval_ms > 0.0 && (max_rr - min_rr) / result.rr_interval_ms > 0.10)) ? 0 : 1;
    }

    result.pr_segment_ms = result.pr_interval_ms - result.p_wave_duration_ms;
    double rr_sec = result.rr_interval_ms / 1000.0;
    result.qtc_bazett_ms = (rr_sec > 0.0 && result.qt_interval_ms > 0.0) ? result.qt_interval_ms / std::sqrt(rr_sec) : 0.0;

    /* ---- T wave duration & ST segment from annotations ---- */
    if (!annotations.empty())
    {
        const pqrst_indxes& a = annotations[0];
        if (a.t[0] >= 0 && a.t[2] >= 0)
            result.t_wave_duration_ms = (a.t[2] - a.t[0]) / sampling_rate * 1000.0;
        if (a.t[0] >= 0 && a.r[2] >= 0)
            result.st_segment_ms = (a.t[0] - a.r[2]) / sampling_rate * 1000.0;
    }

}

ecg_analysis_result analyse_ecg_detect_peaks(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, std::vector<pqrst_indxes>& annotations, std::vector<unsigned int>* peak_indexes, std::string mode, int analysis_ch_indx, int analysis_peak_indx)
{
    ecg_analysis_config c;
    if (analysis_ch_indx < 0) analysis_ch_indx = 0;
    if (analysis_ch_indx >= (int)nr_channels) analysis_ch_indx = nr_channels - 1;

    std::vector<double> peak_signal(nr_samples_per_channel), filt_signal(nr_samples_per_channel), trshld_signal(nr_samples_per_channel);
    bool local_peak_indexes_used = false;
    if (!peak_indexes)
    {
        peak_indexes = new std::vector<unsigned int>;
        local_peak_indexes_used = true;
    }
    peak_detector_offline detector(sampling_rate);
    if (mode == "high_sensitivity")
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == "high_ppv")
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);

    //detector.detect_multichannel(data, nr_channels, nr_samples_per_channel, peak_signal.data(), filt_signal.data(), trshld_signal.data(), peak_indexes);
    detector.detect(data[analysis_ch_indx], nr_samples_per_channel, peak_signal.data(), filt_signal.data(), trshld_signal.data(), peak_indexes);
    ecg_analysis_result result;
    if (true && peak_indexes->size() == 0)
    {
        double min_val, max_val;
        int min_indx = 0, max_indx = 0;
        min_max_idx(data[analysis_ch_indx], nr_samples_per_channel, min_val, max_val, min_indx, max_indx);
        if (peak_indexes->size() == 0)
        {
            double baseline = median(data[analysis_ch_indx], nr_samples_per_channel);
            if ((max_val - baseline) > (baseline - min_val))
                peak_indexes->push_back(max_indx);
            else
                peak_indexes->push_back(min_indx);
        }
    }
    analyse_ecg((double**)data, (unsigned int)nr_channels, nr_samples_per_channel, sampling_rate, *peak_indexes, annotations, result, c, analysis_ch_indx, analysis_peak_indx);

    fill_legacy_fields(result, sampling_rate, *peak_indexes, annotations);

    if (local_peak_indexes_used)
        delete peak_indexes;
    return result;
}

void analyse_ecg_all_beats(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, std::vector<ecg_analysis_result>& results, int analysis_ch_indx)
{
    annotations.clear();
    results.clear();

    if (analysis_ch_indx < 0)
        analysis_ch_indx = 0;

    ecg_analysis_config c;
    annotations.reserve(peak_indexes.size());
    results.reserve(peak_indexes.size());

    for (size_t peak_index = 0; peak_index < peak_indexes.size(); ++peak_index)
    {
        std::vector<pqrst_indxes> beat_annotations;
        ecg_analysis_result beat_result = {};

        analyse_ecg((double**)data,
                    (unsigned int)nr_channels,
                    (unsigned int)nr_samples_per_channel,
                    sampling_rate,
                    peak_indexes,
                    beat_annotations,
                    beat_result,
                    c,
                    (unsigned int)analysis_ch_indx,
                    (int)peak_index);

        fill_legacy_fields(beat_result, sampling_rate, peak_indexes, beat_annotations);

        if (!beat_annotations.empty())
            annotations.push_back(beat_annotations[0]);
        else
            annotations.push_back({});

        results.push_back(beat_result);
    }
}

namespace
{

double missing_value()
{
    return std::numeric_limits<double>::quiet_NaN();
}

bool is_valid_mode(int32_t mode)
{
    return mode == RSPT_MODE_DEFAULT ||
           mode == RSPT_MODE_HIGH_SENSITIVITY ||
           mode == RSPT_MODE_HIGH_PPV;
}

void configure_detector(peak_detector_offline& detector, int32_t mode)
{
    if (mode == RSPT_MODE_HIGH_SENSITIVITY)
        detector.set_mode(peak_detector_offline::Mode::high_sensitivity);
    else if (mode == RSPT_MODE_HIGH_PPV)
        detector.set_mode(peak_detector_offline::Mode::high_ppv);
    else
        detector.set_mode(peak_detector_offline::Mode::def);
}

void copy_message(char* dest, const char* message)
{
    if (!dest)
        return;
    std::snprintf(dest, kStatusMessageSize, "%s", message ? message : "");
}

bool positive_finite(double value)
{
    return std::isfinite(value) && value > 0.0;
}

bool non_negative_finite(double value)
{
    return std::isfinite(value) && value >= 0.0;
}

void set_if_positive(double value, double& field, uint64_t& valid_fields, uint64_t bit)
{
    if (positive_finite(value))
    {
        field = value;
        valid_fields |= bit;
    }
}

void set_if_non_negative(double value, double& field, uint64_t& valid_fields, uint64_t bit)
{
    if (non_negative_finite(value))
    {
        field = value;
        valid_fields |= bit;
    }
}

void initialise_annotation(rspt_pqrst_annotation& annotation)
{
    annotation.p1_onset_sample = -1;
    annotation.p1_peak_sample = -1;
    annotation.p1_offset_sample = -1;
    annotation.p2_onset_sample = -1;
    annotation.p2_peak_sample = -1;
    annotation.p2_offset_sample = -1;
    annotation.qrs_onset_sample = -1;
    annotation.r_peak_sample = -1;
    annotation.qrs_offset_sample = -1;
    annotation.q_peak_sample = -1;
    annotation.s_peak_sample = -1;
    annotation.t_onset_sample = -1;
    annotation.t_peak_sample = -1;
    annotation.t_offset_sample = -1;
    annotation.j_point_sample = -1;
}

void initialise_beat_result(rspt_ecg_beat_result& result)
{
    std::memset(&result, 0, sizeof(result));
    result.status = RSPT_STATUS_OK;
    copy_message(result.status_message, "OK");
    initialise_annotation(result.annotation);

    result.rr_interval_ms = missing_value();
    result.heart_rate_bpm = missing_value();
    result.p_wave_duration_ms = missing_value();
    result.p1_wave_duration_ms = missing_value();
    result.p1_amplitude_input_units = missing_value();
    result.p2_wave_duration_ms = missing_value();
    result.p2_amplitude_input_units = missing_value();
    result.pr_interval_ms = missing_value();
    result.pr_segment_ms = missing_value();
    result.q_duration_ms = missing_value();
    result.q_amplitude_input_units = missing_value();
    result.r_duration_ms = missing_value();
    result.r_amplitude_input_units = missing_value();
    result.s_duration_ms = missing_value();
    result.s_amplitude_input_units = missing_value();
    result.qrs_duration_ms = missing_value();
    result.qt_interval_ms = missing_value();
    result.qtc_bazett_ms = missing_value();
    result.st_segment_ms = missing_value();
    result.t_wave_duration_ms = missing_value();
    result.j_point_amplitude_input_units = missing_value();
    result.st20_amplitude_input_units = missing_value();
    result.st40_amplitude_input_units = missing_value();
    result.st60_amplitude_input_units = missing_value();
    result.st80_amplitude_input_units = missing_value();
    result.t_amplitude_input_units = missing_value();
}

void initialise_metric_statistics(rspt_metric_statistics& stats)
{
    stats.count = 0;
    stats.mean = missing_value();
    stats.standard_deviation = missing_value();
}

void initialise_summary_result(rspt_ecg_summary_result& result)
{
    std::memset(&result, 0, sizeof(result));
    result.status = RSPT_STATUS_OK;
    copy_message(result.status_message, "OK");

    initialise_metric_statistics(result.rr_interval_ms);
    result.rr_variation_ms = missing_value();
    initialise_metric_statistics(result.heart_rate_bpm);
    initialise_metric_statistics(result.p_wave_duration_ms);
    initialise_metric_statistics(result.pr_interval_ms);
    initialise_metric_statistics(result.qrs_duration_ms);
    initialise_metric_statistics(result.qt_interval_ms);
    initialise_metric_statistics(result.qtc_bazett_ms);
    initialise_metric_statistics(result.st_segment_ms);
    initialise_metric_statistics(result.t_wave_duration_ms);
    result.is_sinus_rhythm = -1;
    result.premature_beat_count = -1;
}

int32_t detect_peak_indexes(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    std::vector<unsigned int>& peak_indexes)
{
    peak_indexes.clear();

    if (!signal || sample_count == 0 || !std::isfinite(sampling_rate) || sampling_rate <= 0.0)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (sample_count > static_cast<size_t>(std::numeric_limits<unsigned int>::max()))
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!is_valid_mode(mode))
        return RSPT_STATUS_INVALID_ARGUMENT;

    std::vector<double> peak_signal(sample_count);
    std::vector<double> filtered_signal(sample_count);
    std::vector<double> threshold_signal(sample_count);

    peak_detector_offline detector(sampling_rate);
    configure_detector(detector, mode);
    detector.detect(
        signal,
        static_cast<unsigned int>(sample_count),
        peak_signal.data(),
        filtered_signal.data(),
        threshold_signal.data(),
        &peak_indexes);

    std::sort(peak_indexes.begin(), peak_indexes.end());
    peak_indexes.erase(std::unique(peak_indexes.begin(), peak_indexes.end()), peak_indexes.end());

    return peak_indexes.empty() ? RSPT_STATUS_NO_R_PEAKS : RSPT_STATUS_OK;
}

struct AnalysisData
{
    int32_t analysis_channel_index = 0;
    std::vector<unsigned int> peak_indexes;
    std::vector<pqrst_indxes> annotations;
    std::vector<ecg_analysis_result> results;
};

int32_t run_core_analysis(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    AnalysisData& analysis)
{
    analysis = AnalysisData{};

    if (!channels || channel_count == 0 || samples_per_channel == 0 || !std::isfinite(sampling_rate) || sampling_rate <= 0.0)
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (channel_count > static_cast<size_t>(std::numeric_limits<unsigned int>::max()) ||
        samples_per_channel > static_cast<size_t>(std::numeric_limits<unsigned int>::max()))
        return RSPT_STATUS_INVALID_ARGUMENT;
    if (!is_valid_mode(mode))
        return RSPT_STATUS_INVALID_ARGUMENT;

    for (size_t channel = 0; channel < channel_count; ++channel)
    {
        if (!channels[channel])
            return RSPT_STATUS_NO_CHANNEL_DATA;
    }

    if (analysis_channel_index < 0)
        analysis_channel_index = 0;
    if (static_cast<size_t>(analysis_channel_index) >= channel_count)
        return RSPT_STATUS_INVALID_CHANNEL_INDEX;
    analysis.analysis_channel_index = analysis_channel_index;

    if (r_peak_indexes && r_peak_count > 0)
    {
        analysis.peak_indexes.reserve(r_peak_count);
        for (size_t i = 0; i < r_peak_count; ++i)
        {
            if (r_peak_indexes[i] >= samples_per_channel)
                return RSPT_STATUS_INVALID_PEAK_INDEX;
            analysis.peak_indexes.push_back(static_cast<unsigned int>(r_peak_indexes[i]));
        }
        std::sort(analysis.peak_indexes.begin(), analysis.peak_indexes.end());
        analysis.peak_indexes.erase(std::unique(analysis.peak_indexes.begin(), analysis.peak_indexes.end()), analysis.peak_indexes.end());
    }
    else
    {
        int32_t status = detect_peak_indexes(
            channels[analysis.analysis_channel_index],
            samples_per_channel,
            sampling_rate,
            mode,
            analysis.peak_indexes);
        if (status != RSPT_STATUS_OK)
            return status;
    }

    if (analysis.peak_indexes.empty())
        return RSPT_STATUS_NO_R_PEAKS;

    std::vector<const double*> channel_ptrs(channel_count);
    for (size_t channel = 0; channel < channel_count; ++channel)
        channel_ptrs[channel] = channels[channel];

    analyse_ecg_all_beats(
        channel_ptrs.data(),
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis.peak_indexes,
        analysis.annotations,
        analysis.results,
        analysis.analysis_channel_index);

    return RSPT_STATUS_OK;
}

void copy_annotation(const pqrst_indxes& source, rspt_pqrst_annotation& destination, uint64_t& valid_fields)
{
    initialise_annotation(destination);

    destination.p1_onset_sample = source.p[0];
    destination.p1_peak_sample = source.p[1];
    destination.p1_offset_sample = source.p[2];

    if (source.p[3] >= 0 && source.p[4] >= 0 && source.p[5] >= 0)
    {
        destination.p2_onset_sample = source.p[3];
        destination.p2_peak_sample = source.p[4];
        destination.p2_offset_sample = source.p[5];
    }

    destination.qrs_onset_sample = source.r[0];
    destination.r_peak_sample = source.r[1];
    destination.qrs_offset_sample = source.r[2];
    destination.q_peak_sample = source.r[3];
    destination.s_peak_sample = source.r[4];

    destination.t_onset_sample = source.t[0];
    destination.t_peak_sample = source.t[1];
    destination.t_offset_sample = source.t[2];
    destination.j_point_sample = source.t[3];

    if (source.p[0] >= 0 && source.p[1] >= 0 && source.p[2] >= 0 &&
        source.r[0] >= 0 && source.r[1] >= 0 && source.r[2] >= 0 &&
        source.t[0] >= 0 && source.t[1] >= 0 && source.t[2] >= 0)
    {
        valid_fields |= RSPT_VALID_PQRS_T_ANNOTATION;
    }
}

bool rr_interval_ms_for_index(
    const std::vector<unsigned int>& peak_indexes,
    size_t right_peak_index,
    double sampling_rate,
    double& rr_ms)
{
    if (right_peak_index == 0 || right_peak_index >= peak_indexes.size())
        return false;
    if (peak_indexes[right_peak_index] <= peak_indexes[right_peak_index - 1])
        return false;

    rr_ms = (peak_indexes[right_peak_index] - peak_indexes[right_peak_index - 1]) / sampling_rate * 1000.0;
    return positive_finite(rr_ms);
}

void fill_beat_result(
    const ecg_analysis_result& core_result,
    const pqrst_indxes* annotation,
    const std::vector<unsigned int>& peak_indexes,
    size_t beat_index,
    double sampling_rate,
    uint32_t analysis_channel_index,
    rspt_ecg_beat_result& output)
{
    initialise_beat_result(output);
    output.status = core_result.analysis_status;
    output.analysis_channel_index = analysis_channel_index;
    output.beat_index = static_cast<uint32_t>(beat_index);

    if (beat_index < peak_indexes.size())
    {
        output.r_peak_sample = peak_indexes[beat_index];
        output.valid_fields |= RSPT_VALID_R_PEAK_SAMPLE;
    }

    copy_message(output.status_message, core_result.status_message[0] ? core_result.status_message : rspt::status_message(output.status));

    if (core_result.analysis_status != RSPT_STATUS_OK)
        return;

    if (annotation)
        copy_annotation(*annotation, output.annotation, output.valid_fields);

    set_if_positive(core_result.p_wave_duration_ms, output.p_wave_duration_ms, output.valid_fields, RSPT_VALID_P_WAVE_DURATION_MS);
    output.p1_wave_duration_ms = core_result.p1_wave_duration_ms;
    output.p1_amplitude_input_units = core_result.p1_amplitude_input_units;
    output.p2_wave_duration_ms = core_result.p2_wave_duration_ms;
    output.p2_amplitude_input_units = core_result.p2_amplitude_input_units;

    set_if_positive(core_result.pr_interval_ms, output.pr_interval_ms, output.valid_fields, RSPT_VALID_PR_INTERVAL_MS);
    set_if_non_negative(core_result.pr_segment_ms, output.pr_segment_ms, output.valid_fields, RSPT_VALID_PR_SEGMENT_MS);

    output.q_duration_ms = core_result.q_duration_ms;
    output.q_amplitude_input_units = core_result.q_amplitude_input_units;
    output.r_duration_ms = core_result.r_duration_ms;
    output.r_amplitude_input_units = core_result.r_amplitude_input_units;
    output.s_duration_ms = core_result.s_duration_ms;
    output.s_amplitude_input_units = core_result.s_amplitude_input_units;

    set_if_positive(core_result.qrs_duration_ms, output.qrs_duration_ms, output.valid_fields, RSPT_VALID_QRS_DURATION_MS);
    set_if_positive(core_result.qt_interval_ms, output.qt_interval_ms, output.valid_fields, RSPT_VALID_QT_INTERVAL_MS);
    set_if_non_negative(core_result.st_segment_ms, output.st_segment_ms, output.valid_fields, RSPT_VALID_ST_SEGMENT_MS);
    set_if_positive(core_result.t_wave_duration_ms, output.t_wave_duration_ms, output.valid_fields, RSPT_VALID_T_WAVE_DURATION_MS);

    output.j_point_amplitude_input_units = core_result.j_point_amplitude_input_units;
    output.st20_amplitude_input_units = core_result.st20_amplitude_input_units;
    output.st40_amplitude_input_units = core_result.st40_amplitude_input_units;
    output.st60_amplitude_input_units = core_result.st60_amplitude_input_units;
    output.st80_amplitude_input_units = core_result.st80_amplitude_input_units;
    output.t_amplitude_input_units = core_result.t_amplitude_input_units;

    double previous_rr_ms = 0.0;
    if (rr_interval_ms_for_index(peak_indexes, beat_index, sampling_rate, previous_rr_ms))
    {
        output.rr_interval_ms = previous_rr_ms;
        output.heart_rate_bpm = 60000.0 / previous_rr_ms;
        output.valid_fields |= RSPT_VALID_RR_INTERVAL_MS | RSPT_VALID_HEART_RATE_BPM;
        if (positive_finite(output.qt_interval_ms))
        {
            output.qtc_bazett_ms = output.qt_interval_ms / std::sqrt(previous_rr_ms / 1000.0);
            output.valid_fields |= RSPT_VALID_QTC_BAZETT_MS;
        }
    }
}

void add_to_stats(double value, std::vector<double>& values)
{
    if (positive_finite(value))
        values.push_back(value);
}

void set_stats(const std::vector<double>& values, rspt_metric_statistics& destination, uint64_t& valid_fields, uint64_t bit)
{
    if (values.empty())
        return;

    destination.count = values.size();
    destination.mean = std::accumulate(values.begin(), values.end(), 0.0) / static_cast<double>(values.size());

    double sum_squares = 0.0;
    for (double value : values)
    {
        double diff = value - destination.mean;
        sum_squares += diff * diff;
    }
    destination.standard_deviation = std::sqrt(sum_squares / static_cast<double>(values.size()));
    valid_fields |= bit;
}

}

namespace rspt
{

const char* status_message(int32_t status)
{
    switch (status)
    {
        case RSPT_STATUS_OK:
            return "OK";
        case RSPT_STATUS_NO_R_PEAKS:
            return "No R peaks found";
        case RSPT_STATUS_NO_CHANNEL_DATA:
            return "No channel data given";
        case RSPT_STATUS_INVALID_ARGUMENT:
            return "Invalid argument";
        case RSPT_STATUS_INVALID_CHANNEL_INDEX:
            return "Invalid channel index";
        case RSPT_STATUS_INVALID_PEAK_INDEX:
            return "Invalid peak index";
        case RSPT_STATUS_OUTPUT_BUFFER_TOO_SMALL:
            return "Output buffer too small";
        default:
            return "Unknown status";
    }
}

uint32_t api_version()
{
    return 4;
}

int32_t detect_peaks_double(
    const double* signal,
    size_t sample_count,
    double sampling_rate,
    int32_t mode,
    std::vector<uint32_t>& out_r_peak_indexes)
{
    std::vector<unsigned int> peak_indexes;
    int32_t status = detect_peak_indexes(signal, sample_count, sampling_rate, mode, peak_indexes);

    out_r_peak_indexes.clear();
    out_r_peak_indexes.reserve(peak_indexes.size());
    for (unsigned int peak_index : peak_indexes)
        out_r_peak_indexes.push_back(static_cast<uint32_t>(peak_index));

    return status;
}

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
    std::vector<uint32_t>* out_detected_r_peak_indexes)
{
    AnalysisData analysis;
    int32_t status = run_core_analysis(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis_channel_index,
        r_peak_indexes,
        r_peak_count,
        mode,
        analysis);

    out_beats.clear();

    if (out_detected_r_peak_indexes)
    {
        out_detected_r_peak_indexes->clear();
        out_detected_r_peak_indexes->reserve(analysis.peak_indexes.size());
        for (unsigned int peak_index : analysis.peak_indexes)
            out_detected_r_peak_indexes->push_back(static_cast<uint32_t>(peak_index));
    }

    if (status != RSPT_STATUS_OK)
        return status;

    out_beats.resize(analysis.results.size());
    for (size_t i = 0; i < analysis.results.size(); ++i)
    {
        const pqrst_indxes* annotation = (i < analysis.annotations.size()) ? &analysis.annotations[i] : nullptr;
        fill_beat_result(
            analysis.results[i],
            annotation,
            analysis.peak_indexes,
            i,
            sampling_rate,
            static_cast<uint32_t>(analysis.analysis_channel_index),
            out_beats[i]);
    }

    return RSPT_STATUS_OK;
}

int32_t analyze_ecg_summary_double(
    const double* const* channels,
    size_t channel_count,
    size_t samples_per_channel,
    double sampling_rate,
    int32_t analysis_channel_index,
    const uint32_t* r_peak_indexes,
    size_t r_peak_count,
    int32_t mode,
    rspt_ecg_summary_result& out_summary)
{
    initialise_summary_result(out_summary);

    std::vector<rspt_ecg_beat_result> beats;
    std::vector<uint32_t> detected_r_peak_indexes;
    int32_t status = analyze_ecg_beats_double(
        channels,
        channel_count,
        samples_per_channel,
        sampling_rate,
        analysis_channel_index,
        r_peak_indexes,
        r_peak_count,
        mode,
        beats,
        &detected_r_peak_indexes);

    out_summary.status = status;
    out_summary.analysis_channel_index = analysis_channel_index < 0 ? 0 : static_cast<uint32_t>(analysis_channel_index);
    out_summary.r_peak_count = detected_r_peak_indexes.size();
    copy_message(out_summary.status_message, status_message(status));

    if (status != RSPT_STATUS_OK)
        return status;

    out_summary.analysed_beat_count = 0;
    int32_t first_failure_status = RSPT_STATUS_OK;
    std::vector<double> rr_values, hr_values, p_values, pr_values, qrs_values, qt_values, qtc_values, st_values, t_values;
    double rr_min = std::numeric_limits<double>::max();
    double rr_max = 0.0;

    for (const rspt_ecg_beat_result& beat : beats)
    {
        if (beat.status != RSPT_STATUS_OK)
        {
            if (first_failure_status == RSPT_STATUS_OK)
                first_failure_status = beat.status;
            continue;
        }

        ++out_summary.analysed_beat_count;

        if (positive_finite(beat.rr_interval_ms))
        {
            rr_min = std::min(rr_min, beat.rr_interval_ms);
            rr_max = std::max(rr_max, beat.rr_interval_ms);
        }

        add_to_stats(beat.rr_interval_ms, rr_values);
        add_to_stats(beat.heart_rate_bpm, hr_values);
        add_to_stats(beat.p_wave_duration_ms, p_values);
        add_to_stats(beat.pr_interval_ms, pr_values);
        add_to_stats(beat.qrs_duration_ms, qrs_values);
        add_to_stats(beat.qt_interval_ms, qt_values);
        add_to_stats(beat.qtc_bazett_ms, qtc_values);
        add_to_stats(beat.st_segment_ms, st_values);
        add_to_stats(beat.t_wave_duration_ms, t_values);
    }

    if (out_summary.analysed_beat_count == 0)
    {
        out_summary.status = first_failure_status == RSPT_STATUS_OK ? RSPT_STATUS_INVALID_PEAK_INDEX : first_failure_status;
        copy_message(out_summary.status_message, status_message(out_summary.status));
        return out_summary.status;
    }

    set_stats(rr_values, out_summary.rr_interval_ms, out_summary.valid_fields, RSPT_VALID_RR_INTERVAL_MS);
    set_stats(hr_values, out_summary.heart_rate_bpm, out_summary.valid_fields, RSPT_VALID_HEART_RATE_BPM);
    set_stats(p_values, out_summary.p_wave_duration_ms, out_summary.valid_fields, RSPT_VALID_P_WAVE_DURATION_MS);
    set_stats(pr_values, out_summary.pr_interval_ms, out_summary.valid_fields, RSPT_VALID_PR_INTERVAL_MS);
    set_stats(qrs_values, out_summary.qrs_duration_ms, out_summary.valid_fields, RSPT_VALID_QRS_DURATION_MS);
    set_stats(qt_values, out_summary.qt_interval_ms, out_summary.valid_fields, RSPT_VALID_QT_INTERVAL_MS);
    set_stats(qtc_values, out_summary.qtc_bazett_ms, out_summary.valid_fields, RSPT_VALID_QTC_BAZETT_MS);
    set_stats(st_values, out_summary.st_segment_ms, out_summary.valid_fields, RSPT_VALID_ST_SEGMENT_MS);
    set_stats(t_values, out_summary.t_wave_duration_ms, out_summary.valid_fields, RSPT_VALID_T_WAVE_DURATION_MS);

    if (!rr_values.empty())
    {
        out_summary.rr_variation_ms = rr_max - rr_min;

        int32_t premature_count = 0;
        for (double rr_ms : rr_values)
        {
            if (rr_ms < 0.80 * out_summary.rr_interval_ms.mean)
                ++premature_count;
        }
        out_summary.premature_beat_count = premature_count;
        out_summary.is_sinus_rhythm =
            (premature_count == 0 &&
             out_summary.rr_interval_ms.mean > 0.0 &&
             out_summary.rr_variation_ms / out_summary.rr_interval_ms.mean <= 0.10)
            ? 1
            : 0;
    }

    out_summary.status = RSPT_STATUS_OK;
    copy_message(out_summary.status_message, "OK");
    return RSPT_STATUS_OK;
}

}
