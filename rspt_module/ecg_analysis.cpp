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
#include <limits>

using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "ecg_analysis.h"
#include "filter.h"
#include "iir_filter_opt.h"
#include "peak_detector.h"

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
    return (int)idx_extreme;
}

//static int find_min_from_to(const double* signal, unsigned int from, unsigned int to, unsigned int length)
//{
//    return find_extreme_from_to(signal, from, to, length, std::less<double>());
//}

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

void write_binmx_to_file(const char* filename, const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate)
{
    FILE* f = fopen(filename, "wb");
    if (!f)
    {
        perror("fopen");
        return;
    }

    // --- HEADER ---
    // 0–7: sampling rate
    fwrite(&sampling_rate, sizeof(double), 1, f);

    // 8–11: nr_channels (int32)
    int32_t ch = (int32_t)nr_channels;
    fwrite(&ch, sizeof(int32_t), 1, f);

    // 12–23: padding (zeros)
    uint8_t padding[12] = {0};
    fwrite(padding, sizeof(padding), 1, f);

    // --- DATA ---
    // interleaved: ch0_s0, ch1_s0, ..., chN_s0, ch0_s1, ...
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

void write_binmx_to_file_1ch(const char* filename, const double* data, size_t nr_samples, double sampling_rate)
{
    FILE* f = fopen(filename, "wb");
    if (!f)
    {
        perror("fopen");
        return;
    }

    // --- HEADER ---
    fwrite(&sampling_rate, sizeof(double), 1, f);

    int32_t nr_channels = 1;  // csak 1 csatorna
    fwrite(&nr_channels, sizeof(int32_t), 1, f);

    uint8_t padding[12] = {0};
    fwrite(padding, sizeof(padding), 1, f);

    // --- DATA ---
    // csak lineárisan kiírjuk a mintákat
    fwrite(data, sizeof(double), nr_samples, f);

    fclose(f);
}

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

    //cout << "---------------------------------------------------------------------" << endl;
    for (unsigned int i = start_idx; i <= end_idx; ++i)
    {
        double baseline = m * (double)(i) + b;
        double corrected = signal[i] - baseline;
        //cout << corrected << " \t " << signal[i] << endl;
        if (corrected < min_val)
        {
            min_val = corrected;
            min_idx = i;
        }
    }
    //write_binmx_to_file("/media/sf_SharedFolder/test01.bin", &signal, 1, end_idx - start_idx, 250);
//    static int written = false;
//    if (!written)
//    {
//        const double* tmp = &signal[start_idx];
//        write_binmx_to_file("/media/sf_SharedFolder/test01.bin", &tmp, 1, end_idx - start_idx, 250);
//        //write_binmx_to_file("/media/sf_SharedFolder/test01.bin", &signal, 1, end_idx - start_idx, 250);
//        //write_binmx_to_file_1ch("/media/sf_SharedFolder/test01.bin", &signal[start_idx], end_idx - start_idx, 250);
//    }
//    written = true;
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
    if (start >= end)
    {
        peak_indx = (int)start;
        sign = 1;
        return 0.0;
    }

    const size_t n = end - start;
    if (n < 3)
    {
        peak_indx = (int)start;
        sign = 1;
        return 0.0;
    }

    // 1) Másolás + lineáris detrend (vonal kiegyenesítés)
    std::vector<double> y(n);
    for (size_t i = 0; i < n; ++i) y[i] = sig[start + i];

    const double y0 = y.front();
    const double yN = y.back();
    const double denom = (n > 1) ? double(n - 1) : 1.0;
    for (size_t i = 0; i < n; ++i)
    {
        double baseline = y0 + (yN - y0) * (double(i) / denom);
        y[i] -= baseline;
    }

    // 2) Simítás: kétszeres középre igazított mozgóátlag (gyors komponensek elnyomása)
    auto smooth_centered = [&](const std::vector<double>& in, int win)->std::vector<double>
    {
        std::vector<double> out(n);
        if (win < 3) return in;
        if ((win & 1) == 0) ++win; // legyen páratlan
        const int h = win / 2;

        // prefix összegek a gyors közép-átlaghoz
        std::vector<double> ps(n + 1, 0.0);
        for (size_t i = 0; i < n; ++i) ps[i + 1] = ps[i] + in[i];

        for (size_t i = 0; i < n; ++i)
        {
            int a = (int)i - h;
            if (a < 0) a = 0;
            int b = (int)i + h;
            if (b >= (int)n) b = (int)n - 1;
            double sum = ps[b + 1] - ps[a];
            out[i] = sum / double(b - a + 1);
        }
        return out;
    };

    int w1 = std::max<int>(3, int(n / 20)); // ~5% ablak
    if ((w1 & 1) == 0) ++w1;
    std::vector<double> y1 = smooth_centered(y, w1);

    int w2 = std::max<int>(3, w1 / 2);
    if ((w2 & 1) == 0) ++w2;
    std::vector<double> y2 = smooth_centered(y1, w2);

    // 3) Keresés csak a középső sávban (edge-ek kerülése)
    size_t margin = std::max<size_t>(1, n / 10);           // 10% szélek elhagyása
    size_t lo = margin;
    size_t hi = (n > margin) ? (n - 1 - margin) : (n - 1);
    if (hi <= lo)
    {
        lo = 0;
        hi = n - 1;
    }

    double best_score = -std::numeric_limits<double>::infinity();
    size_t best_i = lo;
    const double center = 0.5 * (lo + hi);

    for (size_t i = lo; i <= hi; ++i)
    {
        // enyhe középre súlyozás, hogy a középsáv preferált legyen
        double centrality = 1.0 - 0.25 * std::fabs(double(i) - center) / (center + 1e-9);
        double score = std::fabs(y2[i]) * centrality;
        if (score > best_score)
        {
            best_score = score;
            best_i = i;
        }
    }

    const double val = y2[best_i];
    peak_indx = int(start + best_i);
    sign = (val >= 0.0) ? 1 : -1;
    return std::fabs(val);
}

int xcor_align(double** leads, size_t nr_samples, size_t ref_ch, size_t ch, size_t t_start, size_t t_end, size_t max_lag_samples)
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

void accumulate_aligned(double* res, const double* src, size_t /**nr_src_samples*/, size_t t_start, size_t t_end, int /**dt*/, double w)
{
//    if (t_end > nr_src_samples) t_end = nr_src_samples;
    size_t len = t_end - t_start;
//    if (len == 0) return;

    for (size_t i = 0; i < len; ++i)
    {
        int si = (int)t_start + (int)i + 0;
//        if (si >= 0 && si < (int)nr_src_samples)
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

int create_ideal_signal(double* res, double** leads, size_t nr_channels, size_t nr_samples, double sampling_rate, size_t t_start, size_t t_end, size_t iso_start, size_t iso_end)
{
    if (t_end > nr_samples) t_end = nr_samples;
    size_t len = t_end - t_start;
    if (len == 0) return -1;

    // 1) SNR becslés és rangsor
    std::vector<LeadScore> scores;
    scores.reserve(nr_channels);
    for (size_t ch = 0; ch < nr_channels; ++ch)
    {
///        cout << "---------------------------------------------------------------------" << endl;
        ///     for (size_t i = t_start; i < t_end; ++i)
        ///      cout << leads[ch][i] << endl;
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
        int max_lag_samples = (int)(0.015 * sampling_rate);
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
        w = 1.0 / (double)K;
        accumulate_aligned(res, leads[scores[k].ch], nr_samples, t_start, t_end, scores[k].dt, scores[k].sign * w);
        //cout << "scores[k].sign: " << scores[k].sign << " leads[scores[k].ch] " << leads[scores[k].ch][0] << " res[0] " << res[0] << " W: " << w << endl;
    }

    return 0;
}

static inline size_t clamp_index(long idx, size_t len)
{
    if (idx < 0) return 0;
    if ((size_t)idx >= len) return len - 1;
    return (size_t)idx;
}

void detect_qrs(double* ecg, size_t len, double fs, size_t orig_r_indx, size_t& q_indx, size_t& r_indx, size_t& s_indx)
{
    // default outputs (safe)
    q_indx = r_indx = s_indx = orig_r_indx;

    if (!ecg || len == 0) return;

    // convert orig index safely
    size_t orig = std::min(orig_r_indx, len - 1);

    // milliseconds -> samples helpers
    auto ms2samps = [&](double ms)->int
    {
        return (int)std::max(1.0, std::round(ms * fs / 1000.0));
    };

    // window sizes (tuned for clinical ECG; can be adjusted)
    int baseline_ms = 300;    // baseline median over ±300 ms
    int local_win_ms = 150;   // search ±150 ms for peaks
    int q_search_ms = 80;     // search up to 80 ms for Q / S
    int min_peak_sep_ms = 20; // minimal separation accepted between biphasic peaks

    int baseline_w = ms2samps(baseline_ms);
    int win = ms2samps(local_win_ms);
    int qwin = ms2samps(q_search_ms);
    int min_peak_sep = ms2samps(min_peak_sep_ms);

    // baseline median computation in [orig - baseline_w, orig + baseline_w]
    int bL = (int)orig - baseline_w;
    int bR = (int)orig + baseline_w;
    bL = std::max(bL, 0);
    bR = std::min(bR, (int)len - 1);

    std::vector<double> tmp;
    tmp.reserve(bR - bL + 1);
    for (int i = bL; i <= bR; ++i) tmp.push_back(ecg[i]);
    // compute median
    std::nth_element(tmp.begin(), tmp.begin() + tmp.size()/2, tmp.end());
    double baseline = tmp[tmp.size()/2];
//    if (tmp.size() >= 2)
//    {
//        // for even size, this is an approximation but fine
//    }

    // create a baseline-corrected view via accessor lambda (avoid copying whole signal)
    auto val = [&](int idx)->double
    {
        idx = std::max(0, std::min((int)len - 1, idx));
        return ecg[idx] - baseline;
    };

    // search interval for peaks
    int L = std::max(0, (int)orig - win);
    int R = std::min((int)len - 1, (int)orig + win);

    // find local extrema in [L,R]
    struct Peak
    {
        int idx;
        double v;
        bool is_max;
    };
    std::vector<Peak> peaks;
    // detect local maxima/minima by simple slope sign changes
    for (int i = std::max(L+1,1); i <= std::min(R-1,(int)len-2); ++i)
    {
        double p = val(i-1);
        double c = val(i);
        double n = val(i+1);
        if (c >= p && c >= n)
        {
            peaks.push_back({i, c, true});
        }
        else if (c <= p && c <= n)
        {
            peaks.push_back({i, c, false});
        }
    }

    // If no peaks found inside window (rare), fall back to taking orig as R and search Q/S there.
    if (peaks.empty())
    {
        r_indx = orig;
        // find Q: local minimum in [orig - qwin, orig]
        int qL = std::max(0, (int)orig - qwin);
        int qR = (int)orig;
        int qmin = qL;
        for (int i = qL; i <= qR; ++i) if (val(i) < val(qmin)) qmin = i;
        q_indx = clamp_index(qmin, len);
        // find S: local minimum in [orig, orig + qwin]
        int sL = (int)orig;
        int sR = std::min((int)len - 1, (int)orig + qwin);
        int smin = sR;
        for (int i = sL; i <= sR; ++i) if (val(i) < val(smin)) smin = i;
        s_indx = clamp_index(smin, len);
        // ensure ordering
        if (!(q_indx < r_indx && r_indx < s_indx))
        {
            if (q_indx >= r_indx) q_indx = (r_indx > 0) ? r_indx - 1 : r_indx;
            if (s_indx <= r_indx) s_indx = std::min((size_t)r_indx + 1, len - 1);
        }
        return;
    }

    // find the peak nearest to orig (most relevant) and also the maximum abs peak in window
    int nearest_peak_idx = 0;
    int maxabs_peak_idx = 0;
    double best_dist = std::numeric_limits<double>::infinity();
    double maxabs = 0.0;
    for (size_t i = 0; i < peaks.size(); ++i)
    {
        double d = std::abs(peaks[i].idx - (int)orig);
        if (d < best_dist)
        {
            best_dist = d;
            nearest_peak_idx = (int)i;
        }
        double a = std::abs(peaks[i].v);
        if (a > maxabs)
        {
            maxabs = a;
            maxabs_peak_idx = (int)i;
        }
    }

    // choose candidate peaks for biphasic detection: look for two opposite-sign peaks within [L,R]
    // pick the two largest-abs peaks (if they have opposite signs and are reasonably close) => biphasic
    // build vector of peaks sorted by abs desc
    std::vector<int> idxs(peaks.size());
    for (size_t i = 0; i < peaks.size(); ++i) idxs[i] = (int)i;
    std::sort(idxs.begin(), idxs.end(), [&](int a, int b)
    {
        return std::abs(peaks[a].v) > std::abs(peaks[b].v);
    });

    bool is_biphasic = false;
    int left_peak = -1, right_peak = -1;
    if (idxs.size() >= 2)
    {
        // examine top few peaks (up to 4) to find a pair with opposite signs
        size_t topN = std::min((size_t)4, idxs.size());
        for (size_t i = 0; i < topN && !is_biphasic; ++i)
        {
            for (size_t j = i + 1; j < topN && !is_biphasic; ++j)
            {
                int ia = idxs[i], ib = idxs[j];
                if (peaks[ia].v * peaks[ib].v < 0.0)   // opposite signs
                {
                    int ia_idx = peaks[ia].idx;
                    int ib_idx = peaks[ib].idx;
                    int sep = std::abs(ia_idx - ib_idx);
                    // require that separation is at least some minimum (not identical samples)
                    if (sep >= min_peak_sep && sep <=  (int)(2*win))
                    {
                        // require both peaks are reasonably large (relative to maxabs in this window)
                        double a1 = std::abs(peaks[ia].v);
                        double a2 = std::abs(peaks[ib].v);
                        // thresholds: both >= 0.35*maxabs and max ratio not too extreme
                        double thr = 0.35;
                        double ratio = (std::max(a1,a2) / std::max(1e-12, std::min(a1,a2)));
                        if (a1 >= thr*maxabs && a2 >= thr*maxabs && ratio < 6.0)
                        {
                            // order them left->right
                            if (ia_idx < ib_idx)
                            {
                                left_peak = ia_idx;
                                right_peak = ib_idx;
                            }
                            else
                            {
                                left_peak = ib_idx;
                                right_peak = ia_idx;
                            }
                            is_biphasic = true;
                        }
                    }
                }
            }
        }
    }

    if (is_biphasic)
    {
        // center R between left_peak and right_peak (midpoint)
        int center = (left_peak + right_peak) / 2;
        r_indx = clamp_index(center, len);

        // find Q: search left from left_peak for baseline crossing (val ~ 0) or local minimum near start
        int q_search_L = std::max(0, left_peak - qwin*2); // allow a larger search just before biphasic start
        int q_search_R = left_peak;
        int q_found = -1;
        // prefer zero-crossing: find last index k in [q_search_L+1, q_search_R] s.t. val(k-1)*val(k) <= 0
        for (int k = q_search_R; k > q_search_L; --k)
        {
            if (val(k) == 0.0 || val(k-1) * val(k) <= 0.0)
            {
                q_found = k;
                break;
            }
        }
        if (q_found < 0)
        {
            // fallback: pick the local extremum (min abs value) near q_search_L..q_search_R
            int best = q_search_L;
            double bestAbs = std::abs(val(best));
            for (int k = q_search_L+1; k <= q_search_R; ++k)
            {
                double a = std::abs(val(k));
                if (a < bestAbs)
                {
                    bestAbs = a;
                    best = k;
                }
            }
            q_found = best;
        }
        q_indx = clamp_index(q_found, len);

        // find S: search right from right_peak for baseline crossing or low abs value
        int s_search_L = right_peak;
        int s_search_R = std::min((int)len - 1, right_peak + qwin*2);
        int s_found = -1;
        for (int k = s_search_L; k < s_search_R; ++k)
        {
            if (val(k) == 0.0 || val(k) * val(k+1) <= 0.0)
            {
                s_found = k;
                break;
            }
        }
        if (s_found < 0)
        {
            int best = s_search_L;
            double bestAbs = std::abs(val(best));
            for (int k = s_search_L+1; k <= s_search_R; ++k)
            {
                double a = std::abs(val(k));
                if (a < bestAbs)
                {
                    bestAbs = a;
                    best = k;
                }
            }
            s_found = best;
        }
        s_indx = clamp_index(s_found, len);

        // ensure ordering Q < R < S, if not fix by nudging within bounds
        if (!(q_indx < r_indx && r_indx < s_indx))
        {
            if (q_indx >= r_indx) q_indx = (r_indx > 0) ? r_indx - 1 : r_indx;
            if (s_indx <= r_indx) s_indx = std::min((size_t)r_indx + 1, len - 1);
        }

        return;
    }

    // Not biphasic: pick dominant peak near orig (prefer nearest or the maxabs)
    int chosen_peak_idx = peaks[nearest_peak_idx].idx;
    // if the nearest peak has small amplitude relative to maxabs, choose the maxabs peak (stability)
    if (std::abs(peaks[nearest_peak_idx].v) < 0.25 * maxabs)
    {
        chosen_peak_idx = peaks[maxabs_peak_idx].idx;
    }
    r_indx = clamp_index(chosen_peak_idx, len);

    // Determine if R is positive-going or negative-going
    double rval = val((int)r_indx);
    bool r_positive = (rval >= 0.0);

    // Q search: search left for a local minimum (for positive R) or local maximum (for negative R),
    // or baseline crossing
    int qL = std::max(0, (int)r_indx - qwin);
    int qR = (int)r_indx;
    int qbest = qR;
    if (r_positive)
    {
        // we want preceding minimum (Q)
        double bestv = val(qL);
        qbest = qL;
        for (int k = qL; k <= qR; ++k)
        {
            double v = val(k);
            if (v < bestv)
            {
                bestv = v;
                qbest = k;
            }
        }
        // if there's a zero crossing nearer to r_indx, pick that (iso-electrical start)
        for (int k = qR; k > qL; --k)
        {
            if (val(k) == 0.0 || val(k-1)*val(k) <= 0.0)
            {
                qbest = k;
                break;
            }
        }
    }
    else
    {
        // negative R -> preceding maximum (Q) (since Q may be small positive)
        double bestv = val(qL);
        qbest = qL;
        for (int k = qL; k <= qR; ++k)
        {
            double v = val(k);
            if (v > bestv)
            {
                bestv = v;
                qbest = k;
            }
        }
        for (int k = qR; k > qL; --k)
        {
            if (val(k) == 0.0 || val(k-1)*val(k) <= 0.0)
            {
                qbest = k;
                break;
            }
        }
    }
    q_indx = clamp_index(qbest, len);

    // S search: symmetrical to the right
    int sL = (int)r_indx;
    int sR = std::min((int)len - 1, (int)r_indx + qwin);
    int sbest = sL;
    if (r_positive)
    {
        // after positive R expect S to be a negative deflection -> local minimum
        double bestv = val(sL);
        sbest = sL;
        for (int k = sL; k <= sR; ++k)
        {
            double v = val(k);
            if (v < bestv)
            {
                bestv = v;
                sbest = k;
            }
        }
        for (int k = sL; k < sR; ++k)
        {
            if (val(k) == 0.0 || val(k)*val(k+1) <= 0.0)
            {
                sbest = k;
                break;
            }
        }
    }
    else
    {
        // after negative R expect a positive rebound -> local maximum
        double bestv = val(sL);
        sbest = sL;
        for (int k = sL; k <= sR; ++k)
        {
            double v = val(k);
            if (v > bestv)
            {
                bestv = v;
                sbest = k;
            }
        }
        for (int k = sL; k < sR; ++k)
        {
            if (val(k) == 0.0 || val(k)*val(k+1) <= 0.0)
            {
                sbest = k;
                break;
            }
        }
    }
    s_indx = clamp_index(sbest, len);

    // Guarantee ordering Q < R < S; if violated, try minor corrections
    if (!(q_indx < r_indx && r_indx < s_indx))
    {
        // try to enforce simple order by nudging
        if (q_indx >= r_indx)
        {
            if (r_indx > 0) q_indx = r_indx - 1;
            else q_indx = r_indx;
        }
        if (s_indx <= r_indx)
        {
            if (r_indx + 1 < len) s_indx = r_indx + 1;
            else s_indx = r_indx;
        }
    }

    // final safety clamp
    q_indx = std::min(q_indx, len - 1);
    r_indx = std::min(r_indx, len - 1);
    s_indx = std::min(s_indx, len - 1);
}

/*
     1) P ablak = Q-tól vissza 120–200 ms
     2) baseline = medián a teljes P ablakban
     3) durva P csúcs lokalizálása (max abs jel)
     4) P-csúcson finom extrémum keresés (lokális max/min)
     5) onset/offset = baseline-közeli metszések + derivált küszöb
*/

static inline size_t clamp_idx(long i, size_t len)
{
    if (i < 0) return 0;
    if ((size_t)i >= len) return len - 1;
    return (size_t)i;
}

int detect_p(double* ecg, size_t len, double fs, size_t q_indx, size_t& p_on_indx, size_t& p_indx, size_t& p_off_indx)
{
    p_on_indx = p_indx = p_off_indx = q_indx;

    if (!ecg || len < 10) return 1;

    // Helper for ms->samples
    auto ms2s = [&](double ms)->int
    {
        return std::max(1, (int)std::floor(ms * fs / 1000.0));
    };

    // Clinical P window
    int Pwin_min = ms2s(80);     // minimal: 80 ms Q előtt
    int Pwin_max = ms2s(220);    // maximális: 220 ms Q előtt
    // Hard limit: typical P width ~80–120 ms

    int R = (int)q_indx;
    int L = R - Pwin_max;
    int M = R - Pwin_min;

    if (M < 0) M = 0;
    if (L < 0) L = 0;
    if (R >= (int)len) R = (int)len - 1;

    if (L >= R) return 2;

    // Compute baseline (median) in the P window
    std::vector<double> buf;
    buf.reserve(R - L + 1);
    for (int i = L; i <= R; ++i) buf.push_back(ecg[i]);

    std::nth_element(buf.begin(), buf.begin() + buf.size()/2, buf.end());
    double baseline = buf[buf.size()/2];

    // baseline corrected accessor
    auto val = [&](int i)->double
    {
        i = std::max(0, std::min((int)len - 1, i));
        return ecg[i] - baseline;
    };

    // ---- 1) Durva P csúcs keresés: legnagyobb abs jel az ablakban ----
    int rough_peak = L;
    double bestAbs = 0.0;
    for (int i = L; i <= M; ++i)
    {
        double a = std::abs(val(i));
        if (a > bestAbs)
        {
            bestAbs = a;
            rough_peak = i;
        }
    }

    // ---- 2) Finom P csúcs: lokális maximum/minimum keresés körül ----
    int search_w = ms2s(20);  // ±20 ms körüli finom keresés
    int fpL = std::max(L, rough_peak - search_w);
    int fpR = std::min(M, rough_peak + search_w);

    int fine_peak = rough_peak;
    double best_ext = std::abs(val(rough_peak));

    // keresd a legjobb lokális extrémumot
    for (int i = fpL + 1; i <= fpR - 1; ++i)
    {
        double p = val(i);
        if ((p >= val(i-1) && p >= val(i+1)) ||
                (p <= val(i-1) && p <= val(i+1)))
        {
            double a = std::abs(p);
            if (a > best_ext)
            {
                best_ext = a;
                fine_peak = i;
            }
        }
    }

    p_indx = clamp_idx(fine_peak, len);

    // Determine polarity (positive or negative P)
    // bool p_positive = (val((int)p_indx) >= 0.0);

    // ---- 3) Onset detection: keresd baseline-közeli pontot balra ----
    int onL = L;
    int onR = (int)p_indx;

    double amp = std::abs(val((int)p_indx));
    double on_thr = amp * 0.08;       // 8% amplitude threshold
    double der_thr = amp * 0.015;     // derivative threshold ~1.5%/sample

    int onset = onL;
    bool found_on = false;

    for (int i = onR; i > onL + 2; --i)
    {
        double dv = std::abs(val(i) - val(i-1));
        if (std::abs(val(i)) < on_thr && dv < der_thr)
        {
            onset = i;
            found_on = true;
            break;
        }
    }

    if (!found_on)
    {
        // fallback: minimum abs value
        int best = onL;
        double bestA = std::abs(val(best));
        for (int i = onL + 1; i <= onR; ++i)
        {
            if (std::abs(val(i)) < bestA)
            {
                bestA = std::abs(val(i));
                best = i;
            }
        }
        onset = best;
    }

    p_on_indx = clamp_idx(onset, len);

    // ---- 4) Offset detection: baseline/derivative jobb oldalon ----
    int offL = (int)p_indx;
    int offR = M;

    double off_thr = amp * 0.08;      // 8% amplitude
    double der_thr2 = amp * 0.015;

    int offset = offR;
    bool found_off = false;

    for (int i = offL + 1; i < offR; ++i)
    {
        double dv = std::abs(val(i) - val(i-1));
        if (std::abs(val(i)) < off_thr && dv < der_thr2)
        {
            offset = i;
            found_off = true;
            break;
        }
    }

    if (!found_off)
    {
        // fallback: minimum abs region
        int best = offL;
        double bestA = std::abs(val(best));
        for (int i = offL + 1; i <= offR; ++i)
        {
            double a = std::abs(val(i));
            if (a < bestA)
            {
                bestA = a;
                best = i;
            }
        }
        offset = best;
    }

    p_off_indx = clamp_idx(offset, len);

    // ---- 6) Biztonság: p_on < p < p_off ----
    if (!(p_on_indx < p_indx && p_indx < p_off_indx))
    {
        p_on_indx = (p_indx > 1 ? p_indx - 1 : p_indx);
        p_off_indx = std::min(p_indx + 1, len - 1);
    }

    return 0;
}

static pqrst_positions detect_pqrst_positions(double** leads, unsigned int nr_ch, unsigned int r_idx, double sampling_rate, unsigned int nr_samples)
{
    //write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window1.bin", &leads[0][0], nr_samples, sampling_rate);
    double** leads_cpy = new double*[nr_ch];
    for (unsigned int ch = 0; ch < nr_ch; ++ch)
    {
        iir_filter_2nd_order bandpass_filter_;
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 1, sampling_rate, 1, 45);
        bandpass_filter_.init_history_values(leads[ch][0], sampling_rate);
        leads_cpy[ch] = new double[nr_samples];
        for (unsigned int i = 0; i < nr_samples; i++)
            leads_cpy[ch][i] = bandpass_filter_.filter(leads[ch][i]);
        for (int i = nr_samples - 1; i >= 0; i--)
            leads_cpy[ch][i] = bandpass_filter_.filter(leads_cpy[ch][i]);
    }
    //write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window2.bin", &leads[0][0], nr_samples, sampling_rate);

    pqrst_positions pos{};
    pos.r_idx = r_idx;

    // Q
    unsigned int window_q = (unsigned int)(0.4 * sampling_rate);
    unsigned int start_q = (r_idx > window_q) ? (r_idx - window_q) : 0;

    size_t q_indx, r_indx, s_indx;
    detect_qrs(&leads_cpy[0][start_q], 0.8 * sampling_rate, sampling_rate, r_idx - start_q, q_indx, r_indx, s_indx);
    pos.q_idx = q_indx + start_q;
    pos.r_idx = r_indx + start_q;
    pos.s_idx = s_indx + start_q;

//    pos.q_idx = find_min_from_to(leads_cpy[0], start_q, r_idx, nr_samples);
//
//    // S
//    unsigned int window_s = (unsigned int)(0.04 * sampling_rate);
//    unsigned int end_s = std::min(r_idx + window_s, nr_samples - 1);
//    pos.s_idx = find_min_from_to(leads_cpy[0], r_idx, end_s, nr_samples);

    // T
    //unsigned int t_search_start = 0.15 * sampling_rate;
    unsigned int t_search_start = pos.s_idx + 0.05 * sampling_rate;
    //cout << "--------------------------t_search_start " << t_search_start << " pos.s_idx " << pos.s_idx << " sampling_rate " << sampling_rate << endl;
    unsigned int t_search_end = std::min(r_idx + (unsigned int)(0.40 * sampling_rate), nr_samples - 1);

    size_t search_samples = t_search_end - t_search_start; ///(size_t)sampling_rate * 2;
    double* sig_window = new double[search_samples];
    create_ideal_signal(sig_window, leads_cpy, nr_ch, nr_samples, sampling_rate, t_search_start, t_search_end, t_search_start, t_search_end);
    if (false)
    {
        write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window1.bin", &leads_cpy[0][t_search_start], search_samples, 250);
        write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window2.bin", &leads_cpy[1][t_search_start], search_samples, 250);
        write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window3.bin", sig_window, search_samples, 250);
        const double* sig_window_[3];
        sig_window_[0] = sig_window;
        sig_window_[1] = &leads_cpy[0][t_search_start];
        sig_window_[2] = &leads_cpy[1][t_search_start];
        write_binmx_to_file("/media/sf_SharedFolder/sig_window.bin", sig_window_, 3, search_samples, 250);
    }

    pos.t_idx = find_max_from_to(sig_window, 1, search_samples, search_samples);
    pos.t_on_idx = find_min_with_baseline_correction(sig_window, 0, pos.t_idx);
    unsigned int t_offset_end = std::min(pos.t_idx + (int)(0.20 * sampling_rate), (int)search_samples - 1);
    pos.t_off_idx = find_min_with_baseline_correction(sig_window, pos.t_idx, t_offset_end);
    pos.t_idx += t_search_start;
    pos.t_on_idx += t_search_start;
    pos.t_off_idx += t_search_start;

    // P
    delete[] sig_window;
    unsigned int p_search_start = (r_idx > (unsigned int)(0.3 * sampling_rate)) ? (r_idx - (unsigned int)(0.3 * sampling_rate)) : 0;
    unsigned int p_search_end = pos.q_idx + 0.01 * sampling_rate;
    search_samples = p_search_end - p_search_start; ///(size_t)sampling_rate * 2;
    sig_window = new double[search_samples];


    if (false)
    {
        //memcpy(sig_window, &leads_cpy[1][p_search_start], search_samples * 8);
        write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window1.bin", &leads_cpy[0][p_search_start], search_samples, 250);
        write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window2.bin", &leads_cpy[1][p_search_start], search_samples, 250);
        write_binmx_to_file_1ch("/media/sf_SharedFolder/sig_window3.bin", sig_window, search_samples, 250);
        const double* sig_window_[3];
        sig_window_[0] = sig_window;
        sig_window_[1] = &leads_cpy[0][p_search_start];
        sig_window_[2] = &leads_cpy[1][p_search_start];
        write_binmx_to_file("/media/sf_SharedFolder/sig_window.bin", sig_window_, 3, p_search_end - p_search_start, 250);
    }

    if (false)
    {
        p_search_end -= 0.03 * sampling_rate;
        search_samples = p_search_end - p_search_start;
        create_ideal_signal(sig_window, leads_cpy, nr_ch, nr_samples, sampling_rate, p_search_start, p_search_end, p_search_start, p_search_end);
        pos.p_idx = find_max_from_to(sig_window, 1, search_samples, search_samples);
        pos.p_on_idx = find_min_with_baseline_correction(sig_window, 0, pos.p_idx);
        unsigned int p_offset_end = std::min(pos.p_idx + (int)(0.20 * sampling_rate), (int)search_samples - 1);
        pos.p_off_idx = find_min_with_baseline_correction(sig_window, pos.p_idx, p_offset_end);
        //pos.p_idx = (pos.p_idx + pos.p_off_idx) / 2;
        //pos.p_idx = find_max_from_to(sig_window, 1, search_samples, search_samples);
    }
    else
    {
        create_ideal_signal(sig_window, leads_cpy, nr_ch, nr_samples, sampling_rate, p_search_start, p_search_end, p_search_start, p_search_end);
        size_t p_on_indx, p_off_indx, p_indx;
        /*int res = */detect_p(sig_window, search_samples, sampling_rate, pos.q_idx - p_search_start, p_on_indx, p_indx, p_off_indx);
        pos.p_on_idx = p_on_indx;
        pos.p_idx = p_indx;
        pos.p_off_idx = p_off_indx;
    }

    pos.p_idx += p_search_start;
    pos.p_on_idx += p_search_start;
    pos.p_off_idx += p_search_start;

    //pos.p_idx = find_max_with_baseline_correction(leads_cpy[1], p_search_start + 1, p_search_end);
    /*
    pos.p_idx = find_max_with_baseline_correction(leads_cpy[1], p_search_start + 1, p_search_end);
    pos.p_on_idx = find_min_with_baseline_correction(leads_cpy[1], p_search_start, pos.p_idx);
    pos.p_off_idx = find_min_with_baseline_correction(leads_cpy[1], pos.p_idx, pos.q_idx);
    */

    for (unsigned int ch = 0; ch < nr_ch; ++ch)
        delete[] leads_cpy[ch];
    delete[] leads_cpy;
    delete[] sig_window;

    return pos;
}

static void fill_analysis_result(ecg_analysis_result& result, const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples, double sampling_rate, pqrst_positions& pos)
{
    result.pr_interval_ms = (pos.q_idx - pos.p_on_idx) / sampling_rate * 1000.0;
    result.p_wave_duration_ms = (pos.p_off_idx - pos.p_on_idx) / sampling_rate * 1000.0;
    unsigned int window_isoelectric_search = (unsigned int)(0.005 * sampling_rate);
    int qrs_start = find_isoelectric_point_before(ecg_signal[1], pos.q_idx - window_isoelectric_search, pos.q_idx);
    window_isoelectric_search = (unsigned int)(0.02 * sampling_rate);
    int qrs_end = find_isoelectric_point_after(ecg_signal[1], pos.s_idx, pos.s_idx + window_isoelectric_search);
    pos.q_idx = qrs_start;
    pos.s_idx = qrs_end;
    result.qrs_duration_ms = (qrs_end - qrs_start) / sampling_rate * 1000.0;
    result.qt_interval_ms = (pos.t_off_idx - pos.q_idx) / sampling_rate * 1000.0;

    double rr_sec = result.rr_interval_ms / 1000.0;
    result.qtc_interval_ms = (rr_sec > 0.0) ? result.qt_interval_ms / std::sqrt(rr_sec) : 0.0;
    result.t_wave_duration_ms = (pos.t_off_idx - pos.t_on_idx) / sampling_rate * 1000.0;

    unsigned int st_point = pos.q_idx + (unsigned int)(0.06 * sampling_rate);
    if (st_point >= nr_samples) st_point = nr_samples - 1;

    for (unsigned int ch = nr_ch; ch < 12; ++ch)
    {
        result.r_peak_amplitude_mV[ch] = 0.0;
        result.s_wave_amplitude_mV[ch] = 0.0;
        result.st_elevation_mV[ch] = 0.0;
        result.st_depression_mV[ch] = 0.0;
    }

    for (unsigned int ch = 0; ch < nr_ch && ch < 12; ++ch)
    {
        const double* lead = ecg_signal[ch];
        result.r_peak_amplitude_mV[ch] = lead[pos.r_idx];

        unsigned int window_s = (unsigned int)(0.04 * sampling_rate);
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

    double net_I = result.r_peak_amplitude_mV[0] - std::fabs(result.s_wave_amplitude_mV[0]);
    double net_aVF = result.r_peak_amplitude_mV[5] - std::fabs(result.s_wave_amplitude_mV[5]);
    result.frontal_plane_axis_deg = std::atan2(net_aVF, net_I) * 180.0 / M_PI;
    result.horizontal_plane_axis_deg = 0.0;
    result.pr_segment_ms = result.pr_interval_ms - result.p_wave_duration_ms;
    result.st_segment_ms = (pos.t_on_idx - pos.s_idx) / sampling_rate * 1000.0;;
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

static void check_sinus_rhythm(const vector<double>& rr_intervals, double rr_mean, double max_rr, double min_rr, double /**sampling_rate*/, ecg_analysis_result& result)
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
    for (double rr_ms : rr_intervals)
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
    result.pathologic_status[0] = 0;
    std::snprintf(result.status_message, sizeof(result.status_message), "OK");

    if (peak_indexes.size() < 1)
    {
        result.analysis_status = 1;
        std::snprintf(result.status_message, sizeof(result.status_message), "Not enough R peaks for RR calculation");
        return;
    }

    if (!ecg_signal || nr_ch < 1)
    {
        result.analysis_status = 2;
        std::snprintf(result.status_message, sizeof(result.status_message), "No channel data given");
        return;
    }

    if (nr_samples_per_ch < sampling_rate)
    {
        result.analysis_status = 3;
        std::snprintf(result.status_message, sizeof(result.status_message), "Not enough data samples given");
        return;
    }

    pqrst_positions pos = detect_pqrst_positions((double**)ecg_signal, nr_ch, peak_indexes[0], sampling_rate, nr_samples_per_ch);
    calculate_rr_statistics(peak_indexes, sampling_rate, result);
    fill_analysis_result(result, ecg_signal, nr_ch, nr_samples_per_ch, sampling_rate, pos);
    fill_annotations(annotations, pos);
}
#include <iomanip>

double median(const double* data, int len)
{
    if (len <= 0) return 0.0;

    std::vector<double> v(data, data + len);
    std::sort(v.begin(), v.end());

    if (len % 2 == 1)
        return v[len / 2];
    else
        return (v[len/2 - 1] + v[len/2]) * 0.5;
}

void search_isoel_bounds(double* signal, int len, double fs, int start_indx, double max_search_ms, double isoel_ms, double perimeter_ms, int& isoel_l, int& isoel_r, double isoel_tolerance = 0.01)
{
    //cout << endl << "------------------------------------------------------------" << endl;
    if (!signal || len <= 0)
        return;

    int nr_max_search_samples = max_search_ms * fs / 1000.0;
    int nr_isoel_samples = isoel_ms * fs / 1000.0;
    int nr_perimeter_samples = perimeter_ms * fs / 1000.0;

    if (nr_isoel_samples <= 0 || nr_max_search_samples <= 0)
        return;

    start_indx -= nr_perimeter_samples + nr_isoel_samples * 3;

    if (start_indx >= len) start_indx = len - 1;
    if (start_indx < 0) start_indx = 0;

    int stop_indx = start_indx - nr_max_search_samples;
    if (stop_indx < 0) stop_indx = 0;
    if (stop_indx >= len) stop_indx = len - 1;

    std::vector<double> deriv(len, 0.0);

    int limit = len - 2 * nr_isoel_samples;
    if (limit < 0) limit = 0;
    for (int i = 0; i < limit; ++i)
        for (int j = 0; j < nr_isoel_samples; ++j)
            deriv[i] += fabs(signal[i + j] - signal[i + j + nr_isoel_samples]);

    int nr_isoel_samples_m = nr_isoel_samples * 2;
    limit = len - 2 * nr_isoel_samples_m;
    if (limit < 0) limit = 0;
    for (int i = 0; i < limit; ++i)
        for (int j = 0; j < nr_isoel_samples_m; ++j)
            deriv[i] += fabs(signal[i + j] - signal[i + j + nr_isoel_samples_m]);

    nr_isoel_samples_m *= 3;
    limit = len - 2 * nr_isoel_samples_m;
    if (limit < 0) limit = 0;
    for (int i = 0; i < limit; ++i)
        for (int j = 0; j < nr_isoel_samples_m; ++j)
            deriv[i] += fabs(signal[i + j] - signal[i + j + nr_isoel_samples_m]);

    double min_val = 1e9, max_val = -1e9;
    for (int i = 0; i < len; ++i)
    {
        if (max_val < deriv[i]) max_val = deriv[i];
        if (min_val > deriv[i]) min_val = deriv[i];
    }

    double qrs_amp = max_val - min_val;
    double max_allowed_dev = qrs_amp * isoel_tolerance;

    std::vector<double> treshold(len, max_allowed_dev);
    //write_binmx_to_file("c:/Tamas/test002.bin", (const double**)&treshold, 1, len, fs);

    isoel_l = -1;
    //cout << "start_indx: " << start_indx << " stop_indx: " << stop_indx << endl;
    for (int i = start_indx; i >= stop_indx; --i)
    {
        //cout << "deriv[i]: " << deriv[i] << " max_allowed_dev: " << max_allowed_dev << " deriv[i] < max_allowed_dev " << (deriv[i] < max_allowed_dev) << endl;
        if (deriv[i] < max_allowed_dev /*&& fabs(signal[i] - signal[start_indx]) < max_allowed_dev * 7.0*/)
        {
            isoel_l = i;
            break;
        }
    }
    if (isoel_l > -1)
    {
        int mid = (start_indx - isoel_l) / 2 + isoel_l;
        if (mid >= len) mid = len - 1;

        bool ascending = signal[isoel_l] < signal[mid];

        for (int i = isoel_l; i < start_indx && i < len; ++i)
            if ((ascending && (signal[i] > signal[isoel_l] + max_allowed_dev)) || (!ascending && (signal[i] < signal[isoel_l] - max_allowed_dev)))
            {
                //isoel_l = i;
                break;
            }
    }

    isoel_r = -1;
    start_indx += nr_perimeter_samples * 2 + nr_isoel_samples * 3;
    stop_indx = start_indx + nr_max_search_samples;
    if (stop_indx < 0) stop_indx = 0;
    if (stop_indx >= len) stop_indx = len - 1;

    //cout << "start_indx: " << start_indx << " stop_indx: " << stop_indx << endl;
    for (int i = start_indx; i < stop_indx; ++i)
    {
        //cout << "deriv[i]: " << deriv[i] << " max_allowed_dev: " << max_allowed_dev << " deriv[i] < max_allowed_dev " << deriv[i] << max_allowed_dev << endl;
        if (deriv[i] < max_allowed_dev /*&& fabs(signal[i] - signal[start_indx]) < max_allowed_dev * 7.0*/)
        {
            //cout << i << ": fabs(signal[i]: " << fabs(signal[i]) << " max_allowed_dev: " << max_allowed_dev << endl;
            isoel_r = i;
            break;
        }
    }
    isoel_l += nr_isoel_samples * 3;
    isoel_r += nr_isoel_samples * 3;
    if (isoel_r > -1)
    {
        int mid = (isoel_r - start_indx) / 2 + start_indx;
        if (mid >= len) mid = len - 1;

        bool ascending = signal[isoel_r] < signal[mid];

        for (int i = isoel_r; i > start_indx; --i)
            if ((ascending && (signal[i] > signal[isoel_r] + max_allowed_dev)) || (!ascending && (signal[i] < signal[isoel_r] - max_allowed_dev)))
            {
                //isoel_r = i;
                break;
            }
    }
    //cout << "isoel_l: " << isoel_l << "  isoel_r: " << isoel_r << endl;

    //write_binmx_to_file("c:/Tamas/test003.bin", (const double**)&deriv, 1, len, fs);
}

void search_p_and_t_peaks(double* signal, /*double* deriv, */int len, double fs, double max_search_ms_l, double max_search_ms_r, int isoel_l, int isoel_r, /*double max_allowed_dev, */int& peak_p, int& peak_t)
{
    int nr_max_search_samples_l = max_search_ms_l * fs / 1000.0;
    int searc_stop_l = isoel_l - nr_max_search_samples_l;
    if (searc_stop_l < 0) searc_stop_l = 0;
    if (isoel_l >= len) isoel_l = len - 1;
    double maxval = -1e9;
    peak_p = -1;
    //cout << "XXXXXX  isoel_l " << isoel_l / 500.0 << " searc_stop_l " << searc_stop_l / 500.0 << endl;
    double base = signal[isoel_l];
    for (int i = isoel_l; i >= searc_stop_l; --i)
    {
        //cout << i << ": " << maxval << " / " << fabs(signal[i] - base) << endl;
        if (maxval < fabs(signal[i] - base))
        {
            maxval = fabs(signal[i] - base);
            peak_p = i;
        }
    }

    int nr_max_search_samples_r = max_search_ms_r * fs / 1000.0;
    int searc_stop_r = nr_max_search_samples_r + isoel_r;
    if (isoel_r < 0) isoel_r = 0;
    if (searc_stop_r >= len) searc_stop_r = len - 1;
    maxval = -1e9;
    peak_t = -1;
    base = signal[isoel_r];
    for (int i = isoel_r; i < searc_stop_r; ++i)
    {
        if (maxval < fabs(signal[i] - base))
        {
            maxval = fabs(signal[i] - base);
            peak_t = i;
        }
    }
}

//    int isoel_count = 0;
//    int isoel_indx = -1;
//    for (int i = start_indx; i >= stop_indx; --i)
//    {
//        if (fabs(avg - signal[i]) < isoel_tolerance)
//        {
//            ++isoel_count;
//            if (!isoel_indx)
//                isoel_indx = i;
//        }
//        else if (isoel_count)
//            break;
//        if (isoel_count == nr_isoel_samples)
//            break;
//    }
//    if (isoel_count == nr_isoel_samples)

//int search_isoel_bounds(double* signal, int len, double fs, int start_indx, double max_search_ms, double isoel_ms, double isoel_tolerance = 0.006)
//{
//    int nr_max_search_samples = max_search_ms * fs / 1000.0;
//    int nr_isoel_samples = isoel_ms * fs / 1000.0;
//
//    double min_val = 1e9, max_val = -1e9, avg = 0;
//
//    for (int i = 0; i < len; ++i)
//    {
//        avg += signal[i];
//        if (max_val < signal[i]) max_val = signal[i];
//        if (min_val > signal[i]) min_val = signal[i];
//    }
//    avg /= len;
//
//    double qrs_amp = max_val - min_val;
//    double max_allowed_dev = qrs_amp * isoel_tolerance / fs;
//
//    int stop_indx = start_indx - nr_max_search_samples;
//    if (stop_indx < nr_isoel_samples)
//        stop_indx = nr_isoel_samples;
//
//    double* test_sig = new double[len];
//    //memset(test_sig, 0, len * sizeof(double));
//    for (int i = 0; i < len; ++i)
//        test_sig[i] = 0;
//    // Visszafelé keresünk az R csúcstól
//    cout << "start_indx: " << start_indx << " stop_indx: " << stop_indx << " nr_isoel_samples: " << nr_isoel_samples << endl;
//    for (int i = 0; i < len - 1; ++i)
//    {
//        double slope = fabs(signal[i] - signal[i + 1]);
//        test_sig[i] = slope;
//        //test_sig[i] *= qrs_amp / (double)nr_isoel_samples;
//    }
//    write_binmx_to_file("c:/Tamas/test003.bin", (const double**)&test_sig, 1, len, fs);
//    delete[] test_sig;
//    return -1; // nem található isoelektromos szakasz
//}

void analyse_ecg(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, ecg_analysis_result& result, unsigned int analysis_ch_indx = 0)
{
    annotations.clear();
    result.analysis_status = 0;
    result.pathologic_status[0] = 0;
    std::snprintf(result.status_message, sizeof(result.status_message), "OK");

    if (peak_indexes.size() < 1)
    {
        result.analysis_status = 1;
        std::snprintf(result.status_message, sizeof(result.status_message), "No R peaks found");
        return;
    }

    if (!ecg_signal || nr_ch < 1)
    {
        result.analysis_status = 2;
        std::snprintf(result.status_message, sizeof(result.status_message), "No channel data given");
        return;
    }

    if (nr_samples_per_ch < sampling_rate)
    {
        result.analysis_status = 3;
        std::snprintf(result.status_message, sizeof(result.status_message), "Not enough data samples given");
        return;
    }

    if (nr_ch <= analysis_ch_indx)
    {
        result.analysis_status = 4;
        std::snprintf(result.status_message, sizeof(result.status_message), "Invalid channel index");
        return;
    }

    unsigned int ch = analysis_ch_indx;
    iir_filter_2nd_order bandpass_filter_;
    create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 1, sampling_rate, 0.5, 35);
    bandpass_filter_.init_history_values(ecg_signal[ch][0], sampling_rate);
    double* lead_cpy = new double[nr_samples_per_ch];
    for (unsigned int i = 0; i < nr_samples_per_ch; i++)
        lead_cpy[i] = bandpass_filter_.filter(ecg_signal[ch][i]);
    bandpass_filter_.init_history_values(lead_cpy[nr_samples_per_ch - 1], sampling_rate);
    for (int i = nr_samples_per_ch - 1; i >= 0; i--)
        lead_cpy[i] = bandpass_filter_.filter(lead_cpy[i]);
    //write_binmx_to_file("c:/Tamas/test001.bin", (const double**)&lead_cpy, 1, nr_samples_per_ch, sampling_rate);
//    std::sort(lead_cpy, lead_cpy + nr_samples_per_ch);
//    for (unsigned int i = 0; i < nr_samples_per_ch; i++)
//        lead_cpy[i] = lead_cpy[nr_samples_per_ch / 2];
//    write_binmx_to_file("c:/Tamas/test002.bin", (const double**)&lead_cpy, 1, nr_samples_per_ch, sampling_rate);

    int isoel_start = -1, isoel_r = -1;
    search_isoel_bounds(lead_cpy, nr_samples_per_ch, sampling_rate, peak_indexes[0], 110, 2, 35, isoel_start, isoel_r, 0.05);
    int peak_p, peak_t;
    search_p_and_t_peaks(lead_cpy, nr_samples_per_ch, sampling_rate, 300, 300, isoel_start, isoel_r, peak_p, peak_t);

    int isoel_p_l = -1, isoel_p_r = -1;
    search_isoel_bounds(lead_cpy, nr_samples_per_ch, sampling_rate, peak_p, 170, 4, 35, isoel_p_l, isoel_p_r, 0.05);

    int isoel_t_l = -1, isoel_t_r = -1;
    search_isoel_bounds(lead_cpy, nr_samples_per_ch, sampling_rate, peak_t, 220, 2, 35, isoel_t_l, isoel_t_r, 0.02);

    //cout << "peak_p: " << peak_p / 500.0 << " peak_t: " << peak_t / 500.0 << endl;
    annotations.push_back({});
    annotations[0].p[0] = isoel_p_l;
    annotations[0].p[1] = peak_p;
    annotations[0].p[2] = isoel_p_r;

    annotations[0].r[0] = isoel_start;
    annotations[0].r[1] = (isoel_r + isoel_start) / 2;
    annotations[0].r[2] = isoel_r;

    annotations[0].t[0] = isoel_t_l;
    annotations[0].t[1] = peak_t;
    annotations[0].t[2] = isoel_t_r;

    //cout << "isoel_start: " << isoel_start / 500.0 << " isoel_r: " << isoel_r / 500.0 << endl;

    delete[] lead_cpy;
}

ecg_analysis_result analyse_ecg_detect_peaks(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, std::vector<pqrst_indxes>& annotations, std::vector<unsigned int>* peak_indexes, std::string mode, int analysis_ch_indx)
{
    if (analysis_ch_indx < 0) analysis_ch_indx = 0;
    if (analysis_ch_indx >= (int)nr_channels) analysis_ch_indx = nr_channels - 1;
    std::vector<double> peak_signal(nr_samples_per_channel), filt_signal(nr_samples_per_channel), threshold_signal(nr_samples_per_channel);
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

    //detector.detect_multichannel(data, nr_channels, nr_samples_per_channel, peak_signal.data(), filt_signal.data(), threshold_signal.data(), peak_indexes);
    detector.detect(data[analysis_ch_indx], nr_samples_per_channel, peak_signal.data(), filt_signal.data(), threshold_signal.data(), peak_indexes);

    ecg_analysis_result result;
    //analyse_ecg_multichannel(data, (unsigned int)nr_channels, nr_samples_per_channel, sampling_rate, *peak_indexes, annotations, result);
    analyse_ecg(data, (unsigned int)nr_channels, nr_samples_per_channel, sampling_rate, *peak_indexes, annotations, result, analysis_ch_indx);

    if (local_peak_indexes_used)
        delete peak_indexes;
    return result;
}
