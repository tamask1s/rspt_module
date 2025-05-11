/***************************************************************************
* Copyright 2024 Tamas Levente Kis - tamkis@gmail.com                      *
*                                                                          *
* Licensed under the Apache License, Version 2.0 (the "License");          *
* you may not use this file except in compliance with the License.         *
* You may obtain a copy of the License at                                  *
*                                                                          *
*     http://www.apache.org/licenses/LICENSE-2.0                           *
*                                                                          *
* Unless required by applicable law or agreed to in writing, software      *
* distributed under the License is distributed on an "AS IS" BASIS,        *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
* See the License for the specific language governing permissions and      *
* limitations under the License.                                           *
***************************************************************************/

/**
 * @brief Class for detecting peaks in a signal.
 *
 * This class implements a peak detection algorithm using a combination of
 * bandpass filtering, integration, and adaptive thresholding.
 * It is designed to work with real-time sampled signals, identifying peaks
 * based on a dynamic threshold derived from the signal itself.
 *
 * The detection logic includes:
 * - Bandpass filtering to isolate relevant frequency components.
 * - Integration to smooth out short-term fluctuations.
 * - Adaptive thresholding to dynamically adjust sensitivity.
 *
 * Once a peak is detected, a marker value is returned, otherwise, zero is returned.
 * Peaks are detected with a delay of 220 ms.
 */
class peak_detector
{
    iir_filter_4th_order bandpass_filter_;
    iir_filter_2nd_order integrative_filter_;
    iir_filter_2nd_order threshold_filter_;

    double previous_peak_amplitude_ = 0;
    double previous_sig_val_ = 0;
    bool searching_for_peaks_ = false;
    int samples_after_peak_count_ = 0;
    int sample_indx_ = 0;

    const double sampling_rate_;
    const double marker_val_;
    const double previous_peak_reference_ratio_ = 0.5;
    const double previous_peak_reference_attenuation_ = 25;
    const double peak_attenuation_;
    const double threshold_ratio_ = 1.5;
    const int nr_slope_samples_;

public:
    /**
     * @brief Constructor for the peak detector.
     *
     * Initializes the detector with a given sampling rate and optional marker value.
     * It also configures the internal filters using predefined Butterworth filter settings.
     *
     * @param sampling_rate The sampling rate of the input signal in Hz.
     * @param marker_val The value returned when a peak is detected (default: 1).
     */
    peak_detector(double sampling_rate, double marker_val = 1)
        : sampling_rate_(sampling_rate),
          marker_val_(marker_val),
          peak_attenuation_(1.0 / (1.0 + previous_peak_reference_attenuation_ / sampling_rate)),
          nr_slope_samples_((100.0 * sampling_rate) / 1000.0)
    {
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 2, sampling_rate, 10, 20);
        create_filter_iir(integrative_filter_.d, integrative_filter_.n, butterworth, low_pass, 2, sampling_rate, 3, 0);
        create_filter_iir(threshold_filter_.d, threshold_filter_.n, butterworth, low_pass, 2, sampling_rate, 0.15, 0);
    }

    /**
     * @brief Processes a new sample and detects peaks.
     *
     * This function applies the filtering and peak detection logic to an incoming sample.
     * The algorithm tracks changes in the signal and determines if a peak has occurred
     * based on an adaptive threshold.
     *
     * @param new_sample The new sample from the signal.
     * @return marker_val_ if a peak is detected, otherwise 0.
     */
    inline double detect(double new_sample, double* peak_sample = 0, double* threshold_sample = 0)
    {
        if (!sample_indx_++)
            bandpass_filter_.init_history_values(new_sample, sampling_rate_);

        double sig_val = bandpass_filter_.filter(new_sample);
        sig_val = integrative_filter_.filter(sig_val * sig_val);
        double threshold = threshold_filter_.filter(sig_val);
        (peak_sample) ? ((*peak_sample) = sig_val) : sample_indx_ = sample_indx_;
        (threshold_sample) ? ((*threshold_sample) = threshold) : sample_indx_ = sample_indx_;

        if (searching_for_peaks_ && (sig_val > threshold * threshold_ratio_) && (previous_sig_val_ > sig_val))
        {
            if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
            {
                previous_peak_amplitude_ = previous_sig_val_;
                samples_after_peak_count_ = 1;
                searching_for_peaks_ = false;
            }
            else
                previous_peak_amplitude_ *= peak_attenuation_;
        }
        else if (previous_sig_val_ < sig_val)
        {
            searching_for_peaks_ = true;
            samples_after_peak_count_ = 0;
        }

        previous_sig_val_ = sig_val;

        if (samples_after_peak_count_)
            ++samples_after_peak_count_;

        if (samples_after_peak_count_ == nr_slope_samples_)
        {
            samples_after_peak_count_ = 0;
            return ((marker_val_ == -1.0) ? sig_val : marker_val_);
        }
        return 0;
    }
};

//class peak_detector_1st_order
//{
//    iir_filter_2nd_order bandpass_filter_;
//    iir_filter_1st_order integrative_filter_;
//    iir_filter_2nd_order threshold_filter_;
//
//    double previous_peak_amplitude_ = 0;
//    double previous_sig_val_ = 0;
//    bool searching_for_peaks_ = false;
//    int samples_after_peak_count_ = 0;
//    int sample_indx_ = 0;
//
//    const double sampling_rate_;
//    const double marker_val_;
//    const double previous_peak_reference_ratio_ = 0.5;
//    const double previous_peak_reference_attenuation_ = 25;
//    const double peak_attenuation_;
//    const double threshold_ratio_ = 1.5;
//    const int nr_slope_samples_;
//
//public:
//    /**
//     * @brief Constructor for the peak detector.
//     *
//     * Initializes the detector with a given sampling rate and optional marker value.
//     * It also configures the internal filters using predefined Butterworth filter settings.
//     *
//     * @param sampling_rate The sampling rate of the input signal in Hz.
//     * @param marker_val The value returned when a peak is detected (default: 1).
//     */
//    peak_detector_1st_order(double sampling_rate, double marker_val = 1)
//        : sampling_rate_(sampling_rate),
//          marker_val_(marker_val),
//          peak_attenuation_(1.0 / (1.0 + previous_peak_reference_attenuation_ / sampling_rate)),
//          nr_slope_samples_((100.0 * sampling_rate) / 1000.0)
//    {
//        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 1, sampling_rate, 10, 20);
//        create_filter_iir(integrative_filter_.d, integrative_filter_.n, butterworth, low_pass, 1, sampling_rate, 3, 0);
//        create_filter_iir(threshold_filter_.d, threshold_filter_.n, butterworth, low_pass, 2, sampling_rate, 0.15, 0);
//    }
//
//    /**
//     * @brief Processes a new sample and detects peaks.
//     *
//     * This function applies the filtering and peak detection logic to an incoming sample.
//     * The algorithm tracks changes in the signal and determines if a peak has occurred
//     * based on an adaptive threshold.
//     *
//     * @param new_sample The new sample from the signal.
//     * @return marker_val_ if a peak is detected, otherwise 0.
//     */
//    inline double detect(double new_sample, double* peak_sample = 0, double* threshold_sample = 0)
//    {
//        if (!sample_indx_++)
//            bandpass_filter_.init_history_values(new_sample, sampling_rate_);
//
//        double sig_val = bandpass_filter_.filter(new_sample);
//        sig_val = integrative_filter_.filter(sig_val * sig_val);
//        double threshold = threshold_filter_.filter(sig_val);
//        (peak_sample) ? ((*peak_sample) = sig_val) : sample_indx_ = sample_indx_;
//        (threshold_sample) ? ((*threshold_sample) = threshold) : sample_indx_ = sample_indx_;
//
//        if (searching_for_peaks_ && (sig_val > threshold * threshold_ratio_) && (previous_sig_val_ > sig_val))
//        {
//            if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
//            {
//                previous_peak_amplitude_ = previous_sig_val_;
//                samples_after_peak_count_ = 1;
//                searching_for_peaks_ = false;
//            }
//            else
//                previous_peak_amplitude_ *= peak_attenuation_;
//        }
//        else if (previous_sig_val_ < sig_val)
//        {
//            searching_for_peaks_ = true;
//            samples_after_peak_count_ = 0;
//        }
//
//        previous_sig_val_ = sig_val;
//
//        if (samples_after_peak_count_)
//            ++samples_after_peak_count_;
//
//        if (samples_after_peak_count_ == nr_slope_samples_)
//        {
//            samples_after_peak_count_ = 0;
//            return ((marker_val_ == -1.0) ? sig_val : marker_val_);
//        }
//        return 0;
//    }
//};

class peak_detector_offline
{
    iir_filter_2nd_order bandpass_filter_;
    iir_filter_1st_order integrative_filter_;
    iir_filter_1st_order baseline_filter_;
    iir_filter_1st_order threshold_filter_;

    double previous_peak_amplitude_ = 0;
    double previous_sig_val_ = 0;
    bool searching_for_peaks_ = false;
    int samples_after_peak_count_ = 0;
    int sample_indx_ = 0;

    const double sampling_rate_;
    const double marker_val_;
    static constexpr double previous_peak_reference_ratio_default_ = 0.0;
    double previous_peak_reference_ratio_ = 0.0;
    double previous_peak_reference_attenuation_ = 150;
    double peak_attenuation_;
    double threshold_ratio_ = 1.5;
    double peak_correction_radius_ms_ = 10.0;
    int nr_slope_samples_;
    int nr_slope_msecs_ = 50;

public:
    enum class Mode { def, high_sensitivity, high_ppv };
    /**
     * @brief Constructor for the peak detector.
     *
     * Initializes the detector with a given sampling rate and optional marker value.
     * It also configures the internal filters using predefined Butterworth filter settings.
     *
     * @param sampling_rate The sampling rate of the input signal in Hz.
     * @param marker_val The value returned when a peak is detected (default: 1).
     */
    peak_detector_offline(double sampling_rate, double marker_val = 1)
        : sampling_rate_(sampling_rate),
          marker_val_(marker_val)
    {
        create_filter_iir(bandpass_filter_.d, bandpass_filter_.n, butterworth, band_pass, 1, sampling_rate, 19, 25);
        create_filter_iir(integrative_filter_.d, integrative_filter_.n, butterworth, low_pass, 1, sampling_rate, 2.7, 0);
        create_filter_iir(baseline_filter_.d, baseline_filter_.n, butterworth, low_pass, 1, sampling_rate, 0.6, 0);
        create_filter_iir(threshold_filter_.d, threshold_filter_.n, butterworth, low_pass, 1, sampling_rate, 0.55, 0);
        set_mode(peak_detector_offline::Mode::def);
    }

    void set_mode(peak_detector_offline::Mode mode)
    {
        if (mode == Mode::high_sensitivity)
        {
            previous_peak_reference_ratio_ = 0.0;
            previous_peak_reference_attenuation_ = 250;
            nr_slope_msecs_ = 15;
            threshold_ratio_ = 0.3;
        }
        else if (mode == Mode::high_ppv)
        {
            previous_peak_reference_ratio_ = 0.35;
            previous_peak_reference_attenuation_ = 50;
            nr_slope_msecs_ = 150;
            threshold_ratio_ = 2.7;
        }
        else
        {
            previous_peak_reference_ratio_ = 0.0;
            previous_peak_reference_attenuation_ = 150;
            nr_slope_msecs_ = 50;
            threshold_ratio_ = 1.5;
        }
        peak_attenuation_ = 1.0 / (1.0 + previous_peak_reference_attenuation_ / sampling_rate_);
        nr_slope_samples_ = (nr_slope_msecs_ * sampling_rate_) / 1000.0;
    }

    inline int detect(double* ecg_signal, unsigned int len, double* peak_signal, double* filt_signal, double* threshold_signal, std::vector<unsigned int>* peak_indexes = 0, double r_reference_value = 0, double previous_peak_reference_ratio = previous_peak_reference_ratio_default_)
    {
        bandpass_filter_.init_history_values(ecg_signal[0], sampling_rate_);
        baseline_filter_.init_history_values(ecg_signal[0], sampling_rate_);
        integrative_filter_.reset();
        threshold_filter_.reset();

        double* baseline = new double[len];
        for (unsigned int i = 0; i < len; ++i)
        {
            baseline[i] = baseline_filter_.filter(ecg_signal[i]);
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        }
        for (int i = len - 1; i > -1; --i)
        {
            baseline[i] = baseline_filter_.filter(baseline[i]);
            filt_signal[i] = bandpass_filter_.filter(filt_signal[i]);
        }
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i] * filt_signal[i]);
        for (int i = len - 1; i > -1; --i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            threshold_signal[i] = threshold_filter_.filter(filt_signal[i]);
        for (int i = len - 1; i > -1; --i)
            threshold_signal[i] = threshold_filter_.filter(threshold_signal[i]);

        for (unsigned int i = 0; i < len; ++i)
        {
            if (searching_for_peaks_ && (filt_signal[i] > threshold_signal[i] * threshold_ratio_) && (previous_sig_val_ > filt_signal[i]))
            {
                if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio))
                {
                    previous_peak_amplitude_ = previous_sig_val_;
                    samples_after_peak_count_ = 1;
                    searching_for_peaks_ = false;
                }
                else
                    previous_peak_amplitude_ *= peak_attenuation_;
                if (r_reference_value)
                    previous_peak_amplitude_ = r_reference_value;
            }
            else if (previous_sig_val_ < filt_signal[i])
            {
                searching_for_peaks_ = true;
                samples_after_peak_count_ = 0;
            }

            previous_sig_val_ = filt_signal[i];

            if (samples_after_peak_count_)
                ++samples_after_peak_count_;

            if (samples_after_peak_count_ == nr_slope_samples_)
            {
                samples_after_peak_count_ = 0;
                peak_signal[i] = ((marker_val_ == -1.0) ? filt_signal[i] : marker_val_);
            }
            else
                peak_signal[i] = 0;
        }
        double average_r = 0;
        unsigned int nr_peaks = 0;
        for (unsigned int i = 0; i < nr_slope_samples_; ++i)
            peak_signal[i] = 0;
        for (unsigned int i = nr_slope_samples_; i < len; ++i)
            if (peak_signal[i])
            {
                double val = peak_signal[i];
                peak_signal[i] = 0;
                peak_signal[i - nr_slope_samples_ + 1] = val;
                average_r += filt_signal[i - nr_slope_samples_ + 1];
                ++nr_peaks;
            }
        average_r /= (double)nr_peaks;
        const int radius = (peak_correction_radius_ms_ * sampling_rate_) / 1000.0;
        unsigned int nr_min_peaks = 0;
        for (unsigned int i = radius; i < len - radius; ++i)
            if (peak_signal[i])
            {
                unsigned int maxindx = 0, minindx = 0;
                double maxval = -2000000, minval = 2000000;
                for (int j = -radius; j < radius; ++j)
                {
                    if (maxval < ecg_signal[i + j] - baseline[i + j])
                    {
                        maxval = ecg_signal[i + j] - baseline[i + j];
                        maxindx = i + j;
                    }
                    if (minval > ecg_signal[i + j] - baseline[i + j])
                    {
                        minval = ecg_signal[i + j] - baseline[i + j];
                        minindx = i + j;
                    }
                }
                double peakval = peak_signal[i];
                peak_signal[i] = 0;
                if (maxval > -minval)
                {
                    nr_min_peaks++;
                    peak_signal[maxindx] = peakval;
                }
                else
                    peak_signal[minindx] = peakval;
            }
        if (peak_indexes)
        {
            peak_indexes->resize(nr_peaks);
            nr_peaks = 0;
            for (unsigned int i = 0; i < len; ++i)
                if (peak_signal[i])
                    (*peak_indexes)[nr_peaks++] = i;
        }
        delete[] baseline;
        return (nr_peaks / 2) > nr_min_peaks ? nr_peaks : (-nr_peaks);
    }

//    inline void cleanup_index_array(std::vector<unsigned int>& final_peak_indexes, std::vector<std::vector<unsigned int>>& peak_index_array, unsigned int tolerance, unsigned int nr_channels)
//    {
//        final_peak_indexes.clear();
//        if (nr_channels == 0) return;
//
//        const std::vector<unsigned int>& ref_peaks = peak_index_array[0];
//
//        for (unsigned int ref_index : ref_peaks)
//        {
//            uint32_t invalid_count = 0;
//
//            // Ellenőrizzük, hogy a többi csatornában van-e megfelelő index
//            for (unsigned int ch = 1; ch < nr_channels; ++ch)
//            {
//                const std::vector<unsigned int>& channel_peaks = peak_index_array[ch];
//
//                // Keresés a tolerance-en belüli indexekre
//                bool found_match = std::any_of(channel_peaks.begin(), channel_peaks.end(), [&](unsigned int idx)
//                {
//                    return std::abs((int)idx - (int)ref_index) <= (int)tolerance;
//                });
//
//                if (!found_match)
//                {
//                    ++invalid_count;
//                    break;
//                }
//            }
//
//            if (invalid_count < nr_channels / 2)
//                final_peak_indexes.push_back(ref_index);
//        }
//    }

    inline void detect_multichannel(double** ecg_signal_multich, unsigned int nr_channels, unsigned int len, double* peak_signal, double* filt_signal, double* threshold_signal, std::vector<unsigned int>* peak_indexes = 0, double r_reference_value = 0, double previous_peak_reference_ratio = previous_peak_reference_ratio_default_)
    {
        double sign_array[nr_channels];
        for (unsigned int i = 0; i < nr_channels; ++i)
            sign_array[i] = 1;
        for (unsigned int i = 0; i < nr_channels; ++i)
        {
            sign_array[i] = detect(ecg_signal_multich[i], len, peak_signal, filt_signal, threshold_signal, 0, r_reference_value, previous_peak_reference_ratio);
            if (sign_array[i] < 0)
                sign_array[i] = -1;
            else
                sign_array[i] = 0;
        }

        double* ecg_signal = new double[len];
        for (size_t i = 0; i < len; ++i)
        {
            ecg_signal[i] = ecg_signal_multich[0][i];
            for (size_t ch = 1; ch < nr_channels; ++ch)
                ecg_signal[i] += ecg_signal_multich[ch][i] * sign_array[ch];
            ecg_signal[i] /= (double)nr_channels;
        }

//        std::vector<std::vector<unsigned int>> peak_detector_offline::peak_index_array;

        detect(ecg_signal, len, peak_signal, filt_signal, threshold_signal, peak_indexes, r_reference_value, previous_peak_reference_ratio);
//        peak_index_array.resize(nr_channels);
//        set_mode(Mode::high_sensitivity);
//        for (unsigned int i = 0; i < nr_channels; ++i)
//            detect(ecg_signal_multich[i], len, peak_signal, filt_signal, threshold_signal, &peak_index_array[i], r_reference_value, previous_peak_reference_ratio);
//        const unsigned int tolerance = (50 * sampling_rate_) / 1000.0;
//        //*peak_indexes = peak_index_array[0];
//        cleanup_index_array(*peak_indexes, peak_index_array, tolerance, nr_channels);
        delete[] ecg_signal;
    }

    /*inline void detect9985_9978_9982(double* ecg_signal, unsigned int len, double* peak_signal, double* filt_signal, double* threshold_signal, std::vector<unsigned int>* peak_indexes = 0, double r_reference_value = 0, double previous_peak_reference_ratio = 0.5)
    {
        const double previous_peak_reference_ratio_ = previous_peak_reference_ratio;
        bandpass_filter_.init_history_values(ecg_signal[0], sampling_rate_);
        baseline_filter_.init_history_values(ecg_signal[0], sampling_rate_);

        double* baseline = new double[len];
        for (unsigned int i = 0; i < len; ++i)
            baseline[i] = baseline_filter_.filter(ecg_signal[i]);
        for (int i = len - 1; i > -1; --i)
            baseline[i] = baseline_filter_.filter(baseline[i]);
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        for (int i = len - 1; i > -1; --i)
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i] * filt_signal[i]);
        for (int i = len - 1; i > -1; --i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            threshold_signal[i] = threshold_filter_.filter(filt_signal[i]);
        for (int i = len - 1; i > -1; --i)
            threshold_signal[i] = threshold_filter_.filter(threshold_signal[i]);

        for (unsigned int i = 0; i < len; ++i)
        {
            if (searching_for_peaks_ && (filt_signal[i] > threshold_signal[i] * threshold_ratio_) && (previous_sig_val_ > filt_signal[i]))
            {
                if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
                {
                    previous_peak_amplitude_ = previous_sig_val_;
                    samples_after_peak_count_ = 1;
                    searching_for_peaks_ = false;
                }
                else
                    previous_peak_amplitude_ *= peak_attenuation_;
                if (r_reference_value)
                    previous_peak_amplitude_ = r_reference_value;
            }
            else if (previous_sig_val_ < filt_signal[i])
            {
                searching_for_peaks_ = true;
                samples_after_peak_count_ = 0;
            }

            previous_sig_val_ = filt_signal[i];

            if (samples_after_peak_count_)
                ++samples_after_peak_count_;

            if (samples_after_peak_count_ == nr_slope_samples_)
            {
                samples_after_peak_count_ = 0;
                peak_signal[i] = ((marker_val_ == -1.0) ? filt_signal[i] : marker_val_);
            }
            else
                peak_signal[i] = 0;
        }
        unsigned int nr_peaks = 0;
        for (unsigned int i = nr_slope_samples_; i < len; ++i)
            if (peak_signal[i])
            {
                double val = peak_signal[i];
                peak_signal[i] = 0;
                peak_signal[i - nr_slope_samples_ + 1] = val;
                ++nr_peaks;
            }
        const int radius = (10.0 * sampling_rate_) / 1000.0;
        for (unsigned int i = radius; i < len - radius; ++i)
            if (peak_signal[i])
            {
                unsigned int maxindx = 0, minindx = 0;
                double maxval = -2000000, minval = 2000000;
                for (int j = -radius; j < radius; ++j)
                {
                    if (maxval < ecg_signal[i + j] - baseline[i + j])
                    {
                        maxval = ecg_signal[i + j] - baseline[i + j];
                        maxindx = i + j;
                    }
                    if (minval > ecg_signal[i + j] - baseline[i + j])
                    {
                        minval = ecg_signal[i + j] - baseline[i + j];
                        minindx = i + j;
                    }
                }
                double peakval = peak_signal[i];
                peak_signal[i] = 0;
                if (maxval > -minval)
                    peak_signal[maxindx] = peakval;
                else
                    peak_signal[minindx] = peakval;
            }
        double average_r = 0;
        for (unsigned int i = nr_slope_samples_; i < len; ++i)
            if (peak_signal[i])
            {
                for (unsigned int j = 0; j < nr_slope_samples_ - 1; ++j)
                    if (peak_signal[i - j] < peak_signal[i - j - 1])
                    {
                        peak_signal[i] = 0;
                        --nr_peaks;
                    }
                if (peak_signal[i])
                {
                    average_r += filt_signal[i];
                }
            }
        average_r /= (double)nr_peaks;
        if (peak_indexes)
        {
            peak_indexes->resize(nr_peaks);
            nr_peaks = 0;
            for (unsigned int i = 0; i < len; ++i)
                if (peak_signal[i])
                    (*peak_indexes)[nr_peaks++] = i;
        }
        delete[] baseline;
        if (!r_reference_value)
            detect(ecg_signal, len, peak_signal, filt_signal, threshold_signal, peak_indexes, average_r, 0.2);
    }*/

    /*inline void detect_fw(double* ecg_signal, unsigned int len, double* peak_signal, double* filt_signal, double* threshold_signal)
    {
        bandpass_filter_.init_history_values(ecg_signal[0], sampling_rate_);

        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = bandpass_filter_.filter(ecg_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            filt_signal[i] = integrative_filter_.filter(filt_signal[i] * filt_signal[i]);
        for (unsigned int i = 0; i < len; ++i)
            threshold_signal[i] = threshold_filter_.filter(filt_signal[i]);

        for (unsigned int i = 0; i < len; ++i)
        {
            if (searching_for_peaks_ && (filt_signal[i] > threshold_signal[i] * threshold_ratio_) && (previous_sig_val_ > filt_signal[i]))
            {
                if ((previous_peak_amplitude_ == 0) || (previous_sig_val_ > previous_peak_amplitude_ * previous_peak_reference_ratio_))
                {
                    previous_peak_amplitude_ = previous_sig_val_;
                    samples_after_peak_count_ = 1;
                    searching_for_peaks_ = false;
                }
                else
                    previous_peak_amplitude_ *= peak_attenuation_;
            }
            else if (previous_sig_val_ < filt_signal[i])
            {
                searching_for_peaks_ = true;
                samples_after_peak_count_ = 0;
            }

            previous_sig_val_ = filt_signal[i];

            if (samples_after_peak_count_)
                ++samples_after_peak_count_;

            if (samples_after_peak_count_ == nr_slope_samples_)
            {
                samples_after_peak_count_ = 0;
                peak_signal[i] = ((marker_val_ == -1.0) ? filt_signal[i] : marker_val_);
            }
            else
                peak_signal[i] = 0;
        }
    }*/
};
