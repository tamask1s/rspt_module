
struct ecg_analysis_result
{
    double heart_rate_bpm;
    double rr_interval_ms;
    double rr_variation_ms;
    double pr_interval_ms;
    double pr_segment_ms;
    double qrs_duration_ms;
    double qt_interval_ms;
    double qtc_interval_ms;
    double st_segment_ms;

    double p_wave_duration_ms;
    double t_wave_duration_ms;

    double r_peak_amplitude_mV[12];
    double s_wave_amplitude_mV[12];

    double st_elevation_mV[12];
    double st_depression_mV[12];

    double frontal_plane_axis_deg;
    double horizontal_plane_axis_deg;

    int is_sinus_rhythm;
    int premature_beat_count;

    int analysis_status;
    char status_message[64];
    char pathologic_status[64];
    char reserved[4];
};

struct pqrst_indxes
{
    int32_t p[3];
    int32_t r[3];
    int32_t t[3];
};

void analyse_ecg_multichannel(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, ecg_analysis_result& result);
ecg_analysis_result analyse_ecg_detect_peaks(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, std::vector<pqrst_indxes>& annotations, std::vector<unsigned int>* peak_indexes, std::string mode = "default");
