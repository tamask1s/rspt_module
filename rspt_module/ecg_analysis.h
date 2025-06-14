
struct ecg_analysis_result
{
    double rr_interval_ms;
    double rr_variation_ms;
    double heart_rate_bpm;
    double pr_interval_ms;
    double qrs_duration_ms;
    double qt_interval_ms;
    double qtc_interval_ms;

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
};

void analyse_ecg_multichannel(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<std::string>& annotations, ecg_analysis_result& result);
