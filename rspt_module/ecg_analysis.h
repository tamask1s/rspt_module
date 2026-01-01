
struct ch_result
{
    int P1_DURATION;
    int P1_AMPLITUDE;
    int P2_DURATION;
    int P2_AMPLITUDE;

    int Q_DURATION;
    int Q_AMPLITUDE;
    int R_DURATION;
    int R_AMPLITUDE;
    int S_DURATION;
    int S_AMPLITUDE;

    int QRS_DURATION;

    int J_AMPLITUDE;
    int ST_20_AMPLITUDE;
    int ST_40_AMPLITUDE;
    int ST_60_AMPLITUDE;
    int ST_80_AMPLITUDE;
    int T_AMPLITUDE;
};

struct ecg_analysis_result
{
    struct ch_result result;
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
    int32_t p[6]; /// p1_on_index, p1_peak_index, p1_off_index, p2_on_index, p2_peak_index, p2_off_index. If p2_on_index, p2_peak_index, p2_off_index are set to 0, p2 is not present.
    int32_t r[3]; /// r_on_index, r_peak_index, r_off_index
    int32_t t[3]; /// t_on_index, t_peak_index, t_off_index
};

void analyse_ecg_multichannel(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, ecg_analysis_result& result);
ecg_analysis_result analyse_ecg_detect_peaks(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, std::vector<pqrst_indxes>& annotations, std::vector<unsigned int>* peak_indexes = nullptr, std::string mode = "default", int analysis_ch_indx = -1);
