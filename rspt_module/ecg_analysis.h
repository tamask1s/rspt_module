#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct ch_result
{
    int P_DURATION; /** In case of a monophasic P wave this value equals to P1_DURATION. In case of biphasic P it is equal to P1_DURATION + P2_DURATION */

    int P1_DURATION;
    int P1_AMPLITUDE;
    int P2_DURATION;
    int P2_AMPLITUDE;

    int PQ_INTERVAL;

    int Q_DURATION;
    int Q_AMPLITUDE;
    int R_DURATION;
    int R_AMPLITUDE;
    int S_DURATION;
    int S_AMPLITUDE;

    int QRS_DURATION;
    int QT_INTERVAL;

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
    int32_t r[5]; /// r_on_index, r_peak_index, r_off_index, q_idx if > 0, s_idx if > 0
    int32_t t[4]; /// t_on_index, t_peak_index, t_off_index, j_point
};

struct ecg_analysis_config
{
    /* ================= PHYSIOLOGICAL LIMITS ================= */
    struct phys_limits_t
    {
        double max_qrs_duration_ms        = 180.0;   // felső élettani határ a QRS teljes szélességére (outlier / hibás detekció levágása)
        double min_qrs_duration_ms        = 40.0;    // alsó élettani határ a QRS-re (zaj / éles spike kizárása)

        double min_qt_interval_ms         = 250.0;   // minimálisan elfogadott QT intervallum
        double max_qt_interval_ms         = 550.0;   // maximálisan elfogadott QT intervallum

        double max_s_fraction_of_qrs      = 0.5;     // S hullám maximális hossza a teljes QRS-hez viszonyítva (túl hosszú S -> zaj)
        double min_wave_presence_ratio    = 0.01;    // abszolút amplitúdó alsó küszöb: ez alatt a hullám nem tekinthető valósnak
    } phys;

    /* ================= Artifact removal ================= */
    struct artifact_t
    {
        double artifact_removal_ratio = 5.0;         // szomszédos mintákhoz viszonyított extrém eltérés aránya (single-sample spike detektálás)

        double global_std_multiplier = 4.0;          // globális szóráshoz viszonyított küszöb (mean-alapú outlier)
        double local_std_multiplier  = 6.0;          // lokális szomszédos különbségek küszöbe (impulzusszerű zaj)

        bool enable_artifact_removal = true;         // artifact eltávolítás ki/bekapcsolása
    } artifact;

    /* ================= Derivative / L1 norm ================= */
    struct
    {
        double l1norm_dt1_multiplier = 2.0;          // rövid időskála lépésköze (QRS meredekség hangsúlyozása)
        double l1norm_dt2_multiplier = 4.0;          // hosszabb időskála lépésköze (T hullám energia)

        double r_wave_deriv_ms = 6.0;                // R-hullám környezetében használt derivált összehasonlítási időablak
    } deriv;

    /* ================= Isoelectrical search ================= */
    struct
    {
        double qrs_isoel_tolerance = 0.035;          // QRS izoelektromos küszöb a derivált amplitúdójára
        double p_isoel_tolerance   = 0.21;           // P-hullám izoelektromos küszöb (lazább, laposabb hullám)
        double t_isoel_tolerance   = 0.25;           // T-hullám izoelektromos küszöb

        double qrs_dt_ms = 2.0;                      // derivált időablak QRS izoel kereséshez
        double p_dt_ms   = 2.0;                      // derivált időablak P-hez
        double t_dt_ms   = 4.0;                      // derivált időablak T-hez (lassabb lecsengés)

        double qrs_perimeter_ms = 25.0;              // QRS körüli kizárt zóna (csúcs körül nem keresünk izoelt)
        double p_perimeter_ms   = 30.0;              // P-hullám kizárt zóna
        double t_perimeter_ms   = 35.0;              // T-hullám kizárt zóna

        double qrs_max_search_ms = 110.0;             // maximális keresési távolság QRS izoelhez
        double p_max_search_ms   = 200.0;             // maximális keresési távolság P-hez
        double t_max_search_ms   = 200.0;             // maximális keresési távolság T-hez

        int qrs_extend_mode = 0;                     // 0: QRS shrink mód (komplex összehúzása)
        int p_extend_mode   = 1;                     // 1: P esetén jobb oldalról való „lecsúszás”
        int t_extend_mode   = 1;                     // 1: T esetén jobb oldalról való „lecsúszás”
    } isoel;

    /* ================= P, Q, R, S thresholds ================= */
    struct
    {
        double wave_presence_ratio  = 0.01;          // Q és S amplitúdó minimális aránya az R-hez képest
        double wave_bound_ratio     = 0.01;          // hullámhatár keresés küszöbe (find_wave_bounds)
        double r_wave_bound_ratio   = 0.05;          // speciális R-hullám határarány (no-R eset felismerés)

        double qrs_shrink_trshld    = 0.05;          // QRS shrink során megengedett amplitúdó eltérés
        double qrs_shrink_bound_ms  = 10.0;          // QRS shrink maximális elmozdulása időben

        double isoel_r_search_ms    = 30.0;          // T-hullám baseline számításakor használt R utáni offset
    } wave;

    /* ================= S-wave validation ================= */
    struct
    {
        double max_s_duration_sec = 0.05;            // S-hullám maximális hossza (másodpercben)
        double rising_edge_ratio  = 0.01;            // jobb oldali visszaemelkedés minimális aránya (zaj vs. valódi S)
    } s_wave;

    /* ================= Timing offsets ================= */
    struct
    {
        double pq_interval_offset_ms = 4.0;          // PQ intervallum korrekció (izoel keresési bias)
        double qt_interval_offset_ms = 8.0;          // QT intervallum korrekció
    } timing;

    /* ================= Bandpass filter ================= */
    struct
    {
        int nr_filter_iterations = 2;                 // másordendű oda-vissza filter ennyiszer szűri meg a jelet
        double low_cut_hz  = 0.1;                     // alsó vágási frekvencia (baseline wander eltávolítás)
        double high_cut_hz = 40.0;                    // felső vágási frekvencia (EMG / zaj csökkentése)
    } filter;

    /* ================= Search windows / morphology ================= */
    struct
    {
        double p_isoel_offset_ms     = 15.0;          // P-hullám izoel keresésének R-től vett eltolása
        double p_search_left_ms      = 250.0;         // maximális balra keresési ablak P-csúcsra
        double t_search_right_ms     = 300.0;         // maximális jobbra keresési ablak T-csúcsra

        double biphasic_p_max_gap_ms = 20.0;          // biphasic P két komponense közti max. távolság
        double biphasic_p_ratio      = 0.3;           // kisebb/nagyobb amplitúdó arány a biphasic P elfogadásához
    } search;

    /* ================= DATASET BIAS (OPTIONAL) ================= */
    struct dataset_bias_t
    {
        double amplitude_scale            = 1.0;
        double duration_scale             = 1.0;

        double noise_sensitivity_scale    = 1.0;
    } bias;
};

void analyse_ecg_multichannel(const double** ecg_signal, unsigned int nr_ch, unsigned int nr_samples_per_ch, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, ecg_analysis_result& result);
ecg_analysis_result analyse_ecg_detect_peaks(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, std::vector<pqrst_indxes>& annotations, std::vector<unsigned int>* peak_indexes = nullptr, std::string mode = "default", int analysis_ch_indx = -1, int analysis_peak_indx = 0);
void analyse_ecg_all_beats(const double** data, size_t nr_channels, size_t nr_samples_per_channel, double sampling_rate, const std::vector<unsigned int>& peak_indexes, std::vector<pqrst_indxes>& annotations, std::vector<ecg_analysis_result>& results, int analysis_ch_indx = -1);
