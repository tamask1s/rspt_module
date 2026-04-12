import pandas as pd

# Fájlnevek
file_a = "ecg_results_rspt.csv"
file_b = "ecg_parsed.csv"
output_file = "ecg_diff_test.csv"

# Beolvasás
df_a = pd.read_csv(file_a)
df_b = pd.read_csv(file_b)

key_cols = ["block", "lead"]

df_a.set_index(key_cols, inplace=True)
df_b.set_index(key_cols, inplace=True)

common_index = df_a.index.intersection(df_b.index)
df_a = df_a.loc[common_index]
df_b = df_b.loc[common_index]

diff_df = df_a - df_b
diff_df.reset_index(inplace=True)

# Numerikus oszlopok
numeric_cols = diff_df.select_dtypes(include="number").columns

# ABS AVG, SD, SIGNED AVG, OUTLIERS, AMP FAILURES, PASS/FAIL számítás
# CSE szabvány szerint: 
# - Duration mérések: 4 szélsőérték eltávolítása
# - Amplitude mérések: 2 szélsőérték eltávolítása
def remove_outliers(series, count=2):
    """Remove extreme values as per CSE standard"""
    if len(series) <= count + 2:  # Ha túl kevés adat
        return series
    series_clean = series.copy()
    for _ in range(count // 2):
        # Eltávolítjuk a legkisebb és legnagyobb értéket
        series_clean = series_clean.drop(series_clean.idxmin())
        series_clean = series_clean.drop(series_clean.idxmax())
    return series_clean

# Statisztikák számítása outlier eltávolítással
abs_avg = {}
std_dev = {}
signed_avg = {}
outliers_count = {}
clean_data_dict = {}  # Tisztított adatok tárolása

for col in numeric_cols:
    data = diff_df[col].dropna()
    if len(data) > 0:
        # Outlier eltávolítás CSE szabvány szerint
        if 'DUR' in col:
            # Duration: 4 outlier eltávolítás
            clean_data = remove_outliers(data, 4)
        elif 'AMP' in col:
            # Amplitude: 2 outlier eltávolítás
            clean_data = remove_outliers(data, 2)
        else:
            clean_data = data
            
        clean_data_dict[col] = clean_data  # Tárolás későbbi használatra
        abs_avg[col] = clean_data.abs().mean()
        std_dev[col] = clean_data.std()
        signed_avg[col] = clean_data.mean()
        
        # Outliers számítás a tisztított adatokon
        outliers_count[col] = ((clean_data < -10) | (clean_data > 10)).sum()
    else:
        clean_data_dict[col] = data
        abs_avg[col] = 0
        std_dev[col] = 0
        signed_avg[col] = 0
        outliers_count[col] = 0

# Amplitúdó szabvány megsértések számolása az EREDETI adatokon
# (Az outlier eltávolítás csak a statisztikákra vonatkozik, nem a failures számításra)
amp_failures = {}
for col in numeric_cols:
    if 'AMP' in col:
        failures = 0
        for idx in diff_df.index[:-8]:  # Csak adatsorok, statisztikai sorok nélkül
            diff_val = diff_df.loc[idx, col]
            if pd.notna(diff_val) and diff_val != "":
                diff_val = float(diff_val)
                # Referencia érték keresése
                block = diff_df.loc[idx, 'block']
                lead = diff_df.loc[idx, 'lead']
                
                # Referencia érték keresése az indexelt df_b-ből
                try:
                    ref_val = abs(df_b.loc[(block, lead), col])
                    
                    # CSE szabvány alkalmazása
                    if ref_val <= 500:
                        limit = 25  # ±25 µV
                    else:
                        limit = ref_val * 0.1  # ±10%
                    
                    if abs(diff_val) > limit:
                        failures += 1
                except (KeyError, IndexError):
                    # Ha nincs referencia érték, kihagyjuk
                    pass
        amp_failures[col] = failures
    else:
        amp_failures[col] = " "  # DUR oszlopokhoz üres space

# Outliers dictionary létrehozása
outliers = outliers_count

# Szabvány szerinti elfogadható értékek
acceptable_mean = {
    'P1_DUR': 10, 'P2_DUR': 10, 'QRS_DUR': 8, 'Q_DUR': 6, 'R_DUR': 6, 'S_DUR': 6
}
acceptable_sd = {
    'P1_DUR': 8, 'P2_DUR': 8, 'QRS_DUR': 5, 'Q_DUR': 5, 'R_DUR': 5, 'S_DUR': 5
}

# PASS/FAIL értékelés
pass_fail_mean = {}
pass_fail_sd = {}
for col in numeric_cols:
    if col in acceptable_mean:
        # Időtartam mérések
        pass_fail_mean[col] = "PASS" if abs(signed_avg[col]) <= acceptable_mean[col] else "FAIL"
        pass_fail_sd[col] = "PASS" if std_dev[col] <= acceptable_sd[col] else "FAIL"
    elif 'AMP' in col:
        # Amplitúdó mérések - pontos szabály implementálása
        # Referencia értékek becslése a mért értékekből (signed_avg + átlagos referencia)
        # Egyszerűsítés: ha |signed_avg| > 500, akkor nagy amplitúdó
        if abs(signed_avg[col]) > 500:
            # Nagy amplitúdó: ±10% szabály
            limit = abs(signed_avg[col]) * 0.1
        else:
            # Kis amplitúdó: ±25 µV szabály
            limit = 25
        pass_fail_mean[col] = "PASS" if abs(signed_avg[col]) <= limit else "FAIL"
        pass_fail_sd[col] = " "  # Szabvány nem definiál szórás limitet amplitúdókra
    else:
        pass_fail_mean[col] = " "
        pass_fail_sd[col] = " "

abs_avg_row = {col: abs_avg[col] for col in numeric_cols}
abs_avg_row["block"] = "ABS AVG"
abs_avg_row["lead"] = ""

std_row = {col: std_dev[col] for col in numeric_cols}
std_row["block"] = "SD"
std_row["lead"] = ""

signed_avg_row = {col: signed_avg[col] for col in numeric_cols}
signed_avg_row["block"] = "SIGNED AVG"
signed_avg_row["lead"] = ""

outliers_row = {col: outliers[col] for col in numeric_cols}
outliers_row["block"] = "OUTLIERS >10"
outliers_row["lead"] = ""

amp_failures_row = {col: amp_failures[col] for col in numeric_cols}
amp_failures_row["block"] = "AMP FAIL >25uV"
amp_failures_row["lead"] = ""

pass_fail_mean_row = {col: pass_fail_mean[col] for col in numeric_cols}
pass_fail_mean_row["block"] = "PASS/FAIL MEAN"
pass_fail_mean_row["lead"] = ""

pass_fail_amp_row = {}
for col in numeric_cols:
    if 'AMP' in col:
        # CSE szabvány szerint ≤2 failure még elfogadható amplitúdó méréseknél
        pass_fail_amp_row[col] = "PASS" if amp_failures[col] <= 2 else "FAIL"
    else:
        pass_fail_amp_row[col] = " "
pass_fail_amp_row["block"] = "PASS/FAIL AMP"
pass_fail_amp_row["lead"] = ""

pass_fail_sd_row = {col: pass_fail_sd[col] for col in numeric_cols}
pass_fail_sd_row["block"] = "PASS/FAIL SD"
pass_fail_sd_row["lead"] = ""

# Statisztikai sorok hozzáadása
diff_df = pd.concat(
    [diff_df, pd.DataFrame([abs_avg_row, std_row, signed_avg_row, outliers_row, amp_failures_row, pass_fail_amp_row, pass_fail_sd_row, pass_fail_mean_row])],
    ignore_index=True
)

# ---- FORMÁZÁS CSV-HEZ ----
stat_rows = diff_df.index[-8:-5]
outliers_row_idx = diff_df.index[-5]
amp_failures_row_idx = diff_df.index[-4]
pass_fail_amp_row_idx = diff_df.index[-3]
pass_fail_rows = diff_df.index[-2:]

for col in numeric_cols:
    # normál sorok: egész szám
    diff_df.loc[diff_df.index[:-8], col] = (
        diff_df.loc[diff_df.index[:-8], col].astype(int).astype(str)
    )
    # statisztikai sorok: 2 tizedes
    for idx in stat_rows:
        diff_df.loc[idx, col] = f"{diff_df.loc[idx, col]:.2f}"
    # outliers és amp failures sorok: egész szám vagy space
    diff_df.loc[outliers_row_idx, col] = str(int(diff_df.loc[outliers_row_idx, col]))
    if diff_df.loc[amp_failures_row_idx, col] == " ":
        diff_df.loc[amp_failures_row_idx, col] = " "
    else:
        diff_df.loc[amp_failures_row_idx, col] = str(int(diff_df.loc[amp_failures_row_idx, col]))
    # pass/fail sorok: szöveg (már megfelelő formátumban vannak)

# Mentés
diff_df.to_csv(output_file, index=False)

print("\nStatisztikák:")
print(diff_df.iloc[-8:].to_string(index=False))

print(f"Kész: {output_file}")
