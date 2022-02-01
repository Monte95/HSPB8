import pandas as pd
import numpy as np


def read_sec_str_data(document_name):
    y = []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 11:
                y.append(cols[6])
    return np.array(y)


secstr_wt_7 = read_sec_str_data("stride_generated_encodermap_wt_7.txt")
secstr_wt_8 = read_sec_str_data("stride_generated_encodermap_wt_8.txt")
secstr_wt_9 = read_sec_str_data("stride_generated_encodermap_wt_9.txt")
secstr_wt_10 = read_sec_str_data("stride_generated_encodermap_wt_10.txt")
secstr_wt_11 = read_sec_str_data("stride_generated_encodermap_wt_11.txt")
secstr_wt_12 = read_sec_str_data("stride_generated_encodermap_wt_12.txt")
secstr_wt_13 = read_sec_str_data("stride_generated_encodermap_wt_13.txt")
secstr_wt_14 = read_sec_str_data("stride_generated_encodermap_wt_14.txt")
secstr_wt_15 = read_sec_str_data("stride_generated_encodermap_wt_15.txt")
secstr_wt_16 = read_sec_str_data("stride_generated_encodermap_wt_16.txt")
secstr_wt_17 = read_sec_str_data("stride_generated_encodermap_wt_17.txt")
secstr_wt_18 = read_sec_str_data("stride_generated_encodermap_wt_18.txt")
secstr_wt_19 = read_sec_str_data("stride_generated_encodermap_wt_19.txt")

secstr_k_1 = read_sec_str_data("stride_generated_encodermap_k_1.txt")
secstr_k_2 = read_sec_str_data("stride_generated_encodermap_k_2.txt")
secstr_k_3 = read_sec_str_data("stride_generated_encodermap_k_3.txt")
secstr_k_4 = read_sec_str_data("stride_generated_encodermap_k_4.txt")
secstr_k_5 = read_sec_str_data("stride_generated_encodermap_k_5.txt")
secstr_k_6 = read_sec_str_data("stride_generated_encodermap_k_6.txt")
secstr_k_20 = read_sec_str_data("stride_generated_encodermap_k_20.txt")
secstr_k_21 = read_sec_str_data("stride_generated_encodermap_k_21.txt")
secstr_k_22 = read_sec_str_data("stride_generated_encodermap_k_22.txt")

secstr_wt_7_df = pd.DataFrame({"Secondary structure": secstr_wt_7})
secstr_wt_8_df = pd.DataFrame({"Secondary structure": secstr_wt_8})
secstr_wt_9_df = pd.DataFrame({"Secondary structure": secstr_wt_9})
secstr_wt_10_df = pd.DataFrame({"Secondary structure": secstr_wt_10})
secstr_wt_11_df = pd.DataFrame({"Secondary structure": secstr_wt_11})
secstr_wt_12_df = pd.DataFrame({"Secondary structure": secstr_wt_12})
secstr_wt_13_df = pd.DataFrame({"Secondary structure": secstr_wt_13})
secstr_wt_14_df = pd.DataFrame({"Secondary structure": secstr_wt_14})
secstr_wt_15_df = pd.DataFrame({"Secondary structure": secstr_wt_15})
secstr_wt_16_df = pd.DataFrame({"Secondary structure": secstr_wt_16})
secstr_wt_17_df = pd.DataFrame({"Secondary structure": secstr_wt_17})
secstr_wt_18_df = pd.DataFrame({"Secondary structure": secstr_wt_18})
secstr_wt_19_df = pd.DataFrame({"Secondary structure": secstr_wt_19})

secstr_k_1_df = pd.DataFrame({"Secondary structure": secstr_k_1})
secstr_k_2_df = pd.DataFrame({"Secondary structure": secstr_k_2})
secstr_k_3_df = pd.DataFrame({"Secondary structure": secstr_k_3})
secstr_k_4_df = pd.DataFrame({"Secondary structure": secstr_k_4})
secstr_k_5_df = pd.DataFrame({"Secondary structure": secstr_k_5})
secstr_k_6_df = pd.DataFrame({"Secondary structure": secstr_k_6})
secstr_k_20_df = pd.DataFrame({"Secondary structure": secstr_k_20})
secstr_k_21_df = pd.DataFrame({"Secondary structure": secstr_k_21})
secstr_k_22_df = pd.DataFrame({"Secondary structure": secstr_k_22})

secstr_wt_df = pd.concat([secstr_wt_7_df, secstr_wt_8_df, secstr_wt_9_df, secstr_wt_10_df, secstr_wt_11_df,
                          secstr_wt_12_df, secstr_wt_13_df, secstr_wt_14_df, secstr_wt_15_df, secstr_wt_16_df,
                          secstr_wt_17_df, secstr_wt_18_df,secstr_wt_19_df], ignore_index=True)

secstr_k_df = pd.concat([secstr_k_1_df, secstr_k_2_df, secstr_k_3_df, secstr_k_4_df, secstr_k_5_df, secstr_k_6_df,
                         secstr_k_20_df, secstr_k_21_df, secstr_k_22_df], ignore_index=True)

# Coil = 0
# Turn = 1
# AlphaHelix = 2
# 310Helix = 3
# Strand = 4
# Bridge = 5

secstr_wt_df.loc[secstr_wt_df['Secondary structure'].str.contains("Coil"), 'SS'] = 0
secstr_wt_df.loc[secstr_wt_df['Secondary structure'].str.contains("Turn"), 'SS'] = 1
secstr_wt_df.loc[secstr_wt_df['Secondary structure'].str.contains("AlphaHelix"), 'SS'] = 2
secstr_wt_df.loc[secstr_wt_df['Secondary structure'].str.contains("310Helix"), 'SS'] = 3
secstr_wt_df.loc[secstr_wt_df['Secondary structure'].str.contains("Strand"), 'SS'] = 4
secstr_wt_df.loc[secstr_wt_df['Secondary structure'].str.contains("Bridge"), 'SS'] = 5

secstr_k_df.loc[secstr_k_df['Secondary structure'].str.contains("Coil"), 'SS'] = 0
secstr_k_df.loc[secstr_k_df['Secondary structure'].str.contains("Turn"), 'SS'] = 1
secstr_k_df.loc[secstr_k_df['Secondary structure'].str.contains("AlphaHelix"), 'SS'] = 2
secstr_k_df.loc[secstr_k_df['Secondary structure'].str.contains("310Helix"), 'SS'] = 3
secstr_k_df.loc[secstr_k_df['Secondary structure'].str.contains("Strand"), 'SS'] = 4
secstr_k_df.loc[secstr_k_df['Secondary structure'].str.contains("Bridge"), 'SS'] = 5

# wt secondary structure percentages
secstr_wt_9_df.loc[secstr_wt_9_df['Secondary structure'].str.contains("Coil"), 'SS'] = 0
secstr_wt_9_df.loc[secstr_wt_9_df['Secondary structure'].str.contains("Turn"), 'SS'] = 1
secstr_wt_9_df.loc[secstr_wt_9_df['Secondary structure'].str.contains("AlphaHelix"), 'SS'] = 2
secstr_wt_9_df.loc[secstr_wt_9_df['Secondary structure'].str.contains("310Helix"), 'SS'] = 3
secstr_wt_9_df.loc[secstr_wt_9_df['Secondary structure'].str.contains("Strand"), 'SS'] = 4
secstr_wt_9_df.loc[secstr_wt_9_df['Secondary structure'].str.contains("Bridge"), 'SS'] = 5

print(secstr_wt_9_df.groupby('Secondary structure').count() / len(secstr_wt_9_df['Secondary structure']) * 100)

# k secondary structure percentages
#print(secstr_k_df.groupby('Secondary structure').count() / len(secstr_k_df['Secondary structure']) * 100)

#Helix, Sheet, Turn, Others