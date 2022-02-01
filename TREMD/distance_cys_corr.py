from distance_cys import *

distance_10_99_df = distance_10_99_df.rename(columns={"Distance": "Distance_10_99"}, errors="raise")
distance_10_99_df_ridotto = distance_10_99_df.loc[distance_10_99_df['Homology Model'] == 'All wt', "Distance_10_99"]

distance_10_195_df = distance_10_195_df.rename(columns={"Distance": "Distance_10_195"})
distance_10_195_df_ridotto = distance_10_195_df.loc[distance_10_195_df['Homology Model'] == 'All wt', "Distance_10_195"]

distance_99_195_df = distance_99_195_df.rename(columns={"Distance": "Distance_99_195"})
distance_99_195_df_ridotto = distance_99_195_df.loc[distance_99_195_df['Homology Model'] == 'All wt', "Distance_99_195"]

distance_cys_df = pd.concat([distance_10_99_df_ridotto, distance_10_195_df_ridotto, distance_99_195_df_ridotto], axis=1)
distance_cys_df = distance_cys_df.rename(columns={'Distance_10_99': "Distance 10-99",
                                                  'Distance_10_195': "Distance 10-195",
                                                  'Distance_99_195': "Distance 99-195"})

print(distance_cys_df.corr())


