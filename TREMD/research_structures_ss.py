from secstr_dssp import ref_structure_wt, ref_structure_k, t_wt, t_k, NUM_RESIDUES
from remd_distributions import read_gyration_data, read_sasa_data
import mdtraj as md
import numpy as np
import pandas as pd

NUM_RESIDUES = 196

# ALPHAFOLD HM
alphafold = md.load("alphafold.pdb")
dssp_alphafold = np.transpose(md.compute_dssp(alphafold, simplified=False))
dssp_alphafold_df = pd.DataFrame(dssp_alphafold)
dssp_alphafold_df = dssp_alphafold_df.replace(['H'], 'Helix')
dssp_alphafold_df = dssp_alphafold_df.replace(['G', 'I'], 'Distorted Helix')
dssp_alphafold_df = dssp_alphafold_df.replace(['E'], 'Strand')
dssp_alphafold_df = dssp_alphafold_df.replace(['B'], 'Distorted Strand')
dssp_alphafold_df = dssp_alphafold_df.replace(['T'], 'Turn')
dssp_alphafold_df = dssp_alphafold_df.replace([' ', 'S'], 'Unordered')
print('ALPHAFOLD: ', (dssp_alphafold_df.apply(pd.value_counts).fillna(0)/NUM_RESIDUES).T)

# I-TASSER HM
it = md.load("b8_itasser.pdb")
dssp_it = np.transpose(md.compute_dssp(it, simplified=False))
dssp_it_df = pd.DataFrame(dssp_it)
dssp_it_df = dssp_it_df.replace(['H'], 'Helix')
dssp_it_df = dssp_it_df.replace(['G', 'I'], 'Distorted Helix')
dssp_it_df = dssp_it_df.replace(['E'], 'Strand')
dssp_it_df = dssp_it_df.replace(['B'], 'Distorted Strand')
dssp_it_df = dssp_it_df.replace(['T'], 'Turn')
dssp_it_df = dssp_it_df.replace([' ', 'S'], 'Unordered')
print('I-TASSER: ', (dssp_it_df.apply(pd.value_counts).fillna(0)/NUM_RESIDUES).T)

# MODELLER HM
m = md.load("b8_modeller.pdb")
dssp_m = np.transpose(md.compute_dssp(m, simplified=False))
dssp_m_df = pd.DataFrame(dssp_m)
dssp_m_df = dssp_m_df.replace(['H'], 'Helix')
dssp_m_df = dssp_m_df.replace(['G', 'I'], 'Distorted Helix')
dssp_m_df = dssp_m_df.replace(['E'], 'Strand')
dssp_m_df = dssp_m_df.replace(['B'], 'Distorted Strand')
dssp_m_df = dssp_m_df.replace(['T'], 'Turn')
dssp_m_df = dssp_m_df.replace([' ', 'S'], 'Unordered')
print('MODELLER: ', (dssp_m_df.apply(pd.value_counts).fillna(0)/NUM_RESIDUES).T)

# ROSETTA HM
ros = md.load("b8_trrosetta_72533_model1.pdb")
dssp_ros = np.transpose(md.compute_dssp(ros, simplified=False))
dssp_ros_df = pd.DataFrame(dssp_ros)
dssp_ros_df = dssp_ros_df.replace(['H'], 'Helix')
dssp_ros_df = dssp_ros_df.replace(['G', 'I'], 'Distorted Helix')
dssp_ros_df = dssp_ros_df.replace(['E'], 'Strand')
dssp_ros_df = dssp_ros_df.replace(['B'], 'Distorted Strand')
dssp_ros_df = dssp_ros_df.replace(['T'], 'Turn')
dssp_ros_df = dssp_ros_df.replace([' ', 'S'], 'Unordered')
print('ROSETTA: ', (dssp_ros_df.apply(pd.value_counts).fillna(0)/NUM_RESIDUES).T)

dssp_wt = np.transpose(md.compute_dssp(t_wt, simplified=False))
dssp_k = np.transpose(md.compute_dssp(t_k, simplified=False))

dssp_wt_df = pd.DataFrame(dssp_wt)
dssp_k_df = pd.DataFrame(dssp_k)

dssp_wt_df = dssp_wt_df.replace(['H'], 'Helix')
dssp_wt_df = dssp_wt_df.replace(['G', 'I'], 'Distorted Helix')
dssp_wt_df = dssp_wt_df.replace(['E'], 'Strand')
dssp_wt_df = dssp_wt_df.replace(['B'], 'Distorted Strand')
dssp_wt_df = dssp_wt_df.replace(['T'], 'Turn')
dssp_wt_df = dssp_wt_df.replace([' ', 'S'], 'Unordered')

dssp_k_df = dssp_k_df.replace(['H'], 'Helix')
dssp_k_df = dssp_k_df.replace(['G', 'I'], 'Distorted Helix')
dssp_k_df = dssp_k_df.replace(['E'], 'Strand')
dssp_k_df = dssp_k_df.replace(['B'], 'Distorted Strand')
dssp_k_df = dssp_k_df.replace(['T'], 'Turn')
dssp_k_df = dssp_k_df.replace([' ', 'S'], 'Unordered')

dssp_wt_df_summary = (dssp_wt_df.apply(pd.value_counts).fillna(0)/NUM_RESIDUES).T
dssp_k_df_summary = (dssp_k_df.apply(pd.value_counts).fillna(0)/NUM_RESIDUES).T

str_wt_cd = dssp_wt_df_summary.loc[
          (dssp_wt_df_summary['Helix'] < 0.02) &
          (dssp_wt_df_summary['Distorted Helix'] < 0.06) &
          (dssp_wt_df_summary['Strand'] < 0.3) & (dssp_wt_df_summary['Strand'] > 0.2) &
          (dssp_wt_df_summary['Distorted Strand'] > 0.02) &
          (dssp_wt_df_summary['Turn'] < 0.25) & (dssp_wt_df_summary['Turn'] > 0.15) &
          (dssp_wt_df_summary['Unordered'] < 0.5) & (dssp_wt_df_summary['Unordered'] > 0.3)]

str_k_cd = dssp_k_df_summary.loc[
          (dssp_k_df_summary['Helix'] < 0.02) &
          (dssp_k_df_summary['Distorted Helix'] < 0.06) &
          (dssp_k_df_summary['Strand'] < 0.3) & (dssp_k_df_summary['Strand'] > 0.2) &
          (dssp_k_df_summary['Distorted Strand'] > 0.02) &
          (dssp_k_df_summary['Turn'] < 0.25) & (dssp_k_df_summary['Turn'] > 0.15) &
          (dssp_k_df_summary['Unordered'] < 0.5) & (dssp_k_df_summary['Unordered'] > 0.3)]

t, rg_wt = read_gyration_data('gyrate_all_wt.xvg')
t, rg_k = read_gyration_data('gyrate_all_k.xvg')

t, sasa_wt, hsasa_wt = read_sasa_data('sasa_tot_all_wt.xvg')
t, sasa_k, hsasa_k = read_sasa_data('sasa_tot_all_k.xvg')

print('Rg wt: {}'.format(rg_wt[478]))
print('Rg k: {}, {}'.format(rg_k[3016], rg_k[3401]))

print('hSASA/SASA wt: {}'.format(hsasa_wt[478]/sasa_wt[478]))
print('hSASA/SASA k: {}, {}'.format(hsasa_k[3016]/sasa_k[3016], hsasa_k[3401]/sasa_k[3401]))


