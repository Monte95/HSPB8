import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ref_structure_wt = "it_wt.gro"
trajectory_wt = "it_wt.xtc"

ref_structure_k = "it_k.gro"
trajectory_k = "all_k.xtc"

t_wt = md.load(trajectory_wt[0:20], top=ref_structure_wt)
# t_k = md.load(trajectory_k, top=ref_structure_k)

trp = [47, 50, 59, 95]
residues = [i for i in range(0, 196)]




