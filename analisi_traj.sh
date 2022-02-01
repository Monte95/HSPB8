#!/bin/bash
PATH=$PATH:"~/Downloads/gromacs-2020/build/bin"

#NAME VARIABLES
#traj_name=("it_wt" "it_k" "m_wt" "m_k" "trros_wt" "trros_k")
traj_name=("m_wt" "m_k")
name_specifier="remd_500ns"

for name in ${traj_name[@]}
do

	#RMSD (backbone)
	printf "4\n4\n" | gmx rms -f ${name}.xtc -s ${name}.tpr -n ${name}.ndx -o rmsd_${name}_${name_specifier}.xvg

	#RMSF
	printf "1\n" | gmx rmsf -f ${name}.xtc -s ${name}.tpr -n ${name}.ndx -res -o rmsf_${name}_${name_specifier}.xvg

	#Rg
	printf "1\n" | gmx gyrate -f ${name}.xtc -s ${name}.tpr -n ${name}.ndx -o gyrate_${name}_${name_specifier}.xvg
#	printf "10\n" | gmx gyrate -f ${name}.xtc -s ${name}.tpr -n ${name}_nter.ndx -o gyrate_${name}_${name_specifier}_nter.xvg
	grep "^[^#@]" gyrate_${name}_${name_specifier}.xvg > gyrate_${name}_${name_specifier}_temp.xvg
	mv gyrate_${name}_${name_specifier}_temp.xvg gyrate_${name}_${name_specifier}.xvg

	#CLUSTERING
#	printf "4\n1\n" | gmx cluster -f ${name}.xtc -s ${name}.tpr -n ${name}.ndx -sz size_cluster_${name}_${name_specifier}_020.xvg -cl clusters_structures_${name}_${name_specifier}_020.pdb -cutoff 0.2 -g cluster_${name}_${name_specifier}_020.log

	#SASA + hSASA
#r 1 | r 2 | r 6 | r 7 | r 8 | r 13 | r 16 | r 20 | r 21 | r 25 | r 26 | r 30 | r 31 | r 35 | r 37 | r 39 | r 40 | r 41 | r 44 | r 46 | r 48 | r 49 | r 51 | r 52 | r 53 | r 54 | r 56 | r 59 | r 60 | r 61 | r 64 | r 68 | r 69 | r 70 | r 73 | r 75 | r 77 | r 79 | r 81 | r 82 | r 83 | r 88 | r 89 | r 90 | r 91 | r 92 | r 95 | r 96 | r 98 | r 100 | r 102 | r 105 | r 107 | r 110 | r 111 | r 112 | r 119 | r 121 | r 134 | r 135 | r 139 | r 143 | r 145 | r 146 | r 147 | r 149 | r 151 | r 152 | r 154 | r 155 | r 156 | r 158 | r 160 | r 163 | r 164 | r 165 | r 166 | r 168 | r 169 | r 171 | r 172 | r 173 | r 177 | r 182 | r 186 | r 187 | r 193

#	gmx make_ndx -f ${name}.gro -o hydrophobic_residues_${name}.ndx
	printf "1\n10\n\n" | gmx sasa -f ${name}.xtc -s ${name}.tpr -n hydrophobic_residues_${name}.ndx -o sasa_tot_${name}_${name_specifier}.xvg -surface -output
	
	grep "^[^#@]" sasa_tot_${name}_${name_specifier}.xvg > sasa_tot_${name}_${name_specifier}_temp.xvg
	mv sasa_tot_${name}_${name_specifier}_temp.xvg sasa_tot_${name}_${name_specifier}.xvg
	awk '{print $1, $2}' sasa_tot_${name}_${name_specifier}.xvg > sasa_${name}_${name_specifier}.xvg
	awk '{print $1, $3}' sasa_tot_${name}_${name_specifier}.xvg > sasa_hydrophobic_${name}_${name_specifier}.xvg

done



