#!/bin/bash
set -xe

# We already have the generated topology
# 100k.gro is equivalent to filamento-chitosa_processed.gro (I guess) 

echo "Simulating a simple aqueous system"
gmx editconf -f 100k.gro -o filamento-chitosa_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp filamento-chitosa_newbox.gro -cs spc216.gro -o filamento-chitosa_solv.gro -p topol.top
echo "------------------------------------------------------------------------"

echo "Add ions..."
if [ ! -e "ions.mdp" ]; then
    curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/ions.mdp
fi
gmx grompp -f ions.mdp -c filamento-chitosa_solv.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o filamento-chitosa_solv_ions.gro -p topol.top -pname NA -nname CL -neutral <<EOF
13
EOF
echo "------------------------------------------------------------------------"

echo "Energy Minimization..."
if [ ! -e "minim.mdp" ]; then
    curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/minim.mdp
fi
gmx grompp -f minim.mdp -c filamento-chitosa_solv_ions.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em
echo "------------------------------------------------------------------------"

echo "Equilibration - Phase 1..."
if [ ! -e "nvt.mdp" ]; then
    curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/nvt.mdp
fi
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt
echo "------------------------------------------------------------------------"

echo "Equilibration - Phase 2..."
if [ ! -e "npt.mdp" ]; then
    curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/npt.mdp
fi
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt
echo "------------------------------------------------------------------------"

echo "Production MD"
if [ ! -e "md.mdp" ]; then
    curl -O http://www.mdtutorials.com/gmx/lysozyme/Files/md.mdp
fi
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr
gmx mdrun -v -deffnm md_0_1
echo "------------------------------------------------------------------------"

### Running GROMACS ON GPU
# gmx mdrun -deffnm md_0_1 -nb gp

## Analysis
# gmx trjconv -s md_0_1.tpr -f md_0_1.xtc -o md_0_1_noPBC.xtc -pbc mol -center
# gmx rms -s md_0_1.tpr -f md_0_1_noPBC.xtc -o rmsd.xvg -tu ns
# gmx rms -s em.tpr -f md_0_1_noPBC.xtc -o rmsd_xtal.xvg -tu ns
