#!/bin/bash
echo "Cleaning target files..."

exclude_files=(
    'Makefile' 
    'README.md'
    'run.sh'
    'run_v2.sh'
    'run_v3.sh'
    'clean.sh'
    'LCQC_oplsaa.ff.zip'
    'filamento-chitosa.pdb'
    '100k.gro'
    'ions.mdp'
    'minim.mdp'
    'nvt.mdp'
    'npt.mdp'
    'md.mdp'
)

exclude_pattern=$(printf "! -name '%s' " "${exclude_files[@]}")
eval "find . -maxdepth 1 $exclude_pattern -exec rm -f {} +"
