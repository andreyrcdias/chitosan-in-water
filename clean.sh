#!/bin/bash
set -x

echo "Cleaning target files..."
find . -maxdepth 1 ! -name 'clean.sh' ! -name 'run.sh' ! -name 'filamento-chitosa.pdb' ! -name 100k.gro ! -name 'README.md' -exec rm -f {} +