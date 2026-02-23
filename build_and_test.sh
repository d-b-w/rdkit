#!/bin/bash
# Build and test script for NewCIPLabeler development

set -e

RDKIT_ROOT=/Users/dbn/builds/26-2/source/rdkit

echo "=== Building with buildinger ==="
cd "$RDKIT_ROOT"
../mmshare/build_tools/buildinger.sh --name 26-2 rdkit

echo ""
echo "=== Running NewCIPLabeler tests ==="
cd /Users/dbn/builds/26-2/build/rdkit/build
env RDBASE=/Users/dbn/builds/26-2/source/rdkit DYLD_LIBRARY_PATH=/Users/dbn/builds/26-2/build/rdkit/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/coordgenlibs-3.0.2/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/boost_1_81_0-py311-BLDMGR-8470/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/zstd-v1.5.5-boost-1_81_0-BLDMGR-7247/lib Code/GraphMol/NewCIPLabeler/testNewCIPLabeler
