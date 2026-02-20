#!/bin/bash
# Helper script to run NewCIPLabeler tests

export RDBASE=/Users/dbn/builds/26-2/source/rdkit
export DYLD_LIBRARY_PATH=/Users/dbn/builds/26-2/build/rdkit/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/coordgenlibs-3.0.2/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/boost_1_81_0-py311-BLDMGR-8470/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/zstd-v1.5.5-boost-1_81_0-BLDMGR-7247/lib

TEST_EXEC=/Users/dbn/builds/26-2/build/rdkit/build/Code/GraphMol/NewCIPLabeler/testNewCIPLabeler

if [ ! -f "$TEST_EXEC" ]; then
    echo "Error: Test executable not found at $TEST_EXEC"
    echo "Please build first with: /Users/dbn/builds/26-2/source/mmshare/build_tools/buildinger.sh --name 26-2 rdkit"
    exit 1
fi

# Run with any arguments passed to script
$TEST_EXEC "$@"
