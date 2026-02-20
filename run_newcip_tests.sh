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

# Run with timeout to prevent infinite loops (30 second default)
# Usage: TIMEOUT=60 ./run_newcip_tests.sh "[phase2]"
TIMEOUT="${TIMEOUT:-30}"

# macOS doesn't have timeout, try gtimeout or use perl workaround
if command -v gtimeout >/dev/null 2>&1; then
    gtimeout "${TIMEOUT}s" $TEST_EXEC "$@"
    EXIT_CODE=$?
elif command -v timeout >/dev/null 2>&1; then
    timeout "${TIMEOUT}s" $TEST_EXEC "$@"
    EXIT_CODE=$?
else
    # Fallback: run with perl-based timeout
    perl -e "alarm ${TIMEOUT}; exec @ARGV" $TEST_EXEC "$@"
    EXIT_CODE=$?
    if [ $EXIT_CODE -eq 142 ]; then  # SIGALRM
        echo ""
        echo "ERROR: Tests timed out after ${TIMEOUT} seconds - likely infinite loop!"
        exit 124
    fi
fi

if [ $EXIT_CODE -eq 124 ]; then
    echo ""
    echo "ERROR: Tests timed out after ${TIMEOUT} seconds - likely infinite loop!"
    exit 124
fi

exit $EXIT_CODE
