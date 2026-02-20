#!/bin/bash
# Helper script to run NewCIPLabeler tests with proper environment and timeout

BUILD_DIR=/Users/dbn/builds/26-2/build/rdkit/build
TEST_EXEC=Code/GraphMol/NewCIPLabeler/testNewCIPLabeler

# Change to build directory first
cd "$BUILD_DIR" || exit 1

if [ ! -f "$TEST_EXEC" ]; then
    echo "Error: Test executable not found at $BUILD_DIR/$TEST_EXEC"
    echo "Please build first with: /Users/dbn/builds/26-2/source/mmshare/build_tools/buildinger.sh --name 26-2 rdkit"
    exit 1
fi

# Run directly with environment variables (matching Makefile exactly)
# Usage: TIMEOUT=60 ./run_newcip_tests.sh "[phase2]"
TIMEOUT="${TIMEOUT:-30}"

# Create timeout wrapper function
run_with_timeout() {
    if command -v gtimeout >/dev/null 2>&1; then
        gtimeout "${TIMEOUT}s" "$@"
    elif command -v timeout >/dev/null 2>&1; then
        timeout "${TIMEOUT}s" "$@"
    else
        # Fallback: run with perl-based timeout
        perl -e "alarm ${TIMEOUT}; exec @ARGV" "$@"
        local EXIT_CODE=$?
        if [ $EXIT_CODE -eq 142 ]; then  # SIGALRM
            return 124
        fi
        return $EXIT_CODE
    fi
}

# Run test with timeout and proper environment
run_with_timeout env RDBASE=/Users/dbn/builds/26-2/source/rdkit \
    DYLD_LIBRARY_PATH=/Users/dbn/builds/26-2/build/rdkit/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/coordgenlibs-3.0.2/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/boost_1_81_0-py311-BLDMGR-8470/lib:/Users/dbn/builds/26-2/build/schrodinger_buildenv_packages/.pixi/envs/schrodinger/zstd-v1.5.5-boost-1_81_0-BLDMGR-7247/lib \
    "$TEST_EXEC" "$@"
EXIT_CODE=$?

if [ $EXIT_CODE -eq 124 ]; then
    echo ""
    echo "ERROR: Tests timed out after ${TIMEOUT} seconds - likely infinite loop!"
    exit 124
fi

exit $EXIT_CODE
