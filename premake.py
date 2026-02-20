#!/usr/bin/env python3
"""
RDKit premake script for buildinger

Move this into rdkit/premake and make it executable.
Tested only on Mac. It will need light changes for Linux.
"""

import argparse
import glob
import os
import platform
import subprocess
import sys
from pathlib import Path


def run_command(cmd):
    """Run a command and return its output."""
    result = subprocess.run(
        cmd,
        check=True,
        text=True,
        capture_output=True
    )
    return result.stdout.strip()


def setup_build_env(os_cpu, version):
    """Source the build_env and figure out what env variables it sets"""

    # Source the mmshare build environment
    schrodinger_src = os.environ.get("SCHRODINGER_SRC")
    if not schrodinger_src:
        print("Error: SCHRODINGER_SRC environment variable not set.")
        sys.exit(1)

    build_env_cmd = f". {schrodinger_src}/mmshare/build_env -b {os_cpu} && printenv"
    env = dict(os.environ)
    env['OS_CPU'] = os_cpu
    if version:
        env['VERSION'] = version

    # build_env_cmd = 'printenv'
    res = subprocess.run(build_env_cmd, env=env, shell=True,
                         stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE,
                        universal_newlines=True)
    vars = set()
    vars.add('SCHRODINGER_LIB=')
    for line in res.stderr.splitlines():
        if line.startswith('Setting'):
            vars.add(line.split()[1] + '=')

    vars = tuple(vars)
    for line in res.stdout.splitlines():
        if line.startswith(vars):
            k, v = line.split('=', maxsplit=1)
            env[k] = v

    return env


def get_library_definitions(env):
    """Add relevant library definitions to env"""
    schrodinger = os.environ.get("SCHRODINGER")
    schrodinger_src = os.environ.get("SCHRODINGER_SRC")

    if not all([schrodinger, schrodinger_src]):
        print("Error: SCHRODINGER and SCHRODINGER_SRC environment variables must be set.")
        sys.exit(1)

    if 'OS_CPU' not in env:
        raise RuntimeError('The build environment needs to be sourced first!')

    sys.path.insert(0, f"{schrodinger_src}/mmshare/build_tools")
    import library_definitions
    library_definitions.OS_CPU = env['OS_CPU']
    env_res = os.environ
    try:
        os.environ = env
        defs = library_definitions.get_library_definitions()
    finally:
        os.environ = env_res

    env['BOOST_PATH'] = defs['BOOST'].base_path
    env['COORDGEN_DIR'] = defs['COORDGENLIBS'].base_path
    env['EIGEN3_INCLUDE_DIR'] = defs['EIGEN'].include_directory
    env['MAEPARSER_DIR'] = defs['COORDGENLIBS'].base_path
    env['ZSTD_DIR'] = defs['ZSTD'].base_path

    # Per platform config specifics
    if sys.platform == 'linux':
        env['CC'] = 'gcc'
        env['CXX'] = 'g++'
        # On Linux we use the system's installed freetype, so don't set any env var for it
    elif sys.platform == 'darwin':
        env['CC'] = 'clang'
        env['CXX'] = 'clang++'
        env['FREETYPE_DIR'] = defs['FREETYPE'].base_path
    else:
        raise RuntimeError('Unsupported platform!')

    # Set ZLIB_PATH - not yet moved to library definitions
    env["ZLIB_PATH"] = f"{env['SCHRODINGER_LIB']}/{env['OS_CPU']}/zlib-1.2.11"


def setup_python_env(env):
    """Setup Python-related environment variables."""
    schrodinger = os.environ.get('SCHRODINGER')

    # Set Python-related variables
    python_exe = f'{schrodinger}/internal/bin/python3'
    env['PYTHON_EXE'] = python_exe

    python_site_path = run_command([python_exe, '-c', 'import site; print(site.getsitepackages()[0])'])
    env['PYTHON_SITE_PATH'] = python_site_path
    python_include_path = run_command([python_exe, '-c', "import sysconfig; print(sysconfig.get_paths()['include'])"])
    env['PYTHON_INCLUDE_PATH'] = python_include_path
    python_root_dir = str(Path(python_exe).parent.parent)
    env['Python3_ROOT_DIR'] = python_root_dir
    env['NUMPY_INCLUDE_PATH'] = run_command([python_exe, '-c', 'import numpy; print(numpy.get_include())'])
    env['PYTHON_LIB'] = glob.glob(f'{python_root_dir}/lib/libpython*.*')[0]


def run_cmake(cmake_build_type, env):
    """Run CMake to configure the build."""
    schrodinger_src = env['SCHRODINGER_SRC']
    build_dir = env['BUILD_DIR']

    # Check if buildvenv has changed and clear CMake cache if needed
    buildvenv_version = env.get('SCHRODINGER_BUILD_ENV_VERSION', '')
    version_file = os.path.join(build_dir, '.buildvenv_version')
    cache_file = os.path.join(build_dir, 'build', 'CMakeCache.txt')

    previous_version = ''
    if os.path.exists(version_file):
        with open(version_file, 'r') as f:
            previous_version = f.read().strip()

    if buildvenv_version and previous_version != buildvenv_version:
        if os.path.exists(cache_file):
            print(f"Build environment changed ({previous_version} -> {buildvenv_version}), removing CMake cache")
            os.remove(cache_file)

        # Store the new version
        with open(version_file, 'w') as f:
            f.write(buildvenv_version)

    python_exe = env['PYTHON_EXE']
    python_lib = env['PYTHON_LIB']
    python_include_path = env['PYTHON_INCLUDE_PATH']
    numpy_include_path = env['NUMPY_INCLUDE_PATH']
    python_root_dir = env['Python3_ROOT_DIR']
    maeparser_dir = env['MAEPARSER_DIR']
    coordgen_dir = env['COORDGEN_DIR']
    eigen3_include_dir = env['EIGEN3_INCLUDE_DIR']
    boost_path = env['BOOST_PATH']
    freetype_dir = env.get('FREETYPE_DIR', '')

    # Handle CPU architecture specific flags
    cmake_extra_flags = "-DRDK_OPTIMIZE_POPCNT=ON"
    build_x86_64_on_mac = ""
    if platform.machine() == "arm64":
        cmake_extra_flags = "-DRDK_OPTIMIZE_POPCNT=OFF"
        build_x86_64_on_mac = "-DCMAKE_OSX_ARCHITECTURES=x86_64"

    dyn_lib_ext = "dylib" if platform.system() == "Darwin" else "so"

    # Prepare CMake command
    cmake_cmd = [
        "cmake", "-Wno-dev", "-G", "Ninja", f"{schrodinger_src}/rdkit"
    ]

    # Add build architecture if specified
    if build_x86_64_on_mac:
        cmake_cmd.append(build_x86_64_on_mac)


    # Add CMake options
    cmake_cmd.extend([
        f"-DCMAKE_BUILD_TYPE={cmake_build_type}",
        f"-DCMAKE_INSTALL_PREFIX={build_dir}",
        "-DRDK_BUILD_PYTHON_WRAPPERS=ON",
        "-DRDK_INSTALL_INTREE=OFF",
        "-DRDK_INSTALL_STATIC_LIBS=OFF",
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        f"-DMAEPARSER_DIR={maeparser_dir}",
        f"-DCOORDGEN_DIR={coordgen_dir}",
        f"-DEIGEN3_INCLUDE_DIR={eigen3_include_dir}",

        f"-DBoost_ROOT={boost_path}",
        "-DBoost_USE_MULTITHREADED=ON",
        "-DBoost_USE_STATIC_RUNTIME=OFF",
        "-DBoost_NO_BOOST_CMAKE=ON",
        "-DFIND_PACKAGE_MESSAGE_DETAILS_Boost=ON",

        f"-DPYTHON_EXECUTABLE={python_exe}",
        f"-DPython3_EXECUTABLE={python_exe}",
        f"-DPYTHON_LIBRARY={python_lib}",
        f"-DPython3_LIBRARY={python_lib}",
        f"-DPython3_LIBRARIES={python_lib}",
        f"-DPython3_ROOT_DIR={python_root_dir}",
        f"-DPYTHON_INCLUDE_DIR={python_include_path}",
        f"-DPython3_INCLUDE_DIR={python_include_path}",
        f"-DPYTHON_NUMPY_INCLUDE_DIR={numpy_include_path}",
        f"-DPython3_NumPy_INCLUDE_DIRS={numpy_include_path}",

        "-DRDK_INSTALL_PYTHON_TESTS=ON",
        "-DRDK_BUILD_CPP_TESTS=ON",
        "-DRDK_BUILD_TEST_GZIP=ON",

        "-D_HAS_AUTO_PTR_ETC=0",

        "-DPy_ENABLE_SHARED=1",
        "-DRDK_BUILD_COORDGEN_SUPPORT=ON",
        "-DRDK_BUILD_AVALON_SUPPORT=ON",
        "-DRDK_BUILD_COMPRESSED_SUPPLIERS=ON",
        "-DRDK_BUILD_CONTRIB=ON",
        "-DRDK_BUILD_FREESASA_SUPPORT=ON",
        "-DRDK_BUILD_INCHI_SUPPORT=ON",
        "-DRDK_BUILD_QT_SUPPORT=ON",
        "-DRDK_BUILD_RPATH_SUPPORT=ON",
        "-DRDK_BUILD_SLN_SUPPORT=ON",
        "-DRDK_BUILD_THREADSAFE_SSS=ON",
        "-DRDK_BUILD_XYZ2MOL_SUPPORT=ON",
        "-DRDK_BUILD_YAEHMOP_SUPPORT=ON",
        "-DRDK_USE_QT6=ON",
        "-DRDK_USE_URF=ON"
    ])

    # Add Freetype options if available
    if freetype_dir:
        cmake_cmd.extend([
            f"-DFREETYPE_INCLUDE_DIRS={freetype_dir}/include",
            f"-DFREETYPE_LIBRARY={freetype_dir}/lib/libfreetype.{dyn_lib_ext}"
        ])

    # Add extra flags
    cmake_cmd.append(cmake_extra_flags)

    # Change to build directory and run CMake
    os.makedirs(build_dir, exist_ok=True)
    res = subprocess.run(cmake_cmd, cwd=build_dir, env=env)
    if res.returncode:
        print("Error: CMake configuration failed.")
        sys.exit(1)


def create_makefile(os_cpu, version, env):
    """Create a Makefile for building and testing."""
    build_dir = env['BUILD_DIR']
    schrodinger_src = env['SCHRODINGER_SRC']
    coordgen_dir = env['COORDGEN_DIR']
    boost_path = env['BOOST_PATH']
    zstd_dir = env['ZSTD_DIR']

    library_path = f"{build_dir}/lib:{coordgen_dir}/lib:{boost_path}/lib:{zstd_dir}/lib"

    makefile_content = f"""
OS_CPU={os_cpu}
VERSION={version}
all:
	ninja -C build install

test:
	cd build; env RDBASE={schrodinger_src}/rdkit DYLD_LIBRARY_PATH={library_path} ctest
"""

    with open(os.path.join(build_dir, "Makefile"), "w") as f:
        f.write(makefile_content)


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="RDKit premake script for Schrodinger build system")
    parser.add_argument("version", nargs="?", help="Version suffix")
    version = parser.parse_args().version

    OS_CPU_OPTS = {
        'linux': 'Linux-x86_64',
        'darwin': 'Darwin-x86_64',
        'win32': 'Windows-x64',
        'cygwin': 'Windows-x64'
    }

    CMAKE_BUILD_TYPES = {
        '-g': 'Debug',
        '.g': 'RelWithDebInfo',
        '': 'Release'
    }

    schrodinger = os.environ['SCHRODINGER']
    build_dir = f'{schrodinger}/rdkit'
    os_cpu = OS_CPU_OPTS[sys.platform]
    if version.lower().startswith(os_cpu.lower()):
        version = version[len(os_cpu):]
    cmake_build_type = CMAKE_BUILD_TYPES[version]

    env = setup_build_env(os_cpu, version)
    env['BUILD_DIR'] = build_dir
    get_library_definitions(env)
    setup_python_env(env)
    run_cmake(cmake_build_type, env)

    # Create Makefile
    create_makefile(os_cpu, version, env)

    # Print summary
    print(f"""Variables:
        BUILD_DIR: {build_dir}
        SRC_DIR: {env['SCHRODINGER_SRC']}/rdkit
        OS_CPU: {os_cpu}
        VERSION: {version}""")
    print(f"CMake configuration complete. You can now run 'make all' in {build_dir} to build RDKit.")


if __name__ == "__main__":
    main()
