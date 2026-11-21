#!/usr/bin/env python
"""
Setup script for zol. Package metadata and dependencies are defined in pyproject.toml.

When installing from source inside a conda environment (e.g. `pip install .` after creating the environment from
zol_env.yml), this script additionally compiles zol's small C++ helper programs into the environment's bin/
directory and sets up activation scripts defining ZOL_DATA_PATH / ZOL_EXEC_PATH. These extra steps are skipped when
building via conda-build (e.g. bioconda), whose recipe handles them itself, or when not in a conda environment.
Problems in these extra steps never cause the Python installation to fail; they are reported in
installation_error.txt instead.
"""
import os
import subprocess
import traceback

from setuptools import setup

CPP_PROGRAMS = {
    "runRBH": "src/zol/orthologs/runRBH.cpp",
    "splitDiamondResults": "src/zol/orthologs/splitDiamondResults.cpp",
}


def conda_post_install_setup(conda_prefix):
    conda_dir = os.path.abspath(conda_prefix)
    zol_exec_path = os.path.join(conda_dir, "bin")
    zol_data_path = os.path.join(conda_dir, "share", "zol", "db")
    conda_act_dir = os.path.join(conda_dir, "etc", "conda", "activate.d")
    conda_dea_dir = os.path.join(conda_dir, "etc", "conda", "deactivate.d")

    compiler = os.environ.get("CXX", "g++")
    for program, source in CPP_PROGRAMS.items():
        subprocess.run(
            [compiler, "-std=c++14", "-o", os.path.join(zol_exec_path, program), source],
            check=True,
        )
        assert os.path.isfile(os.path.join(zol_exec_path, program))

    for directory in (zol_data_path, conda_act_dir, conda_dea_dir):
        os.makedirs(directory, exist_ok=True)

    with open(os.path.join(zol_data_path, "README.txt"), "w") as readme_handle:
        readme_handle.write("Default conda space for downloading annotation databases.\n")

    with open(os.path.join(conda_act_dir, "zol.sh"), "w") as act_handle:
        act_handle.write("#!/usr/bin/env bash\n\n")
        act_handle.write(f"export ZOL_DATA_PATH={zol_data_path}\n")
        act_handle.write(f"export ZOL_EXEC_PATH={zol_exec_path}\n")

    with open(os.path.join(conda_dea_dir, "zol.sh"), "w") as dea_handle:
        dea_handle.write("#!/usr/bin/env bash\n\n")
        dea_handle.write("unset ZOL_DATA_PATH\n")
        dea_handle.write("unset ZOL_EXEC_PATH\n")


if __name__ == "__main__":
    setup()

    conda_prefix = os.environ.get("CONDA_PREFIX")
    running_conda_build = os.environ.get("CONDA_BUILD") == "1"
    if conda_prefix and not running_conda_build:
        try:
            conda_post_install_setup(conda_prefix)
        except Exception:
            with open("installation_error.txt", "w") as error_handle:
                error_handle.write(f"conda dir is: {conda_prefix}\n")
                error_handle.write(
                    "C++ compilation or setup of executables failed or had issues setting up space in conda "
                    "environment for databases!\n"
                )
                error_handle.write(traceback.format_exc())
