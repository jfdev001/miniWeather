#!/usr/bin/env python

from collections import namedtuple
import os
from pathlib import Path
import re

import matplotlib.pyplot as plt

FORTRAN_PROJECT: Path = Path(os.path.realpath(__file__)).parents[2]
BUILD_OUTPUT_DIR: str = os.path.join(
    FORTRAN_PROJECT, "build", "build_output")


def main():
    # get the subdirectories corresponding to the experiments
    subdirs: list[str] = os.listdir(BUILD_OUTPUT_DIR)
    experiment_subdirs: list[str] = [
        os.path.join(BUILD_OUTPUT_DIR, dir) for dir in subdirs if "nx_" in dir]

    # extract parameters based on these subdirectories
    subdir_name_regex = r"nx_(\d*)_nz_(\d*)_sim_time_(\d*)_out_freq_(\d*)"
    SimParams = namedtuple("SimParams", ["nx", "nz", "sim_time", "out_freq"])
    sim_params_set: list[SimParams] = []
    for subdir in experiment_subdirs:
        init_sim_params: list[tuple[str]] = re.findall(
            subdir_name_regex, subdir)
        assert len(init_sim_params) == 1
        sim_params: tuple[int] = SimParams(
            *tuple(map(lambda p: int(p), init_sim_params[0])))
        sim_params_set.append(sim_params)

    print(sim_params_set)

    # extract more parameters based on output directory naming

    # go to the log file of each and get CPU Time: ....

    pass


if __name__ == "__main__":
    main()
