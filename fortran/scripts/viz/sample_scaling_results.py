#!/usr/bin/env python

from collections import namedtuple
from os import listdir
from os.path import join, realpath, basename, getmtime
from pathlib import Path
from re import findall
from typing import NamedTuple

# import matplotlib.pyplot as plt

FORTRAN_PROJECT: Path = Path(realpath(__file__)).parents[2]
BUILD_OUTPUT_DIR: str = join(
    FORTRAN_PROJECT, "build", "build_output")


def main():
    # get the subdirectories corresponding to the experiments
    subdirs: list[str] = listdir(BUILD_OUTPUT_DIR)
    experiment_subdirs: list[str] = [
        join(BUILD_OUTPUT_DIR, dir) for dir in subdirs if "nx_" in dir]

    # extract parameters based on these subdirectories
    # NOTE: you could adapt this based on your naming conventions
    subdir_name_regex = r"nx_(\d*)_nz_(\d*)_sim_time_(\d*)_out_freq_(\d*)"
    SimParams = namedtuple(
        "SimParams", (
            "nx",
            "nz",
            "sim_time",
            "out_freq"
        )
    )
    sim_params_set: list[SimParams] = []
    for subdir in experiment_subdirs:
        find_sim_params: list[tuple[str]] = findall(
            subdir_name_regex, subdir)
        assert len(find_sim_params) == 1, "must find only one match"
        sim_params: SimParams = SimParams(
            *tuple(map(lambda p: int(p), find_sim_params[0])))
        sim_params_set.append(sim_params)

    # extract measurements from output directories and build object holding
    # experiment simulation parameters and measurements
    # NOTE: You could adapt this to your measurements
    Measure = namedtuple(
        "Measure", (
            "nthreads",
            "cpu_time"
        )
    )
    Experiment = namedtuple(
        "Experiment", (
            "sim_params",
            "measures"
        )
    )

    experimental_output_dir_regex = r"output_nthread_(\d*)"
    experiments: list[Experiment] = []
    for subdir, sim_params in zip(experiment_subdirs, sim_params_set):

        experimental_output_dirs = [
            join(subdir, d) for d in listdir(subdir) if "output_" in d]
        measures: list[Measure] = []

        for experimental_output_dir in experimental_output_dirs:
            experimental_output_dir_basename = basename(
                experimental_output_dir)
            find_nthreads = findall(
                experimental_output_dir_regex,
                experimental_output_dir_basename)
            assert len(find_nthreads) == 1, "must find only one match"
            nthreads = int(find_nthreads[0])
            log_files = get_mtime_sorted_log_files(experimental_output_dir)
            most_recent_log_file = log_files[-1]
            cpu_time = get_cpu_time(most_recent_log_file)
            measure = Measure(nthreads, cpu_time)
            measures.append(measure)
        experiment = Experiment(sim_params, measures)
        print(experiment)
        print()
        experiments.append(experiment)

    # TODO:
    # parse expereimnt outputs and plot

    return


def get_mtime_sorted_log_files(experimental_output_dir: str) -> list[str]:
    log_files: list[str] = ["dummy"]

    return log_files


def get_cpu_time(log_file: str) -> float:
    cpu_time: float = 0.0
    return cpu_time


if __name__ == "__main__":
    main()
