# The MiniWeather App

It is a miniature fluid-dynamics code solving the 2-D inviscid Euler equations for stratified fluids. The dynamics themselves are dry compressible, stratified, non-hydrostatic flows dominated by buoyant forces that are relatively small perturbations on a hydrostatic background state. The code uses periodic boundary conditions in the x-direction and solid wall boundary conditions in the z-direction. The MiniWeather App is designed for training in parallel HPC computing. Parallelization approaches currently include:

* MPI (C, Fortran, and C++)
* OpenMP Threading (C and Fortran)
* OpenACC Offload (C and Fortran)

Original Author: Matt Norman, Oak Ridge National Laboratory,
https://mrnorman.github.io

Author for DKRZ Levante: Jared Frazier, Leibniz Institute of Atmospheric
Physics, https://jfdev001.github.io/

# Introduction
There are four main directories in MiniWeather: (1) a Fortran source directory; (2) a C source directory; (3) a C++ source directory; and (4) a documentation directory. We focus on MPI and Open MP in Fortran code, but you can find information on C and C++ and on OpenACC here: https://github.com/mrnorman/miniWeather

## Fluid State Variables

There are four main arrays used in this code: `state`, `state_tmp`, `flux`, and `tend`, and the dimensions for each are given in the code upon declaration in the comments. Each of these arrays is described briefly below:

* `state`: This is the fluid state at the current time step, and it is the only array that persists from one time step to the next. The other four are only used within the calculations to advance the model to the next time step. The fluid state describes the average state over each cell area in the spatial domain. This variable contains four fluid states, which are the traditional mass, momenta, and thermodynamic quantities of most fluid models:
  1. Density (`ID_DENS`): The 2-D density of the fluid, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho" title="\large \rho" />, in <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-2}" title="\large \text{kg}\ \text{m}^{-2}" /> (note this is normally <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-3}" title="\large \text{kg}\ \text{m}^{-3}" />, but this is a 2-D model, not 3-D)
  2. U-momentum (`ID_UMOM`): The momentum per unit area of the fluid in the x-direction calculated as <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho&space;u" title="\large \rho u" />, where u is the x-direction wind velocity. The units are <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-1}\&space;\text{s}^{-1}" title="\large \text{kg}\ \text{m}^{-1}\ \text{s}^{-1}" />. Note that to get true momentum, you must integrate over the cell.
  2. W-momentum (`ID_WMOM`): The momentum per unit area of the fluid in the z-direction calculated as <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho&space;w" title="\large \rho w" />, where w is the z-direction wind velocity. The units are <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{kg}\&space;\text{m}^{-1}\&space;\text{s}^{-1}" title="\large \text{kg}\ \text{m}^{-1}\ \text{s}^{-1}" />. Note that to get true momentum, you must integrate over the cell.
  4. Potential Temperature (`ID_RHOT`): The product of density and potential temperature, <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\rho&space;\theta" title="\large \rho \theta" />, where <img src="https://latex.codecogs.com/svg.latex?\theta=T\left(P_{0}/P\right)^{R_{d}/c_{p}}" /> , <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;P_{0}=10^{5}\,\text{Pa}" title="\large P_{0}=10^{5}\,\text{Pa}" />, T is the true temperature, and <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;R_d" title="\large R_d" /> and<img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;c_p" title="\large c_p" /> are the dry air constant and specific heat at constant pressure for dry air, respectively. The units of this quantity are <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{K}\,\text{kg}\,\text{m}^{-2}" title="\large \text{K}\,\text{kg}\,\text{m}^{-2}" />.
* `state_tmp`: This is a temporary copy of the fluid state used in the Runge-Kutta integration to keep from overwriting the state at the beginning of the time step, and it has the same units and meaning.
* `flux`: This is fluid state at cell boundaries in the x- and z-directions, and the units and meanings are the same as for `state` and `state_tmp`. In the x-direction update, the values of `flux` at indices `i` and `i+1` represents the fluid state at the left- and right-hand boundaries of cell `i`. The indexing is analogous in the z-direction. The fluxes are used to exchange fluid properties with neighboring cells.
* `tend`: This is the time tendency of the fluid state <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\partial\mathbf{q}/\partial&space;t" title="\large \partial\mathbf{q}/\partial t" />, where <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\mathbf{q}" title="\large \mathbf{q}" /> is the the state vector, and as the name suggests, it has the same meaning and units as state, except per unit time (appending <img src="https://latex.codecogs.com/svg.latex?\inline&space;\dpi{300}&space;\large&space;\text{s}^{-1}" title="\large \text{s}^{-1}" /> to the units). In the Finite-Volume method, the time tendency of a cell is equivalent to the divergence of the flux across a cell.

## Numerical Experiments

A number of numerical experiments are in the code for you to play around with. You can set these by changing the `data_spec_int` variable. 

### Rising Thermal

```
data_spec_int = DATA_SPEC_THERMAL
sim_time = 1000
```

This simulates a rising thermal in a neutral atmosphere, which will look something like a “mushroom” cloud (without all of the violence).

Potential Temperature after 500 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/thermal_pt_0500.png" width=400/>

Potential Temperature after 1,000 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/thermal_pt_1000.png" width=400/>

### Colliding Thermals

```
data_spec_int = DATA_SPEC_COLLISION
sim_time = 700
```

This is similar to the rising thermal test case except with a cold bubble at the model top colliding with a warm bubble at the model bottom to produce some cool looking eddies.

Potential Temperature after 200 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/collision_pt_0200.png" width=400/>

Potential Temperature after 400 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/collision_pt_0400.png" width=400/>

Potential Temperature after 700 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/collision_pt_0700.png" width=400/>

### Mountain Gravity Waves

```
data_spec_int = DATA_SPEC_MOUNTAIN
sim_time = 1500
```

This test cases passes a horizontal wind over a faked mountain at the model bottom in a stable atmosphere to generate a train of stationary gravity waves across the model domain.

Potential Temperature after 400 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/mountain_pt_0400.png" width=400/>

Potential Temperature after 1,300 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/mountain_pt_1300.png" width=400/>

### Density Current

```
data_spec_int = DATA_SPEC_DENSITY_CURRENT
sim_time = 600
```

This test case creates a neutrally stratified atmosphere with a strong cold bubble in the middle of the domain that crashes into the ground to give the feel of a weather front (more of a downburst, I suppose).

Potential Temperature after 200 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/density_current_pt_0200.png" width=400/>

Potential Temperature after 600 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/density_current_pt_0600.png" width=400/>

### Injection

```
data_spec_int = DATA_SPEC_INJECTION
sim_time = 1200
```

A narrow jet of fast and slightly cold wind is injected into a balanced, neutral atmosphere at rest from the left domain near the model top. This has nothing to do with atmospheric flows. It's just here for looks. 

Potential Temperature after 300 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/injection_pt_0300.png" width=400/>

Potential Temperature after 1,000 seconds:

<img src="https://github.com/mrnorman/miniWeather/blob/main/documentation/images/injection_pt_1000.png" width=400/>

# Compiling and Running the Code

## Software Dependencies

* Parallel-netcdf: https://github.com/Parallel-NetCDF/PnetCDF
  * This is a dependency for two reasons: (1) NetCDF files are easy to visualize and convenient to work with; (2) The users of this code shouldn't have to write their own parallel I/O.
* Ncview: http://meteora.ucsd.edu/~pierce/ncview_home_page.html
  * This is the easiest way to visualize NetCDF files.
* MPI
* For OpenACC: An OpenACC-capable compiler (PGI / Nvidia, Cray, GNU)
  * A free version of the PGI / Nvidia compiler can be obtained by googling for the "Community Edition"
* For OpenMP: An OpenMP offload capable compiler (Cray, XL, GNU)
* For C++ portability, Nvidia's CUB and AMD's hipCUB and rocPRIM are already included as submodules
* CMake: https://cmake.org

## Basic Setup

```shell
cd /work/bm1233/${USER}  
git clone git@github.com:jfdev001/miniWeather.git
MINIWEATHER_DIR=$(pwd)/miniWeather
cd ${MINIWEATHER_DIR}
git submodule update --init --recursive
```

To find that repository on GitHub, go to  

```text
https://github.com/jfdev001/miniWeather
```

and star it so that you can easily find it later.

If you prefer, you can fork (see [github docs: fork a
repo](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo))
that repo and then clone your own fork of `miniweather`. This is also a good
approach since then you can upload (i.e., `git push`) your code to a repository
on your GitHub. The workflow would look something like the following

```shell
cd /work/bm1233/${USER}  
# assuming you have forked miniweather
git clone git@github.com:YOUR_GITHUB_USER_NAME_HERE/miniWeather.git 
MINIWEATHER_DIR=$(pwd)/miniWeather
cd ${MINIWEATHER_DIR}
git submodule update --init --recursive

# by default, the remote origin (i.e., source of your code on GitHub
# and its destination when pushing) is set to 
# git@github.com:YOUR_GITHUB_USER_NAME_HERE/miniWeather.git ...
# by adding another remote and calling it upstream, you can at any
# time inspect or pull the code provided by jfdev001 without breaking
# you local changes 
git remote add upstream git@github.com:jfdev001/miniWeather.git 
git fetch upstream  # allows you to pull code from jfdev001 in the future
```

## Building the Code and Testing the Workflow

There are four main directories in the mini app: (1) a Fortran source
directory; (2) a C source directory; (3) a C++ source directory; and (4) a
documentation directory. 

`miniWeather` uses the [CMake](https://cmake.org/) build system. While the
focus of this course is the code in the `fortran` directory, you may
look at the `c` and `cpp` directories if you wish to find build scripts matching
different HPC systems and compiler setups. Those build scripts in the `c`
and `cpp` directories will *not* work on Levante and will have to be modified
if you wish to run those codes.

Note that you must source the cmake scripts in the `build/` directores because
they do module loading and set a `TEST_MPI_COMMAND` environment variable
because it will differ from machine to machine.

The first thing you should do is verify that you can compile and run
`miniweather`:

```shell
bash ${MINIWEATHER_DIR}/fortran/build/cmake_levante_test
```

This generates a directory called `${MINIWEATHER_DIR}/fortran/build/build_output/test`
where all configuration (e.g., auto-generated Makefiles) and compilation
artifacts (e.g., executable binaries like `serial`, `openmp`, and `mpi`) are placed.

You should *always* read the usage documentation for any script you run. For
nearly every script provided, you can do the following to get usage
documentation

```shell
# assuming you're in the directory of the script...
# also if you encounter permissions issues, call `chmod +x <script_name>`
./<script_name> -h
```

Pay close attention also to any examples section in usage documentation.

If you are using MobaXTerm you should have two ssh sessions open and connected
to Levante: (1) one just for looking at the uage documentation of scripts so
that you know what the inputs and outputs are, and (2) one for actually
launching scripts and/or editing files. If you are Linux or Mac, you can also
launch multiple terminals and each of them can independently ssh to Levante.

You can check to see what `cmake_levante_test` is doing by typing

```shell
bash ${MINIWEATHER_DIR}/fortran/build/cmake_levante_test -h
```

Note that the `cmake_levante_test` simply wraps the
`cmake_levante_config_and_build` script discussed in the following sections.

In the real world, there may be limited or *no* documentation for software that
you are using. Even worse, there may be documentation but it could be *out of
date*. This is almost worse than having no documentation because you might
think the software is doing one thing while it is in reality doing something
completely unexpected. You need to be prepared to *read* through code to
determine what it does. When doing this, keep mind the following:

* What are the inputs to the code? 
    * Does it take positional arguments? Does it take flags? 
    * Does it have optional arguments? What are the defaults of those optional
    arguments? 
    * Does it read from a file?
    * What assumptions (if any) about paths (e.g., directories) does the code
    make?

* What are the outputs of the code?
    * Does it write to standard output (i.e., the terminal)?
    * Does it produce large files (e.g., netcdf)? Where do those files get 
    written?

* Does the script launch Slurm jobs or execute locally?


## Altering the Code's Configurations

To alter the configuration of the code, you can control the number of cells in
the x- and z-directions, the length of simulation time, the output frequency,
and the initial data to use by calling `cmake_levante_build_and_configure`,
see 

```shell
# assuming in build/ dir
bash cmake_levante_config_and_build -h
```

This script forwards arguments to two calls to `cmake` that configure and build
the script. The generated `cmake` configuration call might look like the below:

* `-DNX=400`: Uses 400 cells in the x-direction
* `-DNZ=200`: Uses 200 cells in the z-direction
* `-DSIM_TIME=1000`: Simulates for 1,000 seconds model time
* `-DOUT_FREQ=10`: Outputs every 10 seconds model time
* `-DDATA_SPEC=DATA_SPEC_THERMAL`: Initializes a rising thermal

It's best if you keep `NX` exactly twice the value of `NZ` since the domain is
20km x 10km. 

The data specifications are `DATA_SPEC_COLLISION`, `DATA_SPEC_THERMAL`,
`DATA_SPEC_MOUNTAIN`, `DATA_SPEC_DENSITY_CURRENT`, and `DATA_SPEC_INJECTION`,
as described above.

## Running the Code

It is recommended to run the `build_output/test/serial_test` code first to get
an idea for the model outputs. Note that a file `output.nc` is always produced 
in the directory from which you call the `miniweather` executables.

As an example:

```shell
# assuming in fortran/ directory... this produces `output.nc` there
cd ${MINIWEATHER_DIR}/fortran
./build/build_output/test/serial_test
```

Since parameters are set in the code itself, you don't need to pass any
parameters to the executables.

This is fine for testing lightweight serial codes; however, we are interested
in parallel codes. Since you are using a compute cluster shared by many people,
jobs requiring more computational resources must be submitted to the Slurm
scheduler.

We provide a script that wraps the generation of a run script which you can
use later for running simulations. You can check out the parameters for this
script here:

```shell
bash ${MINIWEATHER_DIR}/fortran/scripts/templates/make_run_scripts -h
```

This script can be used to generate Slurm scripts specific to your user for
running `miniweather` simulations. These scripts are, by convention, written to
`scripts/run` and are *not* tracked by `git`. If you wish to modify the
`.gitignore` file and remove the line containing `*.run`, `git` will track
your generated run scripts. The run scripts will also be prefixed with the
partition that you have requested. Different partitions on levante (e.g.,
shared, compute, gpu) give the user differ compute resources. By default the
compute partition is used, and this will work with MPI and OpenMP jobs.

You should generate an example run script with the following:

```shell
EMAIL_HERE="put_your_email@gmail.com"
bash ${MINIWEATHER_DIR}/fortran/scripts/templates/make_run_scripts ${EMAIL_HERE}
```

This generates `scripts/run/compute_miniweather.run`. You should inspect what
this script does with

```shell
bash ${MINIWEATHER_DIR}/fortran/scripts/run/compute_miniweather.run -h
```

In particular, you should run each of the `bash` examples in the usage doc to
get an understanding of what *would* be submitted to Slurm. There are lots of
outputs so make sure to read and understand them. When ready to submit jobs to
Slurm, look at the `sbatch` examples in the same usage doc.

Note that the `time` sbatch directive is set to 30 seconds. This is sufficient
for running tests, but may not be sufficient for running larger scale
simulations. If your simulation significantly exceeds the amount of time
allocated by the `time` directive, the simulation will timeout and you the
outputs to `output.nc` may be incomplete. Be aware that increasing the amount
of time you would like to run your job may result in you waiting longer for
the Slurm scheduler to actually launch your job. You should always prototype
any experiments or scripts that you write which involve Slurm such that they
request a very short amount of time (i.e., less than 1 minute).

## Viewing the Output

The file I/O is done in the netCDF format: (https://www.unidata.ucar.edu/software/netcdf). To me, the easiest way to view the data is to use a tool called “ncview” (http://meteora.ucsd.edu/~pierce/ncview_home_page.html). To use it, you can simply type `ncview output.nc`, making sure you have X-forwarding enabled in your ssh session. Further, you can call `ncview -frames output.nc`, and it will dump out all of your frames in the native resolution you're viewing the data in, and you can render a movie with tools like `ffmpeg`. 

# Parallelization

The code is decomposed in two spatial dimensions, x and z, with `nx_glob` and `nz_glob` cells in the global domain and nx and nz cells in the local domain, using straightforward domain decomposition for MPI-level parallelization. The global domain is of size xlen and zlen meters, and hs “halo” cells are appended to both sides of each dimension.

This code was designed to parallelize with MPI first and then OpenMP, OpenACC, OpenMP offload, or `parallel_for` next, but you can always parallelize with OpenMP or OpenACC without MPI if you want. But it is rewarding to be able to run it on multiple nodes at higher resolution for more and sharper eddies in the dynamics.

As you port the code, you'll want to change relatively little code at a time, re-compile, re-run, and look at the output to see that you're still getting the right answer. There are advantages to using a visual tool to check the answer (e.g., `ncview`), as it can sometimes give you clues as to why you're not getting the right answer. 

Note that you only need to make changes code within the first 450 source lines for C and Fortran, and each loop that needs threading is decorated with a `// THREAD ME` comment. Everything below that is initialization and I/O code that doesn't need to be parallelized (unless you want to) for C and Fortran directives-based approaches.

For the C++ code, you will need to work with the initialization and File I/O code, but everything you need to do is explicitly guided via comments in the code.

## Indexing

The code makes room for so-called “halo” cells in the fluid state. This is a common practice in any algorithm that uses stencil-based reconstruction to estimate variation within a domain. In this code, there are `hs` halo cells on either side of each spatial dimension, and I pretty much hard-code `hs=2`.

In the Fortran code's fluid state (`state`), the x- and z-dimensions are dimensioned as multi-dimensional arrays that range from `1-hs:nx+hs`. In the x-direction, `1-hs:0` belong to the MPI task to the left, cells `1:nx` belong to the current MPI task, and `nx+1:nx+hs` belong to the MPI task to the right. In the z-dimension, `1-hs:0` are artificially set to mimic a solid wall boundary condition at the bottom, and `nz+1:nz+hs` are the same for the top boundary. The cell-interface fluxes (`flux`) are dimensioned as `1:nx+1` and `1:nz+1` in the x- and z-directions, and the cell average tendencies (`tend`) are dimensioned `1:nx` and `1:nz` in the x- and z-directions. The cell of index `i` will have left- and right-hand interface fluxes of index `i` and `i+1`, respectively, and it will be evolved by the tendency at index `i`. The analog of this is also true in the z-direction.

## MPI Domain Decomposition

This code was designed to use domain decomposition, where each MPI rank “owns” a certain set of cells in the x-direction and contains two “halo” cells from the left- and right-hand MPI tasks in the x-direction as well. The domain is only decomposed in the x-direction and not the z-direction for simplicity.

**IMPORTANT**: Please be sure to set `nranks`, `myrank`, `nx`, `i_beg`, `left_rank`, and `right_rank`. These are clearly marked in the serial source code. You can set more variables, but these are used elsewhere in the code (particularly in the parallel file I/O), so they must be set.

To parallelize with MPI, there are only two places in the code that need to be altered. The first is the initialization, a subroutine / function named `init`, where you must determine the number of ranks, you process's rank, the beginning index of your rank's first cell in the x-direction, the number of x-direction cells your rank will own, and the MPI rank IDs that are to your left and your right. Because the code is periodic in the x-direction, your left and right neighboring ranks will wrap around. For instance, if your are rank `0`, your left-most rank will be `nranks-1`.

The second place is in the routine that sets the halo values in the x-direction. In this routine, you need to:

1. Create MPI data buffers (at the same place the other arrays are declared) to hold the data that needs to be sent and received, allocate them in the `init()` routine, and deallocate them in the `finalize()` routine.

2. Pack the data you need to send to your left and right MPI neighbors

3. Send the data to your left and right MPI neighbors

4. Receive the data from your left and right MPI neighbors

5. Unpack the data from your left and right neighbors and place the data into your MPI rank's halo cells. 

Once you complete this, the code will be fully parallelized in MPI. Both of the places you need to add code for MPI are marked in the serial code, and there are some extra hints in the `set_halo_values_x()` routine as well.


# MiniWeather Model Scaling
## Details to Consider

If you want to do scaling studies with miniWeather, this section will be important to make sure you're doing an apples-to-apples comparison.

* `sim_time`: The `sim_time` parameter does not mean the wall time it takes to simulate but rather refers amount of model time simulated. As you increase `sim_time`, you should expect the walltime to increase linearly.
* `nx_glob, nz_glob`: As a rule, it's easiest if you always keep `nx_glob = nz_glob * 2` since the domain is always 20km x 10km in the x- and z-directions. As you increase `nx_glob` (and proportionally `nz_glob`) by some factor `f`, the time step automatically reduced by that same factor, `f`. Therefore, increasing `nx_glob` by 2x leads to 8x more work that needs to be done. Thus, with the same amount of parallelism, you should expect a 2x increase in `nx_glob` and `nz_glob` to increase the walltime by 8x (neglecting parallel overhead concerns).
  * More precisely, the time step is directly proportional to the minimum grid spacing. The x- and y-direction grid spacings are: `dx=20km/nx_glob` and `dz=10km/nz_glob`. So as you decrease the minimum grid spacing (by increasing `nx_glob` and/or `nz_glob`), you proportionally decrease the size of the time step and therefore proportionally increase the number of time steps you need to complete the simulation (thus proportionally increasing the expected walltime).
* The larger the problem size, `nx_glob` and `nz_glob`, the lower the relative parallel overheads will be. You can get to a point where there isn't enough work on the accelerator to keep it busy and / or enough local work to amortize parallel overheads. At this point, you'll need to increase the problem size to see better scaling. This is a typical [Amdahl's Law](https://en.wikipedia.org/wiki/Amdahl%27s_law) situation.

Remember that you can control each of these parameters through the CMake configure.

## Running Performance Experiments

You may want to evaluate how the performance of `miniweather` is affected by
increasing the number of threads, increasing the number of MPI processes, or
doing a combination of both. You can inspect a sample bash script that prepares
and launches such experiments:

```shell
# assuming in the fortran/ directory
./scripts/scaling/launch_sample_scaling_experiments -h
```

You can use that script as a template for running your own experiments.

## Visualizing Performance Results

This will also depend heavily on the types of experiments that you wish to run,
however, an example python code that can be launched by:

```shell
# assuming in the fortran/ directory
python scripts/viz/sample_scaling_results.py
```

That script has no `-h` option supported; however, at the top of the file
is a small description of the contents of the script itself and what it's for.

You copy/modify it to accomplish your plotting goals for your experiments.

Below is an example output from the script:

<img width="999" height="799" alt="miniweather_openmp" src="https://github.com/user-attachments/assets/5f2959bf-393a-4ae2-8008-67383dffcc01" />



