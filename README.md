# Xolotl
Xolotl is a continuum simulator under development within two DOE SciDAC
programs, available at https://github.com/ORNL-Fusion/xolotl . It is a
high-performance computing code using advection-reaction-diffusion (ADR) kinetic
rate theory to model the time evolution of the divertor material in next
generation tokamaks like ITER, as well as nuclear fuel in fission reactors.

## Required Dependencies
- C++17 compiler
- CMake >= 3.13
- MPI
- Boost
- Kokkos >= 4.6
- HDF5
- PETSc >= 3.23 (configured with kokkos and kokkos-kernels; see
  scripts/build_petsc.sh)

## Quick-Start Install
Xolotl can now build all dependencies (using PETSc's build system) during the
CMake configure step:
```
cd $HOME
git clone https://github.com/ORNL-Fusion/xolotl $HOME/xolotl-source
mkdir $HOME/xolotl-build
cd $HOME/xolotl-build
cmake \
    -DXolotl_BUILD_PETSC=ON \
    [-DXolotl_ENABLE_CUDA=ON | -DXolotl_ENABLE_OPENMP=ON] \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$HOME/xolotl-build/install \
    $HOME/xolotl-source/
make
make test
make install
```

If Boost, HDF5, or LAPACK are not found during configuration, they will
automatically be requested from PETSc's build system. Or you can specify this
directly via `-DXolotl_BUILD_BOOST=ON`, *etc.*

If you wish to configure and install PETSc separately, use (or reference) our
[PETSc build script](https://github.com/ORNL-Fusion/xolotl/wiki/Build-Configuration#petsc-build-script).

## Wiki
You can find more information about this code on our
[Wiki](https://github.com/ORNL-Fusion/xolotl/wiki) and how to install Xolotl and
its dependencies in [Build
Configuration](https://github.com/ORNL-Fusion/xolotl/wiki/Build-Configuration). 

Once installed, you can use a
[Benchmark](https://github.com/ORNL-Fusion/xolotl/wiki/Benchmark-Problems) file
to run the code:

```
./xolotl <input_parameter_file> 
```

as described in [Running
Xolotl](https://github.com/ORNL-Fusion/xolotl/wiki/Running-Xolotl), with the
description of the different options at [Parameter
File](https://github.com/ORNL-Fusion/xolotl/wiki/Parameter-File) and 
[PETSc Options](https://github.com/ORNL-Fusion/xolotl/wiki/PETSc-Options).

## More Info
If you want to contribute please check
[Guidelines](https://github.com/ORNL-Fusion/xolotl/wiki/Guidelines).

A list of publications using Xolotl is available at 
[Publications](https://github.com/ORNL-Fusion/xolotl/wiki/Publications). 

For any question you can send an email to our mailing list:
xolotl-psi-development@lists.sourceforge.net.
