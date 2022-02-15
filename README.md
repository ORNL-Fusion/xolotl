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
- Kokkos >= 3.2
- HDF5
- PETSc >= 3.14

## Quick-Start Install
If you have the above dependencies installed
```
git clone https://github.com/ORNL-Fusion/xolotl $HOME/xolotl-source
mkdir $HOME/xolotl-build
cd $HOME/xolotl-build
cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DKokkos_DIR=/opt/kokkos-install \
    -DPETSC_DIR=/opt/petsc-install \
    -DCMAKE_INSTALL_PREFIX=$HOME/xolotl-build/install \
    $HOME/xolotl-source/
make
make test
make install
```

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
description of the available options at [Parameter
File](https://github.com/ORNL-Fusion/xolotl/wiki/Parameter-File).

## More Info
If you want to contribute please check
[Guidelines](https://github.com/ORNL-Fusion/xolotl/wiki/Guidelines).

For any question you can send an email to our mailing list:
xolotl-psi-development@lists.sourceforge.net.
