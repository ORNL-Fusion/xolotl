README.txt for Xolotl Plasma Surface Interactions Simulator

This page provides instructions on how to build Xolotl and its dependencies. All the commands were written for a bash shell, if you are using another type of shell you will have to modify the commands accordingly. CMake is required to build Xolotl, although instructions for it are not included below. We recommend that you download CMake using the package manager on your system. On Red Hat-based systems, this is as simple as running "yum install cmake" or "dnf install cmake" as root.

 You should replace any reference to /home/user to the appropriate prefix for your system. All of the packages are also installed to /opt by preference, but you may adjust that prefix as you see fit too.

Building Xolotl
=====

First things first, checkout the repository with

> git clone https://github.com/ORNL-Fusion/xolotl /home/user/xolotl-source

Now, assuming you have all of the dependencies for Xolotl built and a good C++17 compiler (GCC 5 or greater), Xolotl itself is quite easy to build. If you do not have all of the dependencies, check out the next section.

Create a directory in which you want to build Xolotl (for instance /home/user/xolotl-build) and change into it. It can be any directory except for the source directory (which is /home/user/xolotl-source/). From inside the build directory run the following commands:

>PETSC_DIR=/opt/petsc-install cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=/opt/kokkos-install -DCMAKE_CXX_COMPILER=mpicxx ../xolotl-source \
make \
make test 

with the path to your MPI compiler (here mpicxx) in -DCMAKE_CXX_COMPILER and the path to an installed PETSc version in PETSC_DIR (multiple PETSc versions can coexist on the same file-system).

You can also run make in parallel (make -jN, where N is the number of processes) to build Xolotl faster. This is recommended if you have a big machine.

**NOTE:** If CMake fails to find HDF5, add the path to your HDF5 installation directory in the CMake command:

>HDF5_ROOT=/opt/hdf5 PETSC_DIR=/opt/petsc-install cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=/opt/kokkos-install -DCMAKE_CXX_COMPILER=mpicxx ../xolotl-source 

**NOTE (bis):**If you have Boost  but CMake fails to find BOOST, add the path to your Boost installation folder in the CMake command:

>BOOST_ROOT=/opt/boost PETSC_DIR=/opt/petsc-install cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=/opt/kokkos-install -DCMAKE_CXX_COMPILER=mpicxx ../xolotl-source 


Building Xolotl's Dependencies
=====

The following codes must be built and configured as described below in order for Xolotl to compile. If working on a system where modules are available we recommend loading the existing modules.

MPICH 3
-----

[MPICH 3.0](http://www.mpich.org) or higher is required for Xolotl and its dependencies (alternatively, [Open MPI](https://www.open-mpi.org/) can be used). Use the following commands once downloaded:

>./configure --prefix=/opt/mpich --enable-shared=yes \
make \
make install 

Your path needs to be updated to point to the new MPICH3 install, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

>PATH=/opt/mpich/bin:$PATH \
export PATH

Your library path also needs to be updated in the same file:

>LD_LIBRARY_PATH=/opt/mpich/lib:$LD_LIBRARY_PATH \
export LD_LIBRARY_PATH

Kokkos
-----

[Kokkos](https://github.com/kokkos/kokkos/wiki) is used to enable different back-ends for Xolotl (the develop branch):

###SERIAL###

>git clone https://github.com/kokkos/kokkos /opt/kokkos \
cd /opt/kokkos \
git checkout develop \
mkdir build \
cd build \
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_INSTALL_PREFIX=/opt/kokkos-install -DKokkos_ENABLE_SERIAL=ON -DBUILD_SHARED_LIBS=ON .. \
make install 

###OpenMP###

>git clone https://github.com/kokkos/kokkos /opt/kokkos \
cd /opt/kokkos \
git checkout develop \
mkdir build \
cd build \
cmake -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_INSTALL_PREFIX=/opt/kokkos-install -DKokkos_ENABLE_OPENMP=ON -DBUILD_SHARED_LIBS=ON .. \
make install 

###CUDA###

>git clone https://github.com/kokkos/kokkos /opt/kokkos \
cd /opt/kokkos \
git checkout develop \
mkdir build \
cd build \
cmake -DCMAKE_CXX_COMPILER=/opt/kokkos/bin/nvcc_wrapper -DCMAKE_INSTALL_PREFIX=/opt/kokkos-install -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ENABLE_CUDA_CONSTEXPR=ON -DBUILD_SHARED_LIBS=ON .. \
make install 

The Xolotl cmake line has to be updated here:

>PETSC_DIR=/opt/petsc-install cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_DIR=/opt/kokkos-install -DCMAKE_CXX_COMPILER=/opt/kokkos-install/bin/nvcc_wrapper ../xolotl-source \

Boost 
-----

Xolotl uses [Boost](http://www.boost.org). You can get Boost from your package manager or build it from scratch by downloading the latest. Use the following commands in the Boost directory:

>./bootstrap.sh --show-libraries \
./b2 install --prefix=/opt/boost --with-test --with-filesystem --with-program_options --with-timer

This will only install the parts of Boost required for Xolotl.

If Boost is installed to a non-standard prefix such as the above, you can point CMake to it by using the -DBOOST_ROOT option.

You can run Xolotl tests with

>make test

If you want to only run a specific one, run the corresponding executable that is created in your build folder. For instance (from your build directory)

>./test/core/NetworkTester


PETSc
-----

[PETSc](http://www.mcs.anl.gov/petsc/) is used to solve the DR equation in Xolotl. We build PETSc without debugging enabled. Xolotl now uses the "master" version of PETSc. You can download it (for instance in /opt/) using git with the following command:

git clone https://bitbucket.org/petsc/petsc petsc

Enter the petsc directory and build it using the following commands (hint: after the configure step, PETSc will print you what to execute):

>./configure PETSC_DIR=$PWD --prefix=/opt/petsc-install --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif77 --with-debugging=no  --with-shared-libraries
make PETSC_DIR=/opt/petsc PETSC_ARCH=arch-linux2-c-opt all
make PETSC_DIR=/opt/petsc PETSC_ARCH=arch-linux2-c-opt install
make PETSC_DIR=/opt/petsc-install test
make PETSC_DIR=/opt/petsc-install PETSC_ARCH=arch-linux2-c-opt streams NPMAX=<number of MPI processes you intend to use\>

If you don't have Fortran compilers, you can use --with-fc=0 during the configure.

If PETSc complains because it is not finding Blas/Lapack, you have different options: you can install Blas/Lapack yourself (with yum, dnf, or apt-get for instance; you can also download them from [netlib/blas](http://www.netlib.org/blas/) and [netlib/lapack](http://www.netlib.org/lapack/)) and configure PETSc again, or you can tell PETSc to download it with --download-fblaslapack (if Fortran compilers are availalble) or --download-f2cblaslapack (without Fortran compilers).

**Built in place for developers**

When PETSc versions are changing often, it is faster to build PETSc in place instead of installing it. To do so, simply remove the --prefix option from the ./configure command line, and follow PETSc instructions as before:

>make PETSC_DIR=/opt/petsc PETSC_ARCH=arch-linux2-c-opt all

To update PETSc, go to your PETSc directory, do

>git pull

to get the latest version, and rebuild it with

>make PETSC_DIR=/opt/petsc PETSC_ARCH=arch-linux2-c-opt all

To build Xolotl you now need to specify the PETSC_ARCH variable when running cmake additionally to PETSC_DIR.


HDF5 parallel
-----

[HDF5](http://www.hdfgroup.org/HDF5/) is required for I/O in Xolotl and its dependencies. Use the following commands:

>./configure --prefix=/opt/hdf5 --enable-parallel 
make
make install

Your path needs to be updated to point to the new HDF5 install, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

>PATH=/opt/hdf5/bin:$PATH
export PATH

Your library path also needs to be updated in the same file:

>LD_LIBRARY_PATH=/opt/hdf5/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

Performance Libraries
---------------------
The following are the libraries required IF AND ONLY IF you wish to monitor Xolotl's performance. For more information about the performance infrastructure of Xolotl please refer to the [Performance] section of the wiki. 

###PAPI###
[PAPI](http://icl.cs.utk.edu/papi/index.html), Performance Application Programming Interface, provides access to hardware counter information. In Xolotl, this library is used via GPTL to access hardware performance counters.

**FIRST:**  Please refer to the INSTALL.txt file in the root directory to find the installation instructions for your specific operating system. 

The general installation instructions are as follows (configure and make files are in the src directory):

> ./configure
  make
  make test 
  make fulltest (optional step to run all tests)
  make install-all

Your path may need to be updated to set the variable PAPI_PREFIX to point to the root installation of the PAPI library and include files, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

> PAPI_PREFIX=/usr/local
export PAPI_PREFIX

Note that versions other than PAPI 5.3.0 have been known to cause problems, so it is advised to use this version.


###GPTL 5.3###
The General Purpose Timing Library, or [GPTL](http://jmrosinski.github.io/GPTL/), is a lightweight timing library with support for timing arbitrary regions of code, MPI operations, and collecting hardware performance counter information.

See the GPTL web page and INSTALL file in the root of the GPTL source distribution for instructions on building, testing, and installing GPTL.

To be used for building Xolotl, GPTL's code must be compiled so that the resulting static library can be linked into other shared libraries (such as Xolotl's performance library).  For many systems, this means compiling the library using a flag like '-fPIC'.  To configure GPTL to be built this way, add the -fPIC flag to the first definitions of CFLAGS and FFLAGS within the GPTL macros.make configuration file before building.

If you are NOT using version 5.3 of GPTL two other modifications to GPTL might be required, depending on the distribution and version of MPI you are using. GPTL supports MPICH and its variants out of the box, but requires a slight code modification for use with OpenMPI.  If you are using OpenMPI, add -DMPI_STATUS_SIZE=6 to the definition of CFLAGS in the macros.make file.  Also, if you are using a newer release of OpenMPI, you may need to modify the function signatures in pmpi.c so that read-only arguments are declared with 'const'. For instance, with OpenMPI 1.7.4, the signature for MPI_Send is

  int MPI_Send( const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm );

and GPTL's MPI wrappers will not compile without the const qualifier. Workarounds are to either add the missing const keywords in the GPTL pmpi.c file, or use an older version of OpenMPI such as version 1.6.4 that does not declare MPI functions using const qualifiers.

**NOTE:**  Refer to the Configure Requirements section of [Performance] for further instructions on building Xolotl with GPTL (if GPTL is being built with OPENMP=yes in macros.make).

After installing GPTL, your path needs to be updated to set the variable GPTL_PREFIX to point to the new root installation of the GPTL library and include files, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

> GPTL_PREFIX=/usr/local/gptl-$(REVNO)
export GPTL_PREFIX


Visualization Libraries
---------------------
VTK-m is needed in order to plot Xolotl's output.
 
VTK-m
-----

VTK-m is a toolkit of scientific visualization algorithms for emerging processor architectures. VTK-m supports the fine-grained concurrency for data analysis and visualization algorithms required to drive extreme scale computing by providing abstract models for data and execution that can be applied to a variety of algorithms across many different processor architectures. Here is how to get and install it.
 
> mkdir vtkm \
cd vtkm \
git clone https://gitlab.kitware.com/vtk/vtk-m.git \
mkdir build \
cd build \
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/VTKM -DBUILD_SHARED_LIBS=OFF ../vtk-m \
make -j 4 install 
 
Once installed, all that is needed to enable the compilation of the visualization routines in Xolotl is to add the following flag to your cmake (using the absolute path):
-DVTKM_DIR=/opt/VTKM


Running Xolotl
=====

You should now be able to run Xolotl. Check out our page on [Running Xolotl](https://github.com/ORNL-Fusion/xolotl/wiki/Running-Xolotl).
