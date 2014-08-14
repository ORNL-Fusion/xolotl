README.txt for Xolotl Plasma Surface Interactions Simulator, copied directly from
https://sourceforge.net/p/xolotl-psi/wiki/Build%20Configuration/

This page provides instructions on how to build Xolotl and its dependencies. CMake is required to build Xolotl, although instructions for it are not included below. We recommend that you download CMake using the package manager on your system. On Red Hat-based systems, this is as simple as running "yum install cmake" as root.

All of the commands are executed in bash unless stated otherwise. You should replace any reference to /home/user to the appropriate prefix for your system. All of the packages are also installed to /opt by preference, but you may adjust that prefix as you see fit too.

Building Xolotl
===============

Assuming you have all of the dependencies for Xolotl built and a good C++11 compiler (GCC 4.7.2 or greater), Xolotl itself is quite easy to build. If you do not have all of the dependencies, check out the next section.

Create a directory in which you want to build Xolotl and change into it. It can be any directory except for the source directory (which is trunk/xolotl). Most people create a directory called "build" in the source directory (trunk/xolotl/build) or in trunk (trunk/build) and build from there. If your build directory is trunk/xolotl/build, from inside the build directory run the following commands:

cmake .. -DCMAKE_BUILD_TYPE=Release
make

If it is trunk/build run:

cmake ../xolotl/
make

If you want run with custom MPI compilers you can use:

CXX=mpicxx cmake ../xolotl -DCMAKE_BUILD_TYPE=Release

You can also run make in parallel (make -jN, where N is the number of processes) to build Xolotl faster. This is recommended if you have a big machine.

Never use the -DCMAKE_CXX_COMPILER=<C++ compiler> option instead of CXX=<C++ compiler> environment variable. It causes CMake to ignore all of your other options (such as -DBoost_ROOT if you have a special build of Boost).

Building Xolotl's Dependencies
==============================

The following codes must be built and configured as described below in order for Xolotl to compile.

Boost (for tests)
-----------------

Xolotl uses [Boost](http://www.boost.org/) for unit tests. If you don't want to run the tests, then you don't need Boost. However, if you do want to run the tests you should have the latest version of Boost installed for your system. You can get Boost from your package manager or build it from scratch by downloading the latest. Use the following commands in the Boost directory:

>./bootstrap.sh --prefix=/opt/boost-1.55_gcc-4.7.2 --with-libraries=test --show-libraries
./b2 install

This will only install the parts of Boost required for testing.

If Boost is installed to a non-standard prefix such as the above, you can point CMake to it by using the -DBOOST_ROOT option.

MPICH 3.0
---------

[MPICH 3.0](http://www.mpich.org/) is required for Xolotl and its dependencies. Use the following commands:

>./configure --prefix=/opt/mpich-3.0rc1_gcc-4.4.6 --enable-shared=yes 2>&1 | tee c.txt
make 2>&1 | tee m.txt
make install 2>&1 | tee mi.txt

Your path needs to be updated to point to the new MPICH3.0 install, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

>PATH=/opt/mpich-3.0rc1_gcc-4.4.6/bin:$PATH
export PATH

Your library path also needs to be updated in the same file:

>LD_LIBRARY_PATH=/opt/mpich-3.0rc1_gcc-4.4.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

PETSc
-----

[PETSc](http://www.mcs.anl.gov/petsc/) is used to solve the reaction-diffusion equation in Xolotl. We build PETSc without debugging enabled and let it download Hypre, BLAS and LAPACK. Xolotl does not use the "master" or release versions of PETSc because it has special requirements due to its large number of unknowns per grid points. There is a special branch of PETSc available in their repository under the "xolotl" branch. You can download it using git with the following command:

git clone -b xolotl https://bitbucket.org/petsc/petsc petsc-xolotl

Enter the petsc-xolotl directory and build it using the following commands:

>./configure --prefix=/opt/petsc-latest_mpich-3.0.1 --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif77 --with-debugging=no --download-fblaslapack=1 --FOPTFLAGS= --with-shared-libraries=1 --download-hypre=yes --with-debugging=0 --download-superlu_dist --download-parmetis --download-metis
make PETSC_DIR=/home/user/petsc-latest PETSC_ARCH=arch-linux2-c-opt all
make PETSC_DIR=/home/bkj/Programs/petsc-latest PETSC_ARCH=arch-linux2-c-opt install
make PETSC_DIR=/opt/petsc-latest_mpich-3.0rc1 test

Your path needs to be updated to point to the new PETSC_DIR install, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

>PETSC_DIR=/opt/petsc-latest_mpich-3.0rc1_opt
export PETSC_DIR

HDF5
-----

[HDF5](http://www.hdfgroup.org/HDF5/) is required for I/O in Xolotl and its dependencies. Use the following commands:

>./configure --prefix=/opt/hdf5-1.8.10_mpich-3.0rc1_gcc-4.4.6 --enable-parallel 
make
make check
make install
make check-install

Your path needs to be updated to point to the new HDF5 install, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

>PATH=/opt/hdf5-1.8.10_mpich-3.0rc1_gcc-4.4.6/bin:$PATH
export PATH

Your library path also needs to be updated in the same file:

>LD_LIBRARY_PATH=/opt/hdf5-1.8.10_mpich-3.0rc1_gcc-4.4.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

NetCDF
------

The mesh database used by Xolotl requires NetCDF Classic, which is [the C-only version of NetCDF](http://www.unidata.ucar.edu/downloads/netcdf/netcdf-4_2_1_1/index.jsp) configured with the --disable-netcdf-4 flag. It can be built with the following commands:

>./configure --prefix=/opt/netcdf-4.2.1.1-classic_gcc-4.4.6 --disable-netcdf-4
make check install

Your library path also needs to be updated in the same files:

>LD_LIBRARY_PATH=/opt/netcdf-4.1.1.1-classic_gcc-4.4.6/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

MOAB
-----

[MOAB](http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB), short for Mesh-Oriented Database, is used for managing mesh in Xolotl. If it not setup to be configured by default and the Autoconf configuration script needs to be built with "autoreconf" before compiling.

>autorefconf -fi
./configure --prefix=/opt/moab-20121127_mpich-3.0rc1_gcc-4.4.6 --with-mpi=/opt/mpich-3.0rc1_gcc-4.4.6/ --with-hdf5=/opt/hdf5-1.8.10_mpich-3.0rc1_gcc-4.4.6/ --with-netcdf=/opt/netcdf-4.1.1.1-classic_gcc-4.4.6/
make
make check
make install

Performance Libraries
---------------------
The following are the libraries required IF AND ONLY IF you wish to monitor Xolotl's performance.  For more information about the performance infrastructure of Xolotl please refer to the [Performance] section of the wiki. 

###PAPI###
[PAPI](http://icl.cs.utk.edu/papi/index.html), Performance Application Programming Interface, provides access to hardware counter information.  In Xolotl, this library is used via GPTL to access hardware performance counters.

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

If you are NOT using version 5.3 of GPTL two other modifications to GPTL might be required, depending on the distribution and version of MPI you are using.  GPTL supports MPICH and its variants out of the box, but requires a slight code modification for use with OpenMPI.  If you are using OpenMPI, add -DMPI_STATUS_SIZE=6 to the definition of CFLAGS in the macros.make file.  Also, if you are using a newer release of OpenMPI, you may need to modify the function signatures in pmpi.c so that read-only arguments are declared with 'const'.  For instance, with OpenMPI 1.7.4, the signature for MPI_Send is

  int MPI_Send( const void* buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm );

and GPTL's MPI wrappers will not compile without the const qualifier.   Workarounds are to either add the missing const keywords in the GPTL pmpi.c file, or use an older version of OpenMPI such as version 1.6.4 that does not declare MPI functions using const qualifiers.

**NOTE:**  Refer to the Configure Requirements section of [Performance] for further instructions on building Xolotl with GPTL (if GPTL is being built with OPENMP=yes in macros.make).

After installing GPTL, your path needs to be updated to set the variable GPTL_PREFIX to point to the new root installation of the GPTL library and include files, which can be done by adding the following lines in either ~/.bashrc or ~/.bash_profile:

> GPTL_PREFIX=/usr/local/gptl-$(REVNO)
export GPTL_PREFIX

Visualization Libraries
---------------------
The following libraries are needed in order to plot Xolotl's outputs. A wiki page on the visualization will be added soon.

###Mesa###

Mesa is an open-source implementation of the OpenGL specification - a system for rendering interactive 3D graphics. A light version of Mesa-7.5.2 is available [here](https://sourceforge.net/p/xolotl-psi/code/HEAD/tree/trunk/deps/mesa/MesaLib-7.5.2-OSMesa.tar.bz2). Once downloaded, you just have to untar and make it in the repository where you want it

> tar -xf MesaLib-7.5.2-OSMesa.tar.bz2
cd Mesa-7.5.2/
make

After that you simply have to add MESA_PREFIX to your ~/.bashrc or ~/.bash_profile

> MESA_PREFIX=/opt/Mesa-7.5.2/
export MESA_PREFIX

###EAVL###

[EAVL](http://ft.ornl.gov/eavl/) is the Extreme-scale Analysis and Visualization Library. If you have [Git](http://git-scm.com/) you simply have to

> git clone https://github.com/jsmeredith/EAVL EAVL
cd EAVL/

(note: an update this week to EAVL has yet to be reflected in Xolotl's code; this will happen shortly, but in the meantime, after cloning the repository, please revert to an older version with: `git checkout 263414fc4aa416e2376c1af9150c6934146644eb`)

./configure
make
make check

to install it and then add EAVL_PREFIX to your ~/.bashrc or ~/.bash_profile

> EAVL_PREFIX=/opt/EAVL/
export EAVL_PREFIX

Running Xolotl
==============

You should now be able to run Xolotl. Check out our page on [Running Xolotl].

Tricks on Red Hat Enterprise Linux 6
====================================

Getting GCC 4.7.2
-----------------

The default compiler on Red Hat Enterprise Linux (RHEL) 6 is GCC 4.4.7, but Xolotl requires GCC 4.7.2 or greater. RHEL supports GCC 4.7.2 if you have the "Linux Developer" version and you can install it following the instructions on the [Developer Toolset page](https://access.redhat.com/site/documentation/en-US/Red_Hat_Developer_Toolset/1/html/User_Guide/sect-Red_Hat_Developer_Toolset-Install.html).

Most RHEL 6 customers do not have access to the "Linux Developer" version of RHEL 6. The Developer Toolset is open source and can be installed through the CentOS repositories. First, as root, create a file in /etc/yum.repos.d/ called testing-devtools.repo and add the following contents

> [testing-devtools]
name=testing devtools for CentOS $releasever
baseurl=http://people.centos.org/tru/devtools-1.1/6/$basearch/RPMS
gpgcheck=0

Next, also as root, install the Developer Toolset version 1.1 with yum by issuing the following command:

>yum --enablerepo=testing-devtools install devtoolset-1.1

Answer yes at the prompts to install the toolset.

The Developer Toolset and GCC 4.7.2 will be installed to /opt/centos/devtoolset-1.1 and the executables will be available in /opt/centos/devtoolset-1.1/root/usr/bin. This installation will not automatically override the default install of GCC 4.4.7. It will still be on your path and take precedence over 4.7.2, so configure your path or direct CMake to build Xolotl accordingly.

The simplest way to use the new compiler is to rebuild MPI with GCC 4.7.2 (which you kind of have to do anyway) and then rebuild Xolotl with MPI.

Updating Boost
--------------

The default version of Boost on RHEL6 will not work with GCC 4.7.2 because of changes in both packages for the C++11 standard implementation. Download the latest version of Boost from the website and build it with the instructions above. Only the testing libraries are required.

If you build Boost this way, then you will need to add the -DBoost_ROOT=<path to boost> and -DBoost_NO_BOOST_CMAKE=ON to your CMake build configuration as such:

CXX=mpicxx cmake ../xolotl -DCMAKE_BUILD_TYPE=Release -DBoost_ROOT=/opt/boost -DBoost_NO_BOOST_CMAKE=ON

