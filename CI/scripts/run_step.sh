#!/bin/bash

echo -e "\nCheck gcc and clang compilers\n"
gcc --version
clang --version

case "$1" in 

  configure)

    git config --global --add safe.directory ${GITHUB_WORKSPACE}
   
    case "${GH_JOBNAME}" in
      *"clang"*)
        export CC=clang
        export CXX=clang++
        export OMPI_CC=clang
        export OMPI_CXX=clang++
      ;;
      *)
      ;;
    esac

    cd ${GITHUB_WORKSPACE}/..

    git clone https://gitlab.com/petsc/petsc.git -b v3.21.6 petsc
    cd petsc
    bash ${GITHUB_WORKSPACE}/../xolotl/scripts/build_petsc.sh \
        --prefix=${GITHUB_WORKSPACE}/../install \
        --skip-pull

    cd ${GITHUB_WORKSPACE}/..
    mkdir build
    cd build

    cmake \
        -DCMAKE_PREFIX_PATH=${GITHUB_WORKSPACE}/../install \
        ${EXTRA_CMAKE_ARGS} \
        ${GITHUB_WORKSPACE}
   
    ;;

  build)
    cd ${GITHUB_WORKSPACE}/../build
    make -j4
    ;;

  test)
    cd ${GITHUB_WORKSPACE}/../build
    ctest -VV --label-exclude xolotl.tests.system
    ./test/system/SystemTester -- -t
    ;;

  *)
    echo " Invalid step" "$1"
    exit -1
    ;;
esac
