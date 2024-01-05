#!/bin/bash

# Make sure we're sane
set -eu -o pipefail

# Variables
_prefix=$HOME/.local
_dry_run=0
_do_cleanup=1
_do_pull=1
_debug=0
_use_cuda=0
_use_omp=0
_petsc_extra_args=""
_petsc_dir=$PWD
_petsc_dir_arch_set=""
_opt_flags="-O3"

# Read command-line arguments
while [ $# -gt 0 ]
do
    case $1 in
    --dry-run)
        _dry_run=1
        ;;
    --skip-cleanup)
        _do_cleanup=0
        ;;
    --skip-pull)
        _do_pull=0
        ;;
    --prefix=*)
        _prefix="${1:9}" # strip "--prefix="
        ;;
    --petsc-dir=*)
        _petsc_dir="${1:12}" # strip "--petsc-dir="
        _petsc_dir_arch_set="${_petsc_dir_arch_set} PETSC_DIR=${_petsc_dir}"
        ;;
    --petsc-arch=*)
        _petsc_arch="${1:13}" # strip "--petsc-arch="
        _petsc_dir_arch_set="${_petsc_dir_arch_set} PETSC_ARCH=${_petsc_arch}"
        ;;
    --make-np=*)
        _np="${1:10}" # strip "--make-np="
        _petsc_extra_args="${_petsc_extra_args} --with-make-np=${_np}"
        ;;
    --debug|--dbg)
        _debug=1
        _opt_flags="-O0"
        ;;
    --cuda)
        _use_cuda=1
        ;;
    --openmp)
        _use_omp=1
        ;;
    --get-lapack)
        _petsc_extra_args="${_petsc_extra_args} --download-f2cblaslapack"
        ;;
    --get-boost)
        _petsc_extra_args="${_petsc_extra_args} --download-boost"
        ;;
    --get-hdf5)
        _petsc_extra_args="${_petsc_extra_args} --download-hdf5"
        ;;
    *)
        echo "Unsupported argument: $1"
        exit 1
        ;;
    esac
    shift
done

_petsc_extra_args="${_petsc_extra_args} ${_petsc_dir_arch_set}"

# Check for CUDA
if ! [ -x "$(command -v nvidia-smi)" ]; then
    _use_cuda=0
fi

if [ ${_use_cuda} -eq 1 ]; then
    _cuda_sm_ver=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader)
    _cuda_sm=${_cuda_sm_ver//./}
    _petsc_cuda_args=" \
        --with-cuda-arch=${_cuda_sm} \
        --CUDAOPTFLAGS=-O3"
    _petsc_extra_args="${_petsc_extra_args} ${_petsc_cuda_args}"
fi

cd ${_petsc_dir}
echo "Working directory:"
echo "$PWD"

if [ ${_do_cleanup} -eq 1 ]; then
    _clean_cmd="git clean -xdf"
    _reset_cmd="git reset --hard"
    if [ ${_dry_run} -eq 0 ]; then
        ${_clean_cmd}
        ${_reset_cmd}
    else
        echo "Cleanup:"
        echo ${_clean_cmd}
        echo ${_reset_cmd}
    fi
fi

if [ ${_do_pull} -eq 1 ]; then
    _pull_cmd="git pull"
    if [ ${_dry_run} -eq 0 ]; then
        ${_pull_cmd}
    else
        echo "Pull:"
        echo ${_pull_cmd}
    fi
fi

_conf_cmd="./configure \
    ${_petsc_dir_arch_set} \
    --prefix=${_prefix} \
    --with-cc=mpicc \
    --with-cxx=mpicxx \
    --with-fc=0 \
    --with-cuda=${_use_cuda} \
    --with-openmp=${_use_omp} \
    --with-debugging=${_debug} \
    --with-shared-libraries \
    --with-64-bit-indices \
    --download-kokkos \
    --download-kokkos-kernels \
    --COPTFLAGS=${_opt_flags} \
    --CXXOPTFLAGS=${_opt_flags}"

_conf_cmd="${_conf_cmd} ${_petsc_extra_args}"
_build_cmd="make ${_petsc_dir_arch_set} all"
_install_cmd="make ${_petsc_dir_arch_set} install"

if [ ${_dry_run} -eq 0 ]; then
    ${_conf_cmd}
    ${_build_cmd}
    ${_install_cmd}
else
    echo "Configure:"
    echo ${_conf_cmd}
    echo ${_build_cmd}
    echo ${_install_cmd}
fi

