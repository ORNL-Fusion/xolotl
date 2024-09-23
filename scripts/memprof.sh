#!/bin/bash

set -eu -o pipefail
readonly _self_path=$(cd $(dirname "${BASH_SOURCE[0]}"); pwd)

_param_file=""

while [ $# -gt 0 ]
do
    case $1 in
    --kp-root=*)
        _kp_root="${1:10}"
        export KOKKOS_TOOLS_LIBS="${_kp_root}/lib/libkp_memory_events.so"
        ;;
    *)
        _param_file=$1
        ;;
    esac
    shift
done

export PETSC_OPTIONS="-use_gpu_aware_mpi 0"

time ./xolotl/xolotl ${_param_file} &
_time_pid=$!
# wait ${_pid}
_pid=$(pidof -s xolotl)
wait ${_time_pid}

_events_file="$(hostname)-${_pid}.mem_events"
echo "events file: $_events_file"

julia ${_self_path}/../analysis/memprof.jl ${_events_file}
