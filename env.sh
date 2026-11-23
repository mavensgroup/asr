# source this: Intel oneAPI + fpm flags for ASR
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
export FPM_FC=ifx
export ASR_FFLAGS="-O2 -g -traceback -qopenmp -qmkl=parallel -free -standard-semantics -warn all -fp-model=precise -heap-arrays 1024"
export ASR_LDFLAGS="-qopenmp -qmkl=parallel"
asr_build() { [ -x ./tools/gen_build_info.sh ] && ./tools/gen_build_info.sh; fpm build --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_test()  { [ -x ./tools/gen_build_info.sh ] && ./tools/gen_build_info.sh; fpm test  --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_run()   { [ -x ./tools/gen_build_info.sh ] && ./tools/gen_build_info.sh; fpm run   --target asr --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" -- "$@"; }
# Optional local additions (machine-specific helpers), not part of the distribution.
[ -f "$(dirname "${BASH_SOURCE[0]}")/env.local.sh" ] && source "$(dirname "${BASH_SOURCE[0]}")/env.local.sh"
# asr_test      : the fast unit tests (fpm test, test/), about 2 min
# asr_test_slow : the slow self-consistent and enumeration tests (test_slow/, fpm examples); prints the failed ones, nonzero exit if any
asr_test_slow() {
  local f n rc=0 failed=""
  for f in test_slow/test_*.f90; do
    n=$(basename "$f" .f90)
    fpm run --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" --example "$n" || { rc=1; failed="$failed $n"; }
  done
  [ $rc -eq 0 ] && echo "asr_test_slow: all passed" || echo "asr_test_slow: FAILED:$failed"
  return $rc
}
