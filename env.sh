# source this: Intel oneAPI + fpm flags for ASR
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
export FPM_FC=ifx
export ASR_FFLAGS="-O2 -g -traceback -qopenmp -qmkl=parallel -free -standard-semantics -warn all -fp-model=precise -heap-arrays 1024"
export ASR_LDFLAGS="-qopenmp -qmkl=parallel"
asr_build() { fpm build --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_test()  { fpm test  --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_run()   { fpm run   --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
