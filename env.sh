# source this: Intel oneAPI + fpm flags for ASR
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
export FPM_FC=ifx
export ASR_FFLAGS="-O2 -g -traceback -qopenmp -qmkl=parallel -free -standard-semantics -warn all -fp-model=precise -heap-arrays 1024"
export ASR_LDFLAGS="-qopenmp -qmkl=parallel"
asr_build() { fpm build --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_test()  { fpm test  --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_run()   { fpm run   --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
# Code sync between this machine and the desktop goes through the bare hub repository desktop:~/asr.git (remote "hub";
# not "origin": Claude Remote Control has had trouble with origin URLs).  The desktop working copy is ~/asr_remote4.
#   asr_sync : push local main to the hub and fast-forward the desktop working copy (run on the laptop)
#   asr_pull : fetch commits made on the other side (run on either machine)
asr_sync() { git push hub main && ssh desktop 'cd ~/asr_remote4 && git pull --ff-only -q hub main && git log --oneline -1'; }
asr_pull() { git pull --ff-only hub main; }
