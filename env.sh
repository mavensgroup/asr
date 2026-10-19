# source this: Intel oneAPI + fpm flags for ASR
source /opt/intel/oneapi/setvars.sh >/dev/null 2>&1
export FPM_FC=ifx
export ASR_FFLAGS="-O2 -g -traceback -qopenmp -qmkl=parallel -free -standard-semantics -warn all -fp-model=precise -heap-arrays 1024"
export ASR_LDFLAGS="-qopenmp -qmkl=parallel"
asr_build() { fpm build --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_test()  { fpm test  --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
asr_run()   { fpm run   --compiler ifx --flag "$ASR_FFLAGS" --link-flag "$ASR_LDFLAGS" "$@"; }
# Code lives on GitHub: https://github.com/mavensgroup/asr (public).  The laptop pushes over ssh (remote "origin"); the desktop
# cannot reach GitHub over ssh (port 22 reset on the campus network) and pulls read-only over https (remote "github") into its
# working copy ~/asr_remote4.  The old bare hub desktop:~/asr.git is retired.
#   asr_sync : push local main to GitHub and fast-forward the desktop working copy (run on the laptop)
#   asr_pull : fetch the latest main (run on either machine)
asr_sync() { git push origin main && ssh desktop 'cd ~/asr_remote4 && git pull --ff-only -q github main && git log --oneline -1'; }
asr_pull() { if git remote | grep -qx github; then git pull --ff-only github main; else git pull --ff-only origin main; fi; }
