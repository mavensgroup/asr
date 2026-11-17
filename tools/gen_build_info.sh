#!/bin/bash
# Generate src/asr_build_info.f90 from git state
set -e
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# cd instead of git -C: git 1.8 (CentOS 7 clusters) has no -C
REV=$(cd "$ROOT" && git describe --always --dirty 2>/dev/null || echo "unknown")
DIRTY="false"
if [[ "$REV" =~ -dirty$ ]]; then
  DIRTY="true"
  REV="${REV%-dirty}"
fi
TARGET="$ROOT/src/asr_build_info.f90"
TMP=$(mktemp)
cat <<EOF > "$TMP"
!> Build provenance generated at compile time by tools/gen_build_info.sh
module asr_build_info
  implicit none
  character(len=*), parameter :: git_rev = '$REV'
  logical, parameter :: git_dirty = .$DIRTY.
end module asr_build_info
EOF

if [ ! -f "$TARGET" ] || ! cmp -s "$TMP" "$TARGET"; then
  mv "$TMP" "$TARGET"
else
  rm -f "$TMP"
fi

if (cd "$ROOT" && git ls-files --error-unmatch "$TARGET") >/dev/null 2>&1; then
  (cd "$ROOT" && git update-index --skip-worktree "$TARGET") 2>/dev/null || true
fi
