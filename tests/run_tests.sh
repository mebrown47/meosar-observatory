#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Runs the unit tests in BUILD_DIR (default build/) and prints one RESULT line.
set -u
B="${1:-$(dirname "$0")/../build}"
pass=0; fail=0
for t in test_fgb test_sgb test_repeat; do
  echo "== $t"
  if "$B/$t"; then pass=$((pass+1)); else fail=$((fail+1)); echo "FAILED: $t"; fi
done
echo "RESULT: $pass passed, $fail failed"
[ "$fail" -eq 0 ]
