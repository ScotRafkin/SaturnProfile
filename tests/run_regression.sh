#!/usr/bin/env bash
# The regression: runs every accepted suite, or the ones named as arguments, and checks each against its reference count.
# Run from the repository root: `bash tests/run_regression.sh [suite ...]`, a suite named as `step04_5` or `step1/verify_review_changes`.
set -u
export PYTHONIOENCODING=utf-8

# Each accepted suite and its reference count. A count that changes at an acceptance changes here in
# the same commit. SPEC_03 and SPEC_04 counts are STATE.md's; STATE.md carries none for SPEC_01 and
# SPEC_02, so theirs are the full regression STATE.md cites for SPEC_04 Step 5 (25 of 25).
# step03_0 is absent on purpose: it checks its one-time rebuild against `before/` products that
# today's reader refuses (no `epoch`, SPEC_00 v0.17); every later suite covers the same products.
SUITES="
step1/accept_step1 6
step1/verify_review_changes 14
step2/accept_step2 8
step3/accept_step3 5
step4/accept_step4 7
step5/accept_step5 7
step6/accept_step6 6
step7/accept_step7 8
step8/accept_step8 9
step9/accept_step9 6
step02_1/accept_step02_1 6
step02_2/accept_step02_2 7
step02_3/accept_step02_3 9
step02_4/accept_step02_4 9
step02_5/accept_step02_5 7
step02_6/accept_step02_6 7
step03_1/accept_step03_1 9
step03_2/accept_step03_2 8
step03_3/accept_step03_3 16
step03_4/accept_step03_4 13
step04_0/accept_step04_0 15
step04_1/accept_step04_1 13
step04_2/accept_step04_2 13
step04_3/accept_step04_3 9
step04_4/accept_step04_4 11
step04_5/accept_step04_5 11
step05_1/accept_step05_1 11
step05_2/accept_step05_2 3
step05_3/accept_step05_3 7
"

# The suites to run: all of them, or those named (a directory name selects every suite in it).
selected=""
if [ $# -eq 0 ]; then
  selected=$(echo "$SUITES" | awk 'NF == 2')
else
  for name in "$@"; do
    rows=$(echo "$SUITES" | awk -v n="$name" 'NF == 2 && ($1 == n || $1 ~ "^" n "/")')
    if [ -z "$rows" ]; then echo "no suite named $name" >&2; exit 2; fi
    selected="$selected
$rows"
  done
  selected=$(echo "$selected" | awk 'NF == 2 && !seen[$1]++')
fi

D=occul_data/lindal
F=forward/lindal_closure
T=forward/lindal_transfer
R=reports/regression
B=$R/aside
mkdir -p "$R" || exit 1
out=$R/regression.txt; : > "$out"

rm -rf "$B"; mkdir -p "$B" || exit 1
cp "$D/lindal_refractivity.nc" "$B/" && cp -r "$D/figures" "$B/figures" || exit 1
mkdir -p "$B/closure" && cp -r "$F/inputs" "$F/output" "$B/closure/" || exit 1
mkdir -p "$B/transfer" && cp -r "$T/inputs" "$B/transfer/" || exit 1
cp "$T/output/lindal_transfer_profile.nc" "$B/transfer/" || exit 1
before=$(python -c "from casspian.lib import io; print(io.sha256('$B/lindal_refractivity.nc'))")
fbefore=$(cd "$F" && sha256sum inputs/*.nc output/lindal_closure_profile.nc | sha256sum)
tbefore=$(cd "$T" && sha256sum inputs/*.nc output/lindal_transfer_profile.nc | sha256sum)
dirty() { git status --porcelain | sed 's/^/        /'; }
echo "the tree carries $(git status --porcelain | wc -l) paths of its own before the run" >> "$out"
dirty >> "$out"

restore() {
  cp "$B/lindal_refractivity.nc" "$D/lindal_refractivity.nc" || exit 1
  rm -rf "$D/figures" && cp -r "$B/figures" "$D/figures" || exit 1
  rm -rf "$F/inputs" "$F/output" && cp -r "$B/closure/inputs" "$B/closure/output" "$F/" || exit 1
  rm -rf "$T/inputs" && cp -r "$B/transfer/inputs" "$T/" || exit 1
  cp "$B/transfer/lindal_transfer_profile.nc" "$T/output/lindal_transfer_profile.nc" || exit 1
  git checkout -- reports/figures/
}

k=0; n=0
while read -r s ref; do
  [ -z "$s" ] && continue
  n=$((n + 1))
  log="$R/regression_$(basename "$s").txt"
  start=$SECONDS
  python "tests/$s.py" > "$log" 2>&1
  took=$((SECONDS - start))
  line=$(grep -E '[0-9]+ of [0-9]+ checks pass' "$log" | tail -1)
  passed=$(echo "$line" | sed -nE 's/.*\b([0-9]+) of ([0-9]+) checks pass.*/\1/p')
  total=$(echo "$line" | sed -nE 's/.*\b([0-9]+) of ([0-9]+) checks pass.*/\2/p')
  if [ -n "$passed" ] && [ "$passed" = "$total" ] && [ "$total" = "$ref" ]; then
    k=$((k + 1)); mark="        "
  else
    mark="MISMATCH"
  fi
  echo "$mark ${passed:-?} of ${total:-?} checks pass, reference $ref, ${took} s   <- $s" >> "$out"
  restore
  echo "    paths dirty after $s: $(git status --porcelain | wc -l)" >> "$out"
  dirty >> "$out"
done <<< "$selected"

after=$(python -c "from casspian.lib import io; print(io.sha256('$D/lindal_refractivity.nc'))")
fafter=$(cd "$F" && sha256sum inputs/*.nc output/lindal_closure_profile.nc | sha256sum)
tafter=$(cd "$T" && sha256sum inputs/*.nc output/lindal_transfer_profile.nc | sha256sum)
echo "registered kind N restored: sha256 before $before, after $after, equal $([ "$before" = "$after" ] && echo yes || echo NO)" >> "$out"
echo "closure inputs and product restored: equal $([ "$fbefore" = "$fafter" ] && echo yes || echo NO)" >> "$out"
echo "transfer inputs and product restored: equal $([ "$tbefore" = "$tafter" ] && echo yes || echo NO)" >> "$out"
echo "porcelain at the end:" >> "$out"
git status --porcelain >> "$out"
echo "$k of $n suites at their reference counts" >> "$out"
cat "$out"
[ "$k" -eq "$n" ]
