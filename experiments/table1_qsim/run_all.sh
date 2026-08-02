#!/bin/sh
# Full validation sequence. Grid first; timing/transient only if the grid exits 0
# (a failed r=1 gate means nothing downstream is meaningful -- spec 7.3).
cd "$(dirname "$0")"
python3 validate_table1_qsim.py > run_full.log 2>&1
grid=$?
echo "GRID_EXIT=$grid" >> run_full.log
if [ "$grid" -eq 0 ]; then
  python3 validate_table1_qsim.py --time-baseline --transient-check > run_timing.log 2>&1
  echo "TIMING_EXIT=$?" >> run_timing.log
else
  echo "SKIPPED: grid exited $grid, so the timing/transient passes did not run." > run_timing.log
fi
echo "ALL_DONE" >> run_full.log
