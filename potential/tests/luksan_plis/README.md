# LuksanPLIS: the L-BFGS of NLopt with a ring-buffer history

`potential/LuksanPLIS.cpp` is a modified copy of Luksan's PLIS as shipped in NLopt (`src/algs/luksan`), which is
what `NLOPT_LD_LBFGS` runs. Luksan's routines are under the GNU LGPL; the file keeps that notice.

- `generate.py OUTPUT.cpp [NLOPT_SRC]` regenerates the file from an NLopt checkout (v2.10.0-17-gbb0d24b was used):
  the routines are copied verbatim and made `static`, and the edits, marked `LOCAL`, `RING` and `FAST` in the
  output, are applied as exact text replacements that stop the script if NLopt's text has changed.
- `bash run_tests.sh` builds `test_luksan_plis.cpp` against the installed NLopt and checks 96 cases (two
  objectives, 2 to 30,000 variables, history 1/3/40, evaluation limits and stop values): strict mode must give
  NLopt's evaluated points, evaluation count, return value, minimizer and minimum bit for bit; fast mode must
  repeat bitwise on one thread and on a team of four and reach the same minimum when both converge.

Strict mode matches the NLopt version the file was generated from. Another NLopt version may differ in its own
PLIS, in which case the bitwise comparison of this test, not the optimizer, is what needs regenerating.
