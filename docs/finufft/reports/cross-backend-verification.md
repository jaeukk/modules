## Verification verdict

Diagnostic replay evidence: **PASS**. Endpoint agreement objective: **unresolved / FAIL**.

- Executed and repeated (identical captured stdout) the requested replay commands for seed 20260928:
  - `direct` and `finufft12`: indices 1, 57, 149, 800
  - later-pass `direct`: 1038, 1116, 1117, 1141
  - later-pass `finufft12`: 1079, 1124, 1125, 1150  
  All exited 0, passed gates, used actual expected species backends/statuses (FINUFFT status 0; direct no FINUFFT), and reproduced the matching saved backend energy/gradient exactly 0.

- Inspected all 17 non-replay optimization receipts: all exit 0, ground states, and satisfy maxima: energy 1.384e-19, maxS 8.448e-19, force 2.197e-12, overlap 0. FINUFFT receipts have `finufft_calls == rho_calls > 0`, one FINUFFT thread, and zero nonzero-status calls; direct receipts have zero FINUFFT calls. `history1_finufft` is correctly recorded as statuses 5, 1, 1; all other passes are status 1.

- Inputs hash-identical within each seed. All six baseline final files exactly reproduce the earlier validated endpoint files.

- Inspected all 12 replay2 receipts: 90 snapshots across passes 0–2, all successful with no fallback and exact matching saved errors. Aggregate maxima: rho 7.0089e-11, force 3.4564e-14, energy 4.2633e-13.

- Translation-aligned endpoint residuals remain large: baseline max fractional norm 0.04645; epsilon/history controls max 0.06151. Mean wrapped translation is negligible in Cartesian units (baseline ≤2.69e-13; controls ≤5.96e-13). This does not meet the provisional 1e-7 endpoint gate.

- Verified all 26 manifest/replay-manifest source, binary, and reused-object hashes. `trace.cpp` is the history version; the baseline source snapshot is recorded. Hashes establish input consistency, not reproducible compilation. Production code was unchanged and existing dirty fixes were preserved.