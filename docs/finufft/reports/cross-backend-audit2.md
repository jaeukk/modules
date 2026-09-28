**Diagnostic verdict: PASS for the bounded sampled-state conclusion. Endpoint objective: FAILED, unchanged.**

- **Frozen reconstruction:** Replay follows production `GetConfig` Switch 0: preserve the seeded basis, remove particles, insert coordinates divided by the pass scale, and retain particle ordering. Shared setup restores species assignments. The optimizer-gradient conversion uses the appropriate negative force projection divided by scale.

- **Caches and lifetimes:** Each snapshot/backend gets a fresh potential. The configuration and modes remain alive throughout evaluation; returned forces and rho values own their storage. No replay-level stale-cache or dangling-reference defect found. Exact saved-value reproduction also supports consistency with production callbacks at the sampled states.

- **Independent reference:** The reference computes trigonometric sums directly and contact interactions through all particle pairs with minimum-image distances. It does not reuse backend rho or production neighbor enumeration. It shares configuration, modes, and species conventions, so independence is computational rather than a separate specification.

- **Backend evidence:** Replay checks FINUFFT tolerance `1e-12`, sets one FINUFFT thread, and requires actual FINUFFT backend plus status zero for every species. Direct evaluation requires direct backend and `FinufftNotCalled`. Both inspected receipts satisfy these checks.

- **Results and provenance:** The six-entry summary reports 66 passing snapshots and zero matching-backend saved energy/gradient errors. The two inspected `replay2` receipts independently support their corresponding entries, identify trace paths/seeds, and separate stderr from parsed results. The manifest records build commands and source/binary hashes; dependency identity relies on the referenced shared manifest and was not independently authenticated within this review.

**Limits:** Both inspected receipts sample only pass 0; these results do not establish later-pass or endpoint correctness. Saved-value agreement is reported separately, not included in `gates_passed`; the explicit zero-error summary matters. Zero differences establish numeric equality, not bitwise identity. Fresh replay does not exhaustively test history-dependent cache behavior, although matching saved outputs provide evidence against such effects at these callbacks.

**Actionable blockers: None for the proposed conclusion.** The replay sufficiently closes the frozen-coordinate diagnostic gap: numerical trajectory sensitivity is strongly supported, with no state-dependent mapping/force discrepancy found at sampled states. It does not establish exclusive causation, accepted-step correspondence, chaos, or universal bug exclusion. No production change is justified by this diagnostic alone.

Read-only review completed; no optimizations or replay executions rerun.