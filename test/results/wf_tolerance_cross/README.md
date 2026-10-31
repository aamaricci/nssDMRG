# Controlled prediction/tolerance comparison, 2026-10-31

BHZ OBC, 12 sites (LDMRG=6), MDMRG=20 with degeneracy extension,
MH=2, lambda=0.3, (U,Ust,Jh)=(1,0.5,0.25), half filling.
NCV=10*NEV, no additive offset. One process; fresh initial state per run.
All runs check Ritz residuals (NEV extra H*v per step, excluded from NumOp).
Prediction runs additionally diagnose the candidate. Timings include these
costs and construction; they are short individual measurements.

summary.csv records all solver products, products and solver time on the
three eligible prediction steps (8/10/12 sites), total elapsed time, and
maximum energy/site difference across states and steps against the same
NEV strict unseeded run.

IMPORTANT: n2_t8_seed returns a wrong second-lowest level compared with the
strict and unseeded controls, despite small residuals. Its apparent speedup
is NOT accepted. At 8 sites the second energy/site is -4.500785860435260,
versus -4.536921849154052 in the reference. This is consistent with a
single GS seed having insufficient weight in other symmetry sectors; the
exact sector mechanism has not been independently established.

Accepted results: NEV=1 seed reduces eligible-step products 185->75 at
1e-12 and 125->30 at 1e-8. NEV=2 at 1e-12 increases them 577->607.
Without prediction, relaxing tolerance reduces full-run products 910->633
for NEV=2 while energies/site agree to 1.8e-14 on this nondegenerate case.
These statements do not establish safety for near-degenerate multiplets;
see ../lanc_tolerance for the contrasting U8 experiment.

## Intermediate tolerance 1e-10

Same inputs and binary; all four new cases completed successfully.
NEV=1: total NumOp 305->255 on relaxing tolerance without seed;
with seed at 1e-10 the total is 165. Eligible-step counts are 150->60
without/with seed, with energies/site differing from the strict unseeded
reference by at most 1.9e-14.
NEV=2: unseeded total 910->756, energy/site agreement 1.3e-14.
Seeded total 716 is NOT an accepted speedup: at 8 sites the second state
is wrong (energy/site discrepancy 0.03614), although at 10 and 12 sites
it recovers the reference energies. The maximum residual is only 1.9e-9.
