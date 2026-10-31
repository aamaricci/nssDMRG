# ARPACK tolerance scan (2026-10-31)

BHZ OBC, LDMRG=6 (12 sites), MDMRG=20 with degeneracy extension,
MH=2, lambda=0.3, QN_DENSITY=(1,1), one MPI rank, prediction disabled.
Two interactions: (U,Ust,Jh)=(1,0.5,0.25) and (8,4,2).
Default scan: NEV=2, NCV=20, GS_DEG_THRESHOLD=1e-10.
Each run starts fresh; these are full DMRG trajectories, not a frozen-H
comparison. Timing is one short measurement per case. summary.csv lists
solver NumOp (excluding the two extra residual checks per step), wall time,
maximum measured absolute residual, energy/site differences across all
steps/states, final gap, and selected GS count by step.

U1: tolerance 1e-12/1e-8/1e-6/1e-4 costs 910/633/513/343 products.
The largest energy/site changes are 0/1.8e-14/6.1e-11/4.0e-4.
U8: the two-state loose solves fail to retain the low doublet at late steps.
NCV=40 at tolerance 1e-8 does not fix this. NEV=3, NCV=30, tolerance 1e-8
recovers the doublet with 1092 products versus 1610 in the NEV=2 strict
reference, and energies/site agree within 2.4e-13 for the first two states.

Additional cluster experiment: NEV=3, NCV=30, fixed GS cutoff 1e-5.
At tolerance 1e-12: 2148 products, 10.09 s, GS counts 1/2/2/2/2.
At tolerance 1e-8: 1044 products, 6.07 s, GS counts 1/2/2/1/1.
Maximum residual in the loose run is 5.57e-7, below cutoff/10.
Nevertheless the trajectories differ: final GS energy/site differs by
2.8e-5; the loose final splitting is about 4.4e-4. Thus checking only
residual versus cluster cutoff does not validate the multiplet or the
resulting DMRG truncation. No automatic cutoff/tolerance coupling has been
implemented. State overlaps/projector stability and fixed-H comparisons
remain necessary before adopting such a policy in the near-degenerate case.

LANC_CHECK_RESIDUAL=T adds one H*v per returned state and logs absolute
||H psi_i-E_i psi_i|| after the solve. Default false. Real/complex builds
compile; the numerical runs here use complex arithmetic and one MPI rank.

Intermediate tolerance 1e-10, U8, NEV=2, NCV=20, no prediction:
1441 products versus 1610 at 1e-12; wall 7.21 versus 7.72 s.
The GS counts remain 1/1/1/1/2 and both energies/site agree with the
reference within 1.3e-13 across all steps. Maximum absolute residual
is 5.4e-9; the final splitting is approximately 2.1e-11. Thus this
particular run reproduces the doublet but the residual by itself does
not resolve such a small splitting. Inputs and outputs: u8_t10/.
