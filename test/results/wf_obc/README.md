# BHZ OBC prediction experiment, 2026-10-31

Driver: DriverDatabase/nssDMRG/dmrg_bhz_1d.f90, linked against the local
complex MPI build. Executed directly with one MPI rank. Each case is a
fresh process in an empty directory, with its supplied configuration named
DMRG.conf. Timings are individual short runs, not statistical estimates.

Parameters: LDMRG=6 (12 sites), MDMRG=20 (the degeneracy extension can retain
more states), MH=2, LAMBDA=0.3, ULOC=(1,1), UST=0.5, JH=0.25,
QN_DENSITY=(1,1), lazy H*v, tolerance 1e-12, NCV=15.

The first two growth steps do not have a correlated seed. For 8/10/12 sites:

| eigenstates | no seed NumOp | correlated seed NumOp |
|---|---|---|
| 2 | 199 / 199 / 211 | 208 / 208 / 220 |
| 1 | 55 / 55 / 63 | 31 / 23 / 23 |

Two eigenstates: wall 5.49 s without prediction, 5.84 s with prediction.
One eigenstate: wall 2.90 s control, 2.43 s seeded. The one-state control
also constructs/evaluates the candidate to isolate the seed effect.
Diagnostic H*v calls (one per available trial) are additional to NumOp.
Correlated-contraction CPU time is 0.0025--0.0036 s per eligible step;
this excludes the center projection, which is included in total run time.
Maximum energy/site differences: 3.74e-14 (two states), 9.77e-15 (one state).

The ground state is not degenerate in this example. The result does not
establish convergence of an exactly degenerate manifold, MPI scaling, or
performance in the strongly entangled BHZ regime. Only the first old
state is predicted. A single ARPACK starting vector cannot directly carry
a whole predicted subspace.

The historical benchmark used `WF_PREDICT_SEED=T`; the current input name is
`LANC_GS_PREDICTION=T`, restricted to `LANC_NEIGEN=1`. The initial factorized
trial is never supplied. Prediction and diagnostics default to false.
The prototype uses dense temporary contractions and a relative SVD
pseudoinverse cutoff 1e-10; it is intended for small experiments.
