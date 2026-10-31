# Small BHZ PBC ground-state prediction probe

Single process, PBC, effective LDMRG=5 (10 sites), MDMRG=20,
LANC_NEIGEN=1. WF_PREDICT_DIAG=T. The baseline candidate carries the
previous renormalized center and a uniform two-site charge-sector pair,
with explicit n/p tensor order and fermion parity signs. It does not
transport the optimized site-environment tensors across alternating edges.

For the three steps with an available trial (links p,n,p), GS projection
weights are 0.026505, 0.026302, 0.026571. A single Hamiltonian action
provides the tangent direction (H-E_trial) v. Normalized filtered trials
v - tau (H-E_trial)v at tau=0.1 yield 0.074092, 0.073639,
0.074246. The improvement is too small to justify supplying this as an
ARPACK seed. The diagnostic uses the converged state to measure overlaps;
this state is never used to construct the trial.

The probe is intentionally separate from LANC_GS_PREDICTION, which still
rejects PBC. It shows the need for a genuine inter-edge tensor transport.
