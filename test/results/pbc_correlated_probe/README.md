# BHZ PBC correlated ground-state prediction

The BHZ driver was compiled against the complex MPI nssDMRG library and run
with four local ranks. Inputs and full logs for paired runs are included.
`LDMRG=10` is promoted to 11 by the PBC driver. Both pairs request one
eigenstate and use `WF_PREDICT_DIAG=T` to print overlap, residual, and
ARPACK's `NumOp`; the diagnostic adds one H*v outside `NumOp` per eligible
step. `LANC_GS_PREDICTION=T` is present only in the `seed` inputs.

| MDMRG | ARPACK H*v, plain | ARPACK H*v, seeded | Eligible-step H*v, plain -> seeded | Maximum energy/site difference |
|---:|---:|---:|---:|---:|
| 20 | 540 | 380 | 425 -> 265 | 3.2e-14 |
| 80 | 370 | 290 | 245 -> 165 | 1.8e-14 |

The first two steps have no correlated seed. At `M=20` the eligible-step
ground-state weight is about 0.97-0.975. At `M=80` it ranges from 0.526 to
0.973. The streamed projection of the previous superblock coefficients into
the two half tensors took about 0.02 CPU seconds per step at `M=80`; this is
separate from the ARPACK time. `M=500` memory and time have not been measured.

Run an input from its own empty directory as `DMRG.conf`; the driver writes
additional DMRG output files to the current directory. For example:

```sh
mpirun -np 4 /path/to/dmrg_bhz_1d
```

The predictor is enabled by `LANC_GS_PREDICTION=T` and requires
`LANC_NEIGEN=1`, `DMRGTYPE=i`. It does not require `WF_PREDICT_DIAG=T` in
production.
