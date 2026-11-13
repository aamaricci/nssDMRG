# Fixed-impurity Kondo driver

`dmrg_kondo_1d` implements the fixed-spin Kondo problem using one
eight-state electron-spin site and ordinary four-state electron sites.  The
current driver supports the PBC infinite-DMRG construction and places the
impurity on the left seed selected by the growth geometry.

## ED reference

The `ed_reference` directory contains a six-bath-site comparison against
`Kondo_ED`.  Both calculations use a half-filled bath, one spin-1/2 impurity,
periodic boundary conditions, `t=-1`, and isotropic antiferromagnetic exchange
`J=1`.

The two codes use different exchange conventions.  The DMRG coupling

```text
Jz Simp_z sele_z + Jxy/2 (Simp_+ sele_- + Simp_- sele_+)
```

is reproduced in `Kondo_ED` by

```text
JK_Z  = -Jz/2
JK_XY = -Jxy/2
```

because the ED Hamiltonian contains `-2 JK_Z Simp_z sele_z` and the analogous
transverse term.  `ED_FILLING=7` counts the six bath electrons and the
auxiliary fermion used to represent the occupied impurity spin.  The direct
correlation comparison uses `<Simp_z sele_z(i)>`, read from the last row of
`Sz_corr.ed`; this avoids assuming exact spin-rotation symmetry in the
truncated DMRG state.

Load the environment containing `Kondo_ED`, compile `dmrg_kondo_1d`, and run

```sh
module load kondo_ed/gnu/0.1.2
DMRG_EXE=$HOME/.bin/dmrg_kondo_1d \
ED_EXE=$HOME/.bin/1d_Kondo_chain \
./ed_reference/compare_ed_dmrg.sh
```

The executable paths, `MPIEXEC`, `NPROC`, `WORK_ROOT`, `ENERGY_TOLERANCE`, and
`CORR_TOLERANCE` can be overridden through the environment.  The reference
tolerances are `1e-10` for the total ground-state energy and `1e-9` for the
maximum site-resolved correlation error.  `MDMRG=32` preserves the complete
intermediate Kondo-block basis, so the final six-site superblock spans the
same symmetry sector as ED.  The residual correlation difference is limited
by the output precision of `Sz_corr.ed`.
