# DMRG regression data

Each model directory contains a manual reference driver, its model-specific
`DMRG.conf.*` input, and numerical tables produced by a trusted, independent
run.  These files are not consumed directly by CMake.

After validating a reference run, manually copy its `DMRG.conf.*` file and all
required `*.check` tables to the corresponding directory under `test/src`.
CMake stages those approved files in `build/test/work/<model>`, where the CTest
driver starts a complete DMRG run from fresh blocks.

Blank lines and lines beginning with `#` or `!` are ignored.  All remaining
lines must contain only the documented numeric columns.

## Spin1d_SU2

- `energy.check`: the test compares `left_length` and the ground-state energy
  per site; any additional computed eigenvalues are retained in the table but
  are not part of the regression criterion
- `entropy.check`: `left_length entropy_left entropy_right right_length`
- `spin_local.check`: `i <Sz_i> <Sz_i^2>` for `i=1,...,2*Ldmrg`
- `spin_nn.check`: `i i+1 <S_i.S_(i+1)>`
- `spin_1j.check`: `1 j <S_1.S_j>`

## Hubbard1d

- `energy.check`: the test compares `left_length` and the ground-state energy
  per site; any additional computed eigenvalues are retained in the table but
  are not part of the regression criterion
- `entropy.check`: `left_length entropy_left entropy_right right_length`
- `hubbard_local.check`: `i <n_up> <n_down> <n_up*n_down>`
- `density_nn.check`: `i i+1 C(1,1) C(1,2) C(2,1) C(2,2)`
- `density_1j.check`: `1 j C(1,1) C(1,2) C(2,1) C(2,2)`

For the density tables the general ordering is `io` outermost and `jo`
innermost, with `io,jo=1,...,Nspin*Norb`.

The default comparison criterion is

`abs(value-reference) <= 1e-8 + 1e-7*abs(reference)`.

Observable and correlation tables use an absolute tolerance of `1e-6` to
account for the run-to-run variation of truncated eigenvectors.  In
particular, the reference local spin magnetization is zero and its residual
value depends on the numerical eigenvector selected inside a nearly degenerate
subspace.
