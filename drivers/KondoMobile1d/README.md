# Mobile Kondo-chain driver

`dmrg_kondo_mobile_1d` uses the homogeneous twelve-state Kondo site.
Each physical position contains one conduction orbital and one projected
impurity orbital with local states `empty`, `up`, and `down`.

The driver fixes the global impurity population through `NIMP`.  During
the symmetric iDMRG warm-up it uses

```text
Nimp(L) = ceiling(L*NIMP/Lfinal)
Nup(L) + Ndown(L) = L + Nimp(L)
```

where the first term on the right is a half-filled conduction band.  The
spin sector has the smallest non-negative `2*Sz`.  The `ceiling`
ensures that impurities occur before the last growth step, so their
states can enter the renormalized blocks.  The explicit schedule also
avoids independently rounding the three conserved charges, which could
otherwise request an inconsistent intermediate sector.

The conduction band is currently fixed at half filling, `Nc(L)=L`.
A future extension can replace this assumption with a spin-resolved
occupation derived from `QN_DENSITY(1:2)`.

`TS` and `TIMP` are the uniform hopping amplitudes of the conduction and
projected impurity orbitals.  The present driver requires PBC iDMRG.

Positive `PINNING` lowers the impurity energy by that amount on the
distinct physical positions in `PINNING_SITES`.  This selects a spatial
configuration in the disconnected `TIMP=0` limit.  In this case the
target schedule counts the pinned sites already present at each growth
length, so a static impurity enters together with its physical site.
With `IMEASURE=T`, the driver writes the local impurity density.
`ICORRELATIONS=T` also writes the full impurity-electron and
impurity-impurity spin-correlation matrices, both as longitudinal and
scalar products.

## ED reference

The `ed_reference` directory compares the six-site mobile problem with
`Kondo_ED`.  It uses one delocalized impurity, `TS=TIMP=-1`, PBC, and
isotropic `J=1`.  The DMRG benchmark keeps all `12^2=144` states of the
first enlarged block.  This is necessary because sectors with zero
weight at one warm-up length can contribute after the target changes.

The regression checks the total energy, the complete local impurity
density profile, and all longitudinal impurity-bath correlations.  Run
it after loading the `Kondo_ED` module:

```sh
module load kondo_ed/gnu/0.1.2
DMRG_EXE=$HOME/.bin/dmrg_kondo_mobile_1d \
ED_EXE=$HOME/.bin/1d_Kondo_chain \
./ed_reference/compare_ed_dmrg.sh
```
