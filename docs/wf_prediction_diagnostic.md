# Diagnostica della prediction della funzione d'onda

Il primo esperimento lascia invariata la diagonalizzazione. Un candidato
`v0` viene costruito prima di `sb_diag`, nel settore `sb_states` del nuovo
superblocco. Si applica una volta `H` a `v0`, poi ARPACK procede senza `v0`.
Dopo la convergenza si registrano norma, energia di Rayleigh,
`||H v0 - E0 v0||`, e proiezione sul sottospazio dei due autostati richiesti.
Il costo di questa diagnostica va tenuto separato da `NumOp`.

## Convenzione degli indici

`tbasis` usa il secondo fattore come indice veloce. Se `a` e `b` sono
indici a partire da uno, l'indice di `(A.o.B)` e' `b+(a-1)*dim(B)`.
Il superblocco usa l'indice destro come veloce:
`sb_states(k)=j+(i-1)*right%dim`. Per `link='n'` le basi ingrandite
sono `(L.o.sL)` e `(sR.o.R)`; per `link='p'` sono `(sL.o.L)` e
`(R.o.sR)`. La differenza fra i due casi deve essere confinata alla
funzione che calcola questi quattro indici.

Il trasporto usa le rotazioni `U_L`, `U_R` effettivamente selezionate da
`renormalize_block` e tutti i coefficienti del primo autostato del passo
precedente. Con la convenzione della RDM destra, il centro rinormalizzato
e' `C=U_L^dagger Psi U_R^*`. Il candidato iniziale e' il prodotto di questo
centro con una coppia locale uniforme tra gli stati che forniscono la
carica aggiunta. Questa scelta e' una baseline diagnostica, non la
prediction di McCulloch: non trasporta le correlazioni della coppia locale.
Per due stati degeneri, il centro del solo primo stato puo' inoltre
dipendere dalla base scelta dal diagonalizzatore.

Il parametro `WF_PREDICT_DIAG` e' falso per default. La contrazione correlata
proietta i soli coefficienti ammessi del superblocco precedente e forma il
candidato una coppia di siti alla volta. Non alloca la matrice densa `Psi`
ne' il tensore predetto a quattro indici. Il costo aritmetico della
proiezione rimane proporzionale a `dim(SB)*M`; va misurato ai grandi `M`.

## Segni fermionici

Ogni scambio di due fattori con parita' `p` e `q` aggiunge
`(-1)^(p*q)`. Le parita' vanno ricavate dai numeri quantici delle
basi locali, non dalla posizione degli indici. Il test minimo usa stati
di parita' pari e dispari per ciascun fattore e verifica la composizione
dei segni per entrambi i link. Il caso a parita' pari da solo non testa
questa regola.

## Criterio per il passo successivo

Il candidato e' utile se la proiezione sul sottospazio dei due stati
bassi e' alta e il residuo diminuisce rispetto a un seme casuale,
senza un costo di costruzione comparabile alle applicazioni di `H` che
si vorrebbero risparmiare. Solo allora conviene passarlo ad ARPACK.

Nel test BHZ PBC a 10 siti (`LDMRG=5`, `MDMRG=20`, due autostati), la
baseline ha peso nel sottospazio basso 0.034, 0.048 e 0.059 ai tre
passi successivi al primo, con residuo circa 3.6. Le energie finali
coincidono con il run senza diagnostica. Il peso e' troppo piccolo per
usare questo candidato come seme; il prossimo esperimento deve usare
le isometrie di due passi per predire anche la coppia locale.

## Correlated OBC/PBC prediction

For either boundary condition the retained half states are contracted as
`F=UL^dagger Psi`, `G=Psi UR*`, and
`Theta(a,s,t,b)=sum F(a,s,beta) Cprev^+(beta,alpha) G(alpha,t,b)`.
`Cprev` refers to the actual old block bases, with no eigenvector matching.
The pseudoinverse uses SVD and a relative cutoff 1e-10. No explicit extra
fermionic permutation is applied to this same-edge contraction; the
ordered tensor coefficients are carried by the half states.

The previous link selects the tensor slices used to form `F` and `G`;
the current link selects the target coefficient positions. This is identical
for OBC but essential for PBC, where growth alternates between adding sites
at the inner and outer edges. In detail, `n` uses the previous grouping
`(oldL,siteL)/(siteR,oldR)` while `p` uses `(siteL,oldL)/(oldR,siteR)`;
these orderings must not be interchanged when forming the half tensors.

`LANC_GS_PREDICTION=T` supplies this candidate to ARPACK for iDMRG
with `LANC_NEIGEN=1`, in OBC or PBC. It does not require `WF_PREDICT_DIAG=T`.
The first two steps retain the usual solver initialization. See
`test/results/wf_obc/README.md` for the BHZ experiment: excellent overlap
and fewer products for one eigenstate, but no improvement for two.
The OBC experiment predates the streamed implementation. The projected
contraction is now shared by both boundary conditions, with link-aware index
maps and streamed construction of the trial vector.

## PBC single-step filter probe

`WF_PREDICT_DIAG=T` also tests three first-order filtered candidates in PBC,
`v_tau = normalize(v - tau*(Hv-<H>v))`, for tau=0.02,0.05,0.1. This uses
the one H*v already computed for the trial-state residual; no filtered
candidate is passed to ARPACK. The small BHZ PBC test in
`test/results/pbc_prediction_probe/` has GS weights approximately 0.026
before filtering and 0.074 at tau=0.1 in both n and p steps. This is too
weak to serve as the production PBC predictor.

## Correlated PBC prediction

The previous `n`/`p` link determines how `F` and `G` are sliced to recover
their site indices. The current link determines where each predicted
amplitude belongs in `sb_states`. These links alternate in PBC and must
remain distinct. `Cprev^+` uses the same SVD cutoff as OBC. The initial
step without both isometries keeps ARPACK's usual start vector.

With the BHZ driver, `LDMRG=10` (promoted to 11 for PBC), `MDMRG=20`,
`LANC_NEIGEN=1`, four MPI ranks, and `WF_PREDICT_DIAG=T`, the eight seeded
diagonalizations use 30-35 H*v products each; the corresponding unseeded
steps use 50-55. Energies per site agree within 3.2e-14. The first two
steps use the normal initialization. The diagnostic adds one H*v per step,
so its wall times are not the cost of a production run. The candidate's GS
weight is about 0.97-0.975 in these seeded steps. See
`test/results/pbc_correlated_probe/` for the reproducible inputs and logs.

At `MDMRG=80` the full run uses 370 H*v products without prediction and 290
with it. The worst energy/site difference is 1.8e-14. The projected-half
construction takes about 0.02 CPU seconds per step in this case.

The new memory use is `O(dim(SB)+d*M^2)` for the saved state and projected
half tensors, plus the existing superblock/eigensolver storage. The pair
matrix is only `M_L x M_R`. This has been tested through `M=80`; the speed
and peak memory at `M=500` remain unmeasured.
