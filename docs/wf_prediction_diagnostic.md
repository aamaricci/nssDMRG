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

La contrazione diretta del centro costa `O(dim(SB)*M_L*M_R)` e viene
saltata quando la stima supera 50 milioni di termini. Prima di usarla
su taglie grandi occorre una contrazione per settori. Il parametro
`WF_PREDICT_DIAG` e' falso per default.

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

## Correlated OBC prototype

For OBC growth the retained half states are contracted as
`F=UL^dagger Psi`, `G=Psi UR*`, and
`Theta(a,s,t,b)=sum F(a,s,beta) Cprev^+(beta,alpha) G(alpha,t,b)`.
`Cprev` refers to the actual old block bases, with no eigenvector matching.
The pseudoinverse uses SVD and a relative cutoff 1e-10. No explicit extra
fermionic permutation is applied to this same-edge contraction; the
ordered tensor coefficients are carried by the half states.

`LANC_GS_PREDICTION=T` supplies this candidate to ARPACK only for OBC iDMRG
with `LANC_NEIGEN=1`. It does not require `WF_PREDICT_DIAG=T`.
The first two steps retain the usual solver initialization. See
`test/results/wf_obc/README.md` for the BHZ experiment: excellent overlap
and fewer products for one eigenstate, but no improvement for two.
The existing work cap and dense temporary arrays still limit this to
small experiments. Multi-rank execution has not been validated.
