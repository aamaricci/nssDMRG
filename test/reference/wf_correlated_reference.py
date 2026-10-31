"""Exact OBC reference for the half-wavefunction prediction contraction.

This is a reference experiment, not the alternating-edge PBC implementation.
Uses occupation/spin product order with the rightmost site fastest.
"""
import numpy as np


def ground_state(n):
    # XX chain, fixed half filling. Nearest-neighbor fermion hopping has
    # this same matrix in the ordered occupation basis for open boundaries.
    states = [x for x in range(1 << n) if x.bit_count() == n // 2]
    index = {x: k for k, x in enumerate(states)}
    h = np.zeros((len(states), len(states)))
    for k, x in enumerate(states):
        for j in range(n - 1):
            if ((x >> j) & 1) != ((x >> (j + 1)) & 1):
                y = x ^ (3 << j)
                h[index[y], k] = -1.0
    energy, vectors = np.linalg.eigh(h)
    psi = np.zeros(1 << n)
    psi[states] = vectors[:, 0]
    return psi, energy[0], h, states


def schmidt(psi):
    u, s, vh = np.linalg.svd(psi, full_matrices=False)
    return u, np.diag(s), vh.T.conj()


def predict(psi, ul, ur, previous_center, cutoff=1e-12):
    """psi has axes (oldL, siteL, siteR, oldR).

    F = UL^dagger Psi has axes (newL, siteR, oldR).
    G = Psi UR* has axes (oldL, siteL, newR).
    The inverse center joins oldR to oldL; neither basis is matched
    by sorting eigenvalues or by fixing eigenvector signs.
    """
    ml, d, _, mr = psi.shape
    matrix = psi.reshape(ml*d, d*mr)
    f = (ul.conj().T @ matrix).reshape(ul.shape[1], d, mr)
    g = (matrix @ ur.conj()).reshape(ml, d, ur.shape[1])
    inverse = np.linalg.pinv(previous_center, rcond=cutoff)
    return np.einsum('asb,bc,ctd->astd', f, inverse, g)


def main():
    p4, _, _, _ = ground_state(4)
    l4, c4, r4 = schmidt(p4.reshape(4, 4))
    p6, _, _, _ = ground_state(6)
    # Transform the exact six-site state into the retained four-site bases.
    left_embedding = np.kron(l4, np.eye(2))
    right_embedding = np.kron(np.eye(2), r4)
    psi6 = left_embedding.conj().T @ p6.reshape(8, 8) @ right_embedding.conj()
    l6, c6, r6 = schmidt(psi6)
    trial = predict(psi6.reshape(4, 2, 2, 4), l6, r6, c4)
    l6physical = left_embedding @ l6
    r6physical = right_embedding @ r6
    p8, e8, h8, states8 = ground_state(8)
    new_left = np.kron(l6physical, np.eye(2))
    new_right = np.kron(np.eye(2), r6physical)
    physical = (new_left @ trial.reshape(16, 16) @ new_right.T).reshape(-1)
    physical /= np.linalg.norm(physical)
    outside = physical.copy()
    outside[states8] = 0
    assert np.linalg.norm(outside) < 1e-10, 'prediction changed particle number'
    v = physical[states8]
    energy = np.vdot(v, h8 @ v).real
    residual = np.linalg.norm(h8 @ v - energy*v)
    overlap = abs(np.vdot(p8, physical))**2
    print(f'exact E={e8:.12g}; trial E={energy:.12g}')
    print(f'weight={overlap:.12g}; residual={residual:.12g}')
    exact_s = np.linalg.svd(p8.reshape(16,16), compute_uv=False)
    trial_s = np.linalg.svd(physical.reshape(16,16), compute_uv=False)
    print(f'Schmidt probability L1 error={np.sum(abs(exact_s**2-trial_s**2)):.12g}')
    assert overlap > 0.9, 'same-edge prediction failed on the exact reference'

    # Independent gauge test: rotate both OLD virtual bases by arbitrary
    # complex unitaries, consistently transforming psi and the old center.
    rng = np.random.default_rng(41)
    ql = np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))[0]
    qr = np.linalg.qr(rng.normal(size=(4,4))+1j*rng.normal(size=(4,4)))[0]
    tl = np.kron(ql, np.eye(2))
    tr = np.kron(np.eye(2), qr)
    rotated_psi = tl.conj().T @ psi6 @ tr.conj()
    rotated_center = ql.conj().T @ c4 @ qr.conj()
    rotated_trial = predict(rotated_psi.reshape(4,2,2,4),
                            tl.conj().T @ l6, tr.conj().T @ r6,
                            rotated_center)
    err = np.linalg.norm(trial-rotated_trial)/np.linalg.norm(trial)
    print(f'complex gauge covariance error={err:.3e}')
    assert err < 1e-10


if __name__ == '__main__':
    main()
