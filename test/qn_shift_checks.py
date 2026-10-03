"""Small-chain exact spectra, charge measurements, and restart compatibility.

Optional third argument: driver linked to mainSep25, to verify old checkpoints.
All files are confined to the requested scratch directory.
"""
import math
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile

exe = str(Path(sys.argv[1]).resolve())
root = Path(sys.argv[2]).resolve()
old_exe = str(Path(sys.argv[3]).resolve()) if len(sys.argv) > 3 else None
root.mkdir(parents=True, exist_ok=True)


def config(model, density, offset, length=2, mode="normal", legacy=False, run=True):
    qdim = len(density)
    records = (f"QNTYPE={'global' if model == 0 else 'local'}\n"
               f"DMRG_QN={','.join(map(str, offset if model == 0 else density))}\n"
               if legacy else
               f"QN_DENSITY={','.join(map(str, density))}\n"
               f"QN_OFFSET={','.join(map(str, offset))}\n")
    return f"""TEST_MODEL={model}
DMRGTYPE=i
LDMRG={length}
MDMRG=256
EDMRG=0
MSWEEP=256
ESWEEP=0
QNTRUNCATION_ERROR=0
QNTRUNCATION_DIM=0
QNDIM={qdim}
DMRG_MODE={mode}
NORB=1
ULOC=0
HFMODE=F
JX=1
JP=1
LANC_NEIGEN=1
LANC_DIM_THRESHOLD=4096
SPARSE_H=F
SAVE_BLOCK=T
SAVE_ALL_BLOCKS=T
SAVE_UMAT=T
BLOCK_UMAT_CACHE=T
SAVE_MEASURE_STATE=T
IRUN={'T' if run else 'F'}
IMEASURE=T
""" + records


def run_case(name, text, executable=exe, restart=None, fail=False):
    work = root / name
    if work.exists():
        shutil.rmtree(work)
    work.mkdir()
    if restart:
        with tarfile.open(restart) as archive:
            archive.extractall(work, filter="data")
        dirs = list(work.glob("restart_*"))
        assert len(dirs) == 1, dirs
        dirs[0].rename(work / "restart")
    (work / "DMRG.conf").write_text(text)
    result = subprocess.run([executable, "FINPUT=DMRG.conf"], cwd=work,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=90)
    (work / "run.log").write_text(result.stdout)
    if fail:
        assert result.returncode != 0, name
        return result.stdout
    assert result.returncode == 0, f"{name}: {result.stdout[-4000:]}"
    values = [list(map(float, line.split()))
              for line in (work / "qn_result.out").read_text().splitlines()]
    return values, next(work.glob("restart_*.tgz"), None)


def check(values, energy, charges):
    assert abs(values[0][0] - energy) < 1e-8, (values, energy)
    assert len(values[1]) == len(charges)
    assert all(abs(a-b) < 1e-8 for a, b in zip(values[1], charges)), (values, charges)


def free_energy(n, up, down):
    levels = [-2*math.cos(k*math.pi/(n+1)) for k in range(1, n+1)]
    return sum(levels[:up]) + sum(levels[:down])


def spin_energy(n, magnetization):
    # Independent ED of H=sum_i S_i.S_{i+1}, in a bit basis.
    states = [bits for bits in range(1 << n)
              if bits.bit_count() == n//2 + magnetization]
    index = {bits: i for i, bits in enumerate(states)}
    dim = len(states)
    h = [[0.0]*dim for _ in states]
    for bits, i in index.items():
        for site in range(n-1):
            opposite = ((bits >> site) ^ (bits >> (site+1))) & 1
            h[i][i] += -.25 if opposite else .25
            if opposite:
                h[i][index[bits ^ (3 << site)]] += .5
    # Jacobi diagonalization keeps the test independent of external Python packages.
    for _ in range(100*dim*dim):
        if dim == 1:
            break
        p, q = max(((i, j) for i in range(dim) for j in range(i+1, dim)),
                   key=lambda pair: abs(h[pair[0]][pair[1]]))
        if abs(h[p][q]) < 1e-13:
            break
        angle = .5*math.atan2(2*h[p][q], h[q][q]-h[p][p])
        c, t = math.cos(angle), math.sin(angle)
        app, aqq, apq = h[p][p], h[q][q], h[p][q]
        for k in range(dim):
            if k in (p, q):
                continue
            hp, hq = h[k][p], h[k][q]
            h[k][p] = h[p][k] = c*hp-t*hq
            h[k][q] = h[q][k] = t*hp+c*hq
        h[p][p] = c*c*app+t*t*aqq-2*c*t*apq
        h[q][q] = t*t*app+c*c*aqq+2*c*t*apq
        h[p][q] = h[q][p] = 0
    else:
        raise AssertionError("ED did not converge")
    return min(h[i][i] for i in range(dim))


for shift in [0, 1, -1, 2]:
    values, _ = run_case(f"spin_{shift}", config(0, [0], [shift]))
    check(values, spin_energy(4, shift), [shift])

values, _ = run_case("spin_density", config(0, [-.25], [0]))
check(values, spin_energy(4, -1), [-1])
values, _ = run_case("spin_one", config(2, [0], [4]))
check(values, 3, [4])

for name, rho, offset, charges in [
        ("added_up", [.5, .5], [1, 0], [3, 2]),
        ("removed_down", [.5, .5], [0, -1], [2, 1]),
        ("spin_transfer", [.5, .5], [1, -1], [3, 1]),
        ("quarter", [.25, .25], [0, 0], [1, 1]),
        ("third_truncated", [.33333333, .33333333], [0, 0], [1, 1]),
        ("absolute", [0, 0], [2, 1], [2, 1])]:
    values, _ = run_case(name, config(1, rho, offset))
    check(values, free_energy(4, *charges), charges)

values, _ = run_case("third_float", config(1, [.33333333]*2, [0, 0], length=3))
check(values, free_energy(6, 1, 1), [1, 1])

for mode, rho, shift in [("superc", 0, 2), ("nonsu2", 1, 1)]:
    values, _ = run_case(mode, config(1, [rho], [shift], mode=mode))
    target = int(4*rho)+shift
    sectors = [(free_energy(4, u, d), u, d) for u in range(5) for d in range(5)
               if (u-d if mode == "superc" else u+d) == target]
    expected = min(sectors)[0]
    assert abs(values[0][0]-expected) < 1e-8, values
    actual = values[1][0]
    assert abs(actual-target) < 1e-8, values

# Check the target against the actual final superblock length in a finite sweep.
finite, _ = run_case("finite", config(0, [.25], [0], length=3).replace("DMRGTYPE=i", "DMRGTYPE=f"))
actual_length = int(finite[2][0])
assert abs(finite[1][0]-int(actual_length*.25)) < 1e-8, finite

values, _ = run_case("warmup_spin", config(0, [0], [3], length=3))
check(values, 1.25, [3])
assert "Warm-up QN clipping at length:" in (root / "warmup_spin" / "run.log").read_text()
values, _ = run_case("warmup_fermion", config(1, [.5, .5], [3, 0], length=3))
check(values, free_energy(6, 6, 3), [6, 3])
assert "Warm-up QN clipping at length:" in (root / "warmup_fermion" / "run.log").read_text()
message = run_case("impossible_final", config(0, [0], [4], length=3), fail=True)
assert "requested QN sector is empty" in message

message = run_case("empty", config(1, [.5, .5], [3, 0]), fail=True)
assert "requested QN sector is empty" in message
message = run_case("fractional", config(1, [.5, .5], [.5, 0]), fail=True)
assert "requested QN sector is empty" in message
message = run_case("obsolete_input", config(0, [0], [0], legacy=True), fail=True)
assert "migrate the old quantum-number input" in message

for model, rho, offset in [(0, [0], [1]), (1, [.5, .5], [0, 0])]:
    seed_exe = old_exe or exe
    def restart_config(**kwargs):
        text = config(model, rho, offset, **kwargs)
        return text.replace("ULOC=0", "ULOC=2") if model == 1 else text
    seed, archive = run_case(f"seed_{model}",
        restart_config(length=3, legacy=bool(old_exe)), seed_exe)
    post, _ = run_case(f"post_{model}",
        restart_config(length=3, run=False), restart=archive)
    check(post, seed[0][0], seed[1])
    resumed, _ = run_case(f"resume_{model}",
        restart_config(length=4), restart=archive)
    fresh, _ = run_case(f"fresh_{model}", restart_config(length=4))
    check(resumed, fresh[0][0], fresh[1])

print("QN spectra, charges, rejected inputs, post-processing and restart checks passed"
      + (" (checkpoints generated by mainSep25)" if old_exe else ""))
