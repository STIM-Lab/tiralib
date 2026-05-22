"""
Stage 3 benchmark: tk_vote2_fast vs ts.tk_vote2 (N=0).

1. Correctness check on a small field.
2. Timing comparison on a 64×64 field.
"""

import sys, os, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import tensors as ts
from tensors_lib.optimize import tk_vote2_fast

np.random.seed(42)

def make_field(H, W, theta=np.pi / 4):
    """Uniform stick field pointing at angle theta."""
    T = np.zeros((H, W, 2, 2))
    qx, qy = np.cos(theta), np.sin(theta)
    T[:, :, 0, 0] = qx * qx
    T[:, :, 0, 1] = qx * qy
    T[:, :, 1, 0] = qy * qx
    T[:, :, 1, 1] = qy * qy
    return T


def make_random_field(H, W):
    """Random PSD tensor field."""
    T = np.zeros((H, W, 2, 2))
    for i in range(H):
        for j in range(W):
            v = np.random.randn(2)
            T[i, j] = np.outer(v, v)
    return T


# ─────────────────────────────────────────────────────
# CORRECTNESS — small field, deterministic
# ─────────────────────────────────────────────────────
print("=" * 60)
print("STAGE 3: Correctness check (10×10 uniform field)")
print("=" * 60)

T_small = make_field(10, 10, theta=np.pi / 3)
sigma1, sigma2, power = 3.0, 1.0, 2

ref  = ts.tk_vote2(T_small, sigma1=sigma1, sigma2=sigma2, power=power, plate=True, N=0)
fast = tk_vote2_fast(T_small, sigma1=sigma1, sigma2=sigma2, power=power)

diff     = np.abs(ref - fast)
denom    = np.maximum(np.abs(ref), 1e-12)
reldiff  = diff / denom

print(f"sigma1={sigma1}, sigma2={sigma2}, power={power}")
print(f"Max absolute error : {diff.max():.4e}")
print(f"Max relative error : {reldiff.max():.4e}")
print(f"Mean absolute error: {diff.mean():.4e}")
cx, cy = 5, 5
print(f"\nCentral pixel ref:\n{np.round(ref[cx,cy], 6)}")
print(f"Central pixel fast:\n{np.round(fast[cx,cy], 6)}")

tol_abs = 1e-8
tol_rel = 1e-6
if diff.max() < tol_abs and reldiff.max() < tol_rel:
    print("\nCorrectness: PASS")
else:
    print(f"\nCorrectness: FAIL (abs {diff.max():.2e}, rel {reldiff.max():.2e})")

# ─────────────────────────────────────────────────────
# CORRECTNESS — mixed random field
# ─────────────────────────────────────────────────────
print()
print("=" * 60)
print("STAGE 3: Correctness check (15×15 random PSD field)")
print("=" * 60)

T_rand = make_random_field(15, 15)
ref_r   = ts.tk_vote2(T_rand, sigma1=2.0, sigma2=0.5, power=1, plate=True, N=0)
fast_r  = tk_vote2_fast(T_rand, sigma1=2.0, sigma2=0.5, power=1)

diff_r   = np.abs(ref_r - fast_r)
denom_r  = np.maximum(np.abs(ref_r), 1e-12)
reldiff_r = diff_r / denom_r

print(f"Max absolute error : {diff_r.max():.4e}")
print(f"Max relative error : {reldiff_r.max():.4e}")
if diff_r.max() < tol_abs and reldiff_r.max() < tol_rel:
    print("Correctness: PASS")
else:
    print(f"Correctness: FAIL (abs {diff_r.max():.2e}, rel {reldiff_r.max():.2e})")

# ─────────────────────────────────────────────────────
# BENCHMARK — 64×64 field
# ─────────────────────────────────────────────────────
print()
print("=" * 60)
print("STAGE 3: Benchmark (64×64 field)")
print("=" * 60)

T64 = make_field(64, 64)
sigma1_b, sigma2_b, power_b = 5.0, 1.0, 2

# Warm up
_ = ts.tk_vote2(T64, sigma1=sigma1_b, sigma2=sigma2_b, power=power_b, plate=True, N=0)
_ = tk_vote2_fast(T64, sigma1=sigma1_b, sigma2=sigma2_b, power=power_b)

# Time reference (N=0 original)
N_REPS = 3
t_ref_total = 0.0
for _ in range(N_REPS):
    t0 = time.perf_counter()
    out_ref = ts.tk_vote2(T64, sigma1=sigma1_b, sigma2=sigma2_b, power=power_b, plate=True, N=0)
    t_ref_total += time.perf_counter() - t0
t_ref = t_ref_total / N_REPS

# Time fast version
t_fast_total = 0.0
for _ in range(N_REPS):
    t0 = time.perf_counter()
    out_fast = tk_vote2_fast(T64, sigma1=sigma1_b, sigma2=sigma2_b, power=power_b)
    t_fast_total += time.perf_counter() - t0
t_fast = t_fast_total / N_REPS

speedup = t_ref / t_fast if t_fast > 0 else float('inf')
diff64  = np.abs(out_ref - out_fast)

print(f"sigma1={sigma1_b}, sigma2={sigma2_b}, power={power_b}")
print(f"tk_vote2     (original, N=0): {t_ref*1000:.1f} ms  (mean over {N_REPS} reps)")
print(f"tk_vote2_fast (optimized)   : {t_fast*1000:.1f} ms  (mean over {N_REPS} reps)")
print(f"Speedup: {speedup:.2f}×")
print(f"Max abs error (64×64 check): {diff64.max():.4e}")

# ─────────────────────────────────────────────────────
# BENCHMARK — 32×32 (typical optimizer call size)
# ─────────────────────────────────────────────────────
print()
print("=" * 60)
print("STAGE 3: Benchmark (32×32 field — optimizer typical)")
print("=" * 60)

T32 = make_field(32, 32)
t_ref_total = 0.0
for _ in range(N_REPS):
    t0 = time.perf_counter()
    _ = ts.tk_vote2(T32, sigma1=sigma1_b, sigma2=sigma2_b, power=power_b, plate=True, N=0)
    t_ref_total += time.perf_counter() - t0

t_fast_total = 0.0
for _ in range(N_REPS):
    t0 = time.perf_counter()
    _ = tk_vote2_fast(T32, sigma1=sigma1_b, sigma2=sigma2_b, power=power_b)
    t_fast_total += time.perf_counter() - t0

t_ref32  = t_ref_total / N_REPS
t_fast32 = t_fast_total / N_REPS
speedup32 = t_ref32 / t_fast32 if t_fast32 > 0 else float('inf')
print(f"tk_vote2     (original, N=0): {t_ref32*1000:.1f} ms")
print(f"tk_vote2_fast (optimized)   : {t_fast32*1000:.1f} ms")
print(f"Speedup: {speedup32:.2f}×")
