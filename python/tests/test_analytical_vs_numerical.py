"""
Stage 1: Sanity check — analytical (N=0) vs numerical (N>0) plate vote in tk_vote2.

Tests whether the closed-form plate field (N=0) agrees with the numerical
integration over N stick votes (N>0) for identical sigma1, sigma2, power.
"""

import sys
import os
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import tensors as ts

# ──────────────────────────────────────────────
# 1. Build a small, controlled synthetic tensor field
# ──────────────────────────────────────────────
H, W = 10, 10
T = np.zeros((H, W, 2, 2))

# All tensors pointing in the same direction: theta = pi/4
theta = np.pi / 4
qx, qy = np.cos(theta), np.sin(theta)
# Stick tensor: outer product of q with itself (eigenvalue 1 in direction q)
T[:, :, 0, 0] = qx * qx
T[:, :, 0, 1] = qx * qy
T[:, :, 1, 0] = qy * qx
T[:, :, 1, 1] = qy * qy

# ──────────────────────────────────────────────
# 2. Parameters
# ──────────────────────────────────────────────
sigma1 = 3.0
sigma2 = 1.0
power  = 2

# ──────────────────────────────────────────────
# 3. Run both modes
# ──────────────────────────────────────────────
t0 = time.perf_counter()
out_analytical = ts.tk_vote2(T, sigma1=sigma1, sigma2=sigma2, power=power, plate=True, N=0)
t_analytical = time.perf_counter() - t0

t0 = time.perf_counter()
out_numerical_20 = ts.tk_vote2(T, sigma1=sigma1, sigma2=sigma2, power=power, plate=True, N=20)
t_numerical_20 = time.perf_counter() - t0

t0 = time.perf_counter()
out_numerical_50 = ts.tk_vote2(T, sigma1=sigma1, sigma2=sigma2, power=power, plate=True, N=50)
t_numerical_50 = time.perf_counter() - t0

# ──────────────────────────────────────────────
# 4. Compare
# ──────────────────────────────────────────────
diff_20  = np.abs(out_analytical - out_numerical_20)
diff_50  = np.abs(out_analytical - out_numerical_50)

# Relative error: |A - B| / max(|A|, eps)
denom = np.maximum(np.abs(out_analytical), 1e-12)
reldiff_20 = diff_20 / denom
reldiff_50 = diff_50 / denom

print("=" * 60)
print("STAGE 1: Analytical (N=0) vs Numerical (N>0) Plate Vote")
print("=" * 60)
print(f"Field size:  {H}x{W}")
print(f"sigma1={sigma1}, sigma2={sigma2}, power={power}")
print()
print(f"N=0  (analytical): {t_analytical*1000:.1f} ms")
print(f"N=20 (numerical) : {t_numerical_20*1000:.1f} ms")
print(f"N=50 (numerical) : {t_numerical_50*1000:.1f} ms")
print()
print(f"Max absolute error (N=20): {diff_20.max():.6e}")
print(f"Max relative error (N=20): {reldiff_20.max():.6e}")
print(f"Mean absolute error(N=20): {diff_20.mean():.6e}")
print()
print(f"Max absolute error (N=50): {diff_50.max():.6e}")
print(f"Max relative error (N=50): {reldiff_50.max():.6e}")
print(f"Mean absolute error(N=50): {diff_50.mean():.6e}")

# ──────────────────────────────────────────────
# 5. Spot-check: print the central pixel tensor for both
# ──────────────────────────────────────────────
cx, cy = H // 2, W // 2
print()
print(f"Central pixel [{cx},{cy}] — out_analytical:")
print(np.round(out_analytical[cx, cy], 6))
print(f"Central pixel [{cx},{cy}] — out_numerical (N=20):")
print(np.round(out_numerical_20[cx, cy], 6))
print(f"Central pixel [{cx},{cy}] — diff:")
print(np.round(diff_20[cx, cy], 6))

# ──────────────────────────────────────────────
# 6. Verdict
# ──────────────────────────────────────────────
print()
threshold_abs = 1e-3
threshold_rel = 1e-2
ok = diff_20.max() < threshold_abs and reldiff_20.max() < threshold_rel

if ok:
    print(f"VERDICT: PASS — N=0 and N=20 agree to max abs err {diff_20.max():.2e} "
          f"and max rel err {reldiff_20.max():.2e}.")
    print("  => Safe to use N=0 (analytical) in the optimizer.")
else:
    print(f"VERDICT: FAIL — discrepancy exceeds threshold.")
    print(f"  max abs err = {diff_20.max():.2e} (threshold {threshold_abs:.0e})")
    print(f"  max rel err = {reldiff_20.max():.2e} (threshold {threshold_rel:.0e})")
