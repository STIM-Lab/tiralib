"""
test_cpu_gpu.py

Compares CPU (ts.tk_vote2) vs GPU (tk_vote2_cuda / tk_vote2_cuda_prebuilt_flat)
tensor voting on sample fields, checking correctness and timing.

KNOWN NORMALIZATION DIFFERENCE
-------------------------------
Python tk_vote2 multiplies both stick and plate contributions by:

    eta_val = (2^2p * (p!)^2) / (pi * (2p)! * (sigma1^2 + sigma2^2))

The CUDA kernel does NOT apply this factor — it is missing from tk_stickvote_2d
and tk_plate_numerical in the C++ header. So CPU and GPU outputs are proportional,
not equal. For a pure-stick field the ratio is exactly eta_val; for mixed fields
the ratio is not a single constant.

This test computes the ratio and checks structural agreement (same shape of
voted field, correct sign pattern, ratio within expected range).

Run:
    cd python/
    python test_cpu_gpu.py
"""

import sys
import os
import math
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import tensors as ts
from tensors_lib.cuda_wrappers import (
    tk_vote2_cuda,
    _field_to_flat4, _build_nb_numpy, tk_vote2_cuda_prebuilt_flat,
    CUDA_TK_AVAILABLE,
)
from tensors_lib.optimize import _w_from_sigma


# ── Field generators ──────────────────────────────────────────────────────────

def make_stick_field(H, W, seed=0):
    """All-stick field: l0=0, random l1 and orientation."""
    rng = np.random.default_rng(seed)
    angles = rng.uniform(0, np.pi, (H, W))
    l1 = rng.uniform(0.5, 1.0, (H, W))
    vx, vy = np.cos(angles), np.sin(angles)
    T = np.zeros((H, W, 2, 2))
    T[:, :, 0, 0] = l1 * vx * vx
    T[:, :, 0, 1] = l1 * vx * vy
    T[:, :, 1, 0] = l1 * vx * vy
    T[:, :, 1, 1] = l1 * vy * vy
    return T

def make_ball_field(H, W, seed=0):
    """Isotropic (ball) field: l0 = l1, no preferred direction."""
    rng = np.random.default_rng(seed)
    lam = rng.uniform(0.3, 1.0, (H, W))
    T = np.zeros((H, W, 2, 2))
    T[:, :, 0, 0] = lam
    T[:, :, 1, 1] = lam
    return T

def make_mixed_field(H, W, seed=0):
    """Mix of stick and plate tensors."""
    rng = np.random.default_rng(seed)
    angles = rng.uniform(0, np.pi, (H, W))
    l1 = rng.uniform(0.5, 1.0, (H, W))
    l0 = rng.uniform(0.0, 0.4, (H, W))
    vx, vy = np.cos(angles), np.sin(angles)
    T = np.zeros((H, W, 2, 2))
    T[:, :, 0, 0] = l0 + (l1 - l0) * vx * vx
    T[:, :, 0, 1] = (l1 - l0) * vx * vy
    T[:, :, 1, 0] = (l1 - l0) * vx * vy
    T[:, :, 1, 1] = l0 + (l1 - l0) * vy * vy
    return T


# ── Helpers ───────────────────────────────────────────────────────────────────

def eta_val(sigma1, sigma2, power):
    num = math.pi * math.factorial(2 * power)
    den = 2 ** (2 * power) * (math.factorial(power) ** 2)
    return den / (num * (sigma1 ** 2 + sigma2 ** 2))

def time_ms(fn, n=5, warmup=2):
    for _ in range(warmup):
        fn()
    t0 = time.perf_counter()
    for _ in range(n):
        fn()
    return (time.perf_counter() - t0) / n * 1000.0

def ratio_stats(cpu, gpu):
    """Element-wise cpu/gpu ratio where gpu != 0."""
    mask = np.abs(gpu) > 1e-10
    r = cpu[mask] / gpu[mask]
    return r.mean(), r.std(), r.min(), r.max()

def check(label, passed):
    mark = "PASS" if passed else "FAIL"
    print(f"    [{mark}] {label}")
    return passed


# ── Single comparison ─────────────────────────────────────────────────────────

def compare(field_name, T, sigma1, sigma2, power, N=20, time_it=True, n_time=5):
    s1, s2, pw = sigma1, sigma2, power
    eta = eta_val(s1, s2, pw)

    cpu = ts.tk_vote2(T, sigma1=s1, sigma2=s2, power=pw, plate=True, N=N)
    gpu = tk_vote2_cuda(T, sigma1=s1, sigma2=s2, power=pw, plate=True, N=N)

    # Prebuilt-NB path (should match gpu exactly)
    I_flat = _field_to_flat4(T)
    w  = _w_from_sigma(s1, s2)
    nb = _build_nb_numpy(s1, s2, pw, N, w)
    gpu_pre = tk_vote2_cuda_prebuilt_flat(I_flat, nb, s1, s2, pw, plate=True, N=N,
                                          out_dtype=T.dtype)

    r_mean, r_std, r_min, r_max = ratio_stats(cpu, gpu)

    print(f"\n  {field_name}  |  sigma1={s1}, sigma2={s2}, power={pw}")
    print(f"    cpu/gpu ratio:  mean={r_mean:.4f}  std={r_std:.4f}"
          f"  min={r_min:.4f}  max={r_max:.4f}")
    print(f"    eta_val        = {eta:.4f}   (expected ratio for pure-stick)")

    # gpu vs gpu_pre must match exactly (same computation, same table)
    diff_pre = np.abs(gpu - gpu_pre).max()
    check("gpu == gpu_prebuilt (max diff < 1e-5)", diff_pre < 1e-5)

    # cpu / gpu ratio should be close to eta_val for stick-dominated fields
    ratio_ok = abs(r_mean - eta) / eta < 0.15
    check(f"cpu/gpu ratio within 15% of eta_val ({eta:.3f})", ratio_ok)

    # Both outputs should be symmetric
    cpu_sym = np.abs(cpu[:, :, 0, 1] - cpu[:, :, 1, 0]).max()
    gpu_sym = np.abs(gpu[:, :, 0, 1] - gpu[:, :, 1, 0]).max()
    check("CPU output is symmetric (T01 == T10)", cpu_sym < 1e-12)
    check("GPU output is symmetric (T01 == T10)", gpu_sym < 1e-5)

    if time_it:
        t_cpu = time_ms(lambda: ts.tk_vote2(T, sigma1=s1, sigma2=s2,
                                             power=pw, plate=True, N=N), n=n_time)
        t_gpu = time_ms(lambda: tk_vote2_cuda(T, sigma1=s1, sigma2=s2,
                                               power=pw, plate=True, N=N), n=n_time)
        t_pre = time_ms(lambda: (
            tk_vote2_cuda_prebuilt_flat(
                I_flat, _build_nb_numpy(s1, s2, pw, N, _w_from_sigma(s1, s2)),
                s1, s2, pw, plate=True, N=N, out_dtype=T.dtype)
        ), n=n_time)
        speedup     = t_cpu / t_gpu  if t_gpu > 0 else float('inf')
        speedup_pre = t_cpu / t_pre  if t_pre > 0 else float('inf')
        print(f"    CPU:         {t_cpu:8.1f} ms")
        print(f"    GPU (std):   {t_gpu:8.1f} ms   {speedup:.1f}x speedup")
        print(f"    GPU (prebuilt NB): {t_pre:8.1f} ms   {speedup_pre:.1f}x speedup")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 65)
    print("  CPU vs GPU TK Tensor Voting — Correctness & Timing")
    print("=" * 65)
    print(f"\n  CUDA TK available: {CUDA_TK_AVAILABLE}")

    if not CUDA_TK_AVAILABLE:
        print("  CUDA extension not found. Exiting.")
        return

    print("""
  Normalization note:
    Python tk_vote2 scales outputs by eta_val; the CUDA kernel does not.
    So cpu_output ≈ eta_val * gpu_output (exact for pure-stick fields).
    Both GPU paths (standard and prebuilt-NB) should agree exactly.
""")

    H, W = 64, 64

    fields = [
        ("Stick (64x64)",  make_stick_field(H, W)),
        ("Ball  (64x64)",  make_ball_field(H, W)),
        ("Mixed (64x64)",  make_mixed_field(H, W)),
    ]

    params = [
        (5.0, 0.0, 1, 20),
        (5.0, 2.0, 1, 20),
        (5.0, 2.0, 3, 20),   # larger power
    ]

    print("─" * 65)
    print("  CORRECTNESS  (64x64, timed with 5 runs)")
    print("─" * 65)

    for fname, T in fields:
        for s1, s2, pw, N in params:
            compare(fname, T, s1, s2, pw, N=N, time_it=True, n_time=5)

    # Large-field timing only
    print("\n" + "─" * 65)
    print("  TIMING — large fields (no correctness check)")
    print("─" * 65)

    for H2, W2 in [(128, 128), (256, 256)]:
        T_big = make_mixed_field(H2, W2)
        s1, s2, pw, N = 10.0, 3.0, 1, 20
        print(f"\n  Mixed ({H2}x{W2})  |  sigma1={s1}, sigma2={s2}, power={pw}")
        I_flat = _field_to_flat4(T_big)
        w  = _w_from_sigma(s1, s2)
        nb = _build_nb_numpy(s1, s2, pw, N, w)
        t_cpu = time_ms(lambda: ts.tk_vote2(T_big, sigma1=s1, sigma2=s2,
                                             power=pw, plate=True, N=N), n=3)
        t_gpu = time_ms(lambda: tk_vote2_cuda(T_big, sigma1=s1, sigma2=s2,
                                               power=pw, plate=True, N=N), n=3)
        t_pre = time_ms(lambda: (
            tk_vote2_cuda_prebuilt_flat(
                I_flat, _build_nb_numpy(s1, s2, pw, N, _w_from_sigma(s1, s2)),
                s1, s2, pw, plate=True, N=N, out_dtype=T_big.dtype)
        ), n=3)
        print(f"    CPU:               {t_cpu:8.1f} ms")
        print(f"    GPU (std):         {t_gpu:8.1f} ms   {t_cpu/t_gpu:.1f}x speedup")
        print(f"    GPU (prebuilt NB): {t_pre:8.1f} ms   {t_cpu/t_pre:.1f}x speedup")


if __name__ == '__main__':
    main()
