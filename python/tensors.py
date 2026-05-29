import tensors_lib.phantoms as phantoms
import tensors_lib.noise as noise
import tensors_lib.metrics as metrics
import tensors_lib.plots as plots
import numpy as np
import math

def eigmag(T):
    
    eigenValues, eigenVectors = np.linalg.eigh(T)
    magValues = np.abs(eigenValues)
    
    idx = np.argsort(magValues, -1)
    
    sortedValues = np.take_along_axis(eigenValues, idx, -1)
    sortedVectors = np.zeros_like(eigenVectors)
    sortedVectors[..., 0, :] = np.take_along_axis(eigenVectors[..., 0, :], idx, -1)
    sortedVectors[..., 1, :] = np.take_along_axis(eigenVectors[..., 1, :], idx, -1)
    
    return sortedValues, sortedVectors

# calculate the structure tensor for the input array, optionally adding noise to the tensor components
def structure(image, noise=0.0):
    # Calculate the image gradient
    dIdy = np.gradient(image, axis=0, edge_order=2)
    dIdx = np.gradient(image, axis=1, edge_order=2)
    
    # Create the structure tensor
    T = np.zeros((image.shape[0], image.shape[1], 2, 2))
    
    T[:, :, 0, 0] = dIdx * dIdx        
    T[:, :, 1, 1] = dIdy * dIdy
    T[:, :, 0, 1] = dIdx * dIdy        

    if noise > 0:
        T[:, :, 0, 0] = T[:, :, 0, 0] + np.abs(np.random.normal(0.0, noise, T[:, :, 0, 0].shape))
        T[:, :, 1, 1] = T[:, :, 1, 1] + np.abs(np.random.normal(0.0, noise, T[:, :, 1, 1].shape))
        T[:, :, 0, 1] = T[:, :, 0, 1] + np.random.normal(0.0, noise, T[:, :, 0, 1].shape)

    T[:, :, 1, 0] = T[:, :, 0, 1]

    return T

def _eccentricity(evals):
    """Ellipse eccentricity sqrt(1 - (b/a)^2) from eigh-sorted eigenvalues (ascending)."""
    a = np.abs(evals[..., 1])       # largest eigenvalue
    b = np.abs(evals[..., 0])       # smallest eigenvalue
    safe_a = np.where(a > 1e-12, a, 1.0)
    ratio  = np.clip(b / safe_a, 0.0, 1.0)
    return np.where(a > 1e-12, np.sqrt(1.0 - ratio ** 2), 0.0) # or one?


def blur(T, sigma):
    """Gaussian blur of tensor field T. sigma=0 returns T unchanged."""
    if sigma == 0:
        return T.copy()
    import scipy.ndimage
    return scipy.ndimage.gaussian_filter(T, sigma=(sigma, sigma, 0, 0), mode="constant")


def atv_vote2(T, sigma=3):
    """Executes ATV over an entire input tensor field T."""
    if sigma == 0:
        return T.copy()
    # Perform eigendecomposition of the field
    evals, evecs = np.linalg.eigh(T)
    
    # Calculate optimal window size
    pad = int(3 * sigma)
    w = 2 * pad + 1
    x = np.linspace(-pad, pad, w)
    X0, X1 = np.meshgrid(x, x)
    
    # Precompute distance and attenuation (c) for the grid
    dist_sq = X0**2 + X1**2
    dist = np.sqrt(dist_sq)
    
    c = np.zeros_like(dist)
    if sigma != 0:
        c = np.exp(-dist_sq / (sigma**2))
    else:
        c[dist < 1e-12] = 1.0

    # Precompute normalized direction vectors (dx, dy)
    dx = np.divide(X0, dist, out=np.zeros_like(X0), where=dist!=0)
    dy = np.divide(X1, dist, out=np.zeros_like(X1), where=dist!=0)

    # Precompute Rotation Matrix R components and their combinations
    R_a = 1.0 - 2.0 * dx**2
    R_b = -2.0 * dx * dy
    R_c = 1.0 - 2.0 * dy**2

    Ra2, Rb2, Rc2 = R_a**2, R_b**2, R_c**2
    Rab, Rbc, Rac = R_a * R_b, R_b * R_c, R_a * R_c

    # Create padded vote field
    VF = np.pad(np.zeros(T.shape), ((pad, pad), (pad, pad), (0, 0), (0, 0)))
    
    inv3 = 1.0 / 3.0
    inv4 = 1.0 / 4.0

    for x0 in range(T.shape[0]):
        for x1 in range(T.shape[1]):
            lambdas = evals[x0, x1]
            # Optimization: Skip completely empty tensors
            if np.sum(np.abs(lambdas)) > 0:
                evec = evecs[x0, x1]
                
                # Sort eigenvalues to ensure l0 <= l1
                idx = np.argsort(lambdas)
                l0, l1 = lambdas[idx[0]], lambdas[idx[1]]
                e0, e1 = evec[:, idx[0]], evec[:, idx[1]]

                # Base normalized elementary tensors
                T1_a, T1_b, T1_c = e1[0]**2, e1[0]*e1[1], e1[1]**2
                T0_a, T0_b, T0_c = e0[0]**2, e0[0]*e0[1], e0[1]**2
                T2_a, T2_b, T2_c = T0_a + T1_a, T0_b + T1_b, T0_c + T1_c

                # Internal terms: v^T * T * v
                vT1v = dx * (dx * T1_a + dy * T1_b) + dy * (dx * T1_b + dy * T1_c)
                vT1vT1_a, vT1vT1_b, vT1vT1_c = T1_a * vT1v, T1_b * vT1v, T1_c * vT1v
                
                vT2v = dx * (dx * T2_a + dy * T2_b) + dy * (dx * T2_b + dy * T2_c)
                vT2vT2_a, vT2vT2_b, vT2vT2_c = T2_a * vT2v, T2_b * vT2v, T2_c * vT2v

                # Internal terms: T * v * v^T * T
                t1vx = T1_a * dx + T1_b * dy
                t1vy = T1_b * dx + T1_c * dy
                T1vvT1_a, T1vvT1_b, T1vvT1_c = t1vx**2, t1vx*t1vy, t1vy**2

                t2vx = T2_a * dx + T2_b * dy
                t2vy = T2_b * dx + T2_c * dy
                T2vvT2_a, T2vvT2_b, T2vvT2_c = t2vx**2, t2vx*t2vy, t2vy**2

                # Assemble the H matrix
                s1 = l1 - l0
                H1_a = s1 * (T1_a - inv3 * (vT1vT1_a + 2.0 * T1vvT1_a))
                H1_b = s1 * (T1_b - inv3 * (vT1vT1_b + 2.0 * T1vvT1_b))
                H1_c = s1 * (T1_c - inv3 * (vT1vT1_c + 2.0 * T1vvT1_c))

                H2_a = l0 * (T2_a - inv4 * (vT2vT2_a + 2.0 * T2vvT2_a))
                H2_b = l0 * (T2_b - inv4 * (vT2vT2_b + 2.0 * T2vvT2_b))
                H2_c = l0 * (T2_c - inv4 * (vT2vT2_c + 2.0 * T2vvT2_c))

                H_a, H_b, H_c = H1_a + H2_a, H1_b + H2_b, H1_c + H2_c

                # Similarity Transform: V = c * (R * H * R^T)
                V_a = c * (Ra2 * H_a + 2.0 * Rab * H_b + Rb2 * H_c)
                V_b = c * (Rb2 * H_b + Rab * H_a + Rac * H_b + Rbc * H_c)
                V_c = c * (Rb2 * H_a + 2.0 * Rbc * H_b + Rc2 * H_c)

                # Direct Accumulation
                vf_slice = VF[x0:x0 + w, x1:x1 + w]
                vf_slice[:, :, 0, 0] += V_a
                vf_slice[:, :, 0, 1] += V_b
                vf_slice[:, :, 1, 0] += V_b
                vf_slice[:, :, 1, 1] += V_c
                
    return VF[pad:-pad, pad:-pad, :, :]


def tk_vote2(T, sigma1=3, sigma2=0, power=1, plate=True, N=10):
    """
    Executes tensor kernels voting over an entire input tensor field T.
    """
    # 1. Eigendecomposition and magnitude sorting
    evals, evecs = np.linalg.eigh(T)
    magValues = np.abs(evals)
    
    idx = np.argsort(magValues, -1)
    sorted_evals = np.take_along_axis(evals, idx, -1)
    sorted_mags = np.take_along_axis(magValues, idx, -1)
    
    E = np.zeros_like(evecs)
    E[..., 0, :] = np.take_along_axis(evecs[..., 0, :], idx, -1)
    E[..., 1, :] = np.take_along_axis(evecs[..., 1, :], idx, -1)

    # 2. Grid & Precomputations
    sigmax = max(sigma1, sigma2)
    pad = int(3 * sigmax)
    w = 2 * pad + 1
    x = np.linspace(-pad, pad, w)
    X0, X1 = np.meshgrid(x, x)
    
    L_sq = X0**2 + X1**2
    L = np.sqrt(L_sq)
    
    Dx = np.divide(X0, L, out=np.zeros_like(X0), where=L!=0)
    Dy = np.divide(X1, L, out=np.zeros_like(X1), where=L!=0)
    
    # Decay attenuations c1, c2 at every offset. Convention (VotingMath2.ipynb):
    #   exp(-ℓ²/σ²) = 0 if ℓ > 0 and σ = 0
    #   exp(-ℓ²/σ²) = 1 if ℓ = 0 (any σ, including σ = 0)
    # So when σ=0 we use zeros for ℓ>0, then overwrite the center to 1 below.
    d1 = np.exp(-L_sq / sigma1**2) if sigma1 > 0 else np.zeros_like(L)
    d2 = np.exp(-L_sq / sigma2**2) if sigma2 > 0 else np.zeros_like(L)
    d1[pad, pad] = 1.0
    d2[pad, pad] = 1.0

    # Normalization factor (eta)
    num = np.pi * math.factorial(2 * power)
    den = 2**(2 * power) * (math.factorial(power)**2)
    eta_val = den / (num * (sigma1**2 + sigma2**2))

    # Base rotation matrix components
    R11 = 1.0 - 2.0 * Dx**2
    R12 = -2.0 * Dx * Dy
    R22 = 1.0 - 2.0 * Dy**2

    # 3. Precompute Plate Field (PF) entirely outside the loop
    if plate:
        PF = np.zeros((w, w, 2, 2))
        if N == 0:
            # Analytical Plate Field
            c = 1.0 / (sigma1**2 + sigma2**2)
            ALPHA = np.arctan2(X1, X0)
            TWO_ALPHA = 2 * ALPHA
            COS_2A = np.cos(TWO_ALPHA)
            SIN_2A = np.sin(TWO_ALPHA)

            M00 = 0.25 * (COS_2A + 2)
            M01 = 0.25 * SIN_2A
            M11 = 0.25 * (2 - COS_2A)

            PF[:, :, 0, 0] = c * (d1 * (1 - M00) + d2 * M00)
            PF[:, :, 0, 1] = c * (d1 * (0 - M01) + d2 * M01)
            PF[:, :, 1, 0] = PF[:, :, 0, 1]
            PF[:, :, 1, 1] = c * (d1 * (1 - M11) + d2 * M11)
        else:
            # Numerical Plate Field
            dtheta = np.pi / N
            for n in range(N):
                nx = np.cos(n * dtheta)
                ny = np.sin(n * dtheta)

                nqTd = nx * Dx + ny * Dy
                ncos_2_theta = nqTd**2
                nsin_2_theta = 1.0 - ncos_2_theta
                
                DECAY = d1 * np.power(nsin_2_theta, power) + d2 * np.power(ncos_2_theta, power)
                nRq_x = R11 * nx + R12 * ny
                nRq_y = R12 * nx + R22 * ny

                V00 = eta_val * DECAY * (nRq_x * nRq_x)
                V01 = eta_val * DECAY * (nRq_x * nRq_y)
                V11 = eta_val * DECAY * (nRq_y * nRq_y)

                # At the center (L=0): R=I, qTd=0, so DECAY=d1=1 and (Rq)(Rq)^T = qq^T.
                # The general expressions above already produce eta_val * (nx*nx, nx*ny, ny*ny)
                # there since d1[pad,pad]=1 was forced above. No center-patch needed.

                # Riemann weight dβ = π/N for integral over [0, π]
                w_beta = np.pi / N
                PF[:, :, 0, 0] += V00 * w_beta
                PF[:, :, 0, 1] += V01 * w_beta
                PF[:, :, 1, 0] += V01 * w_beta
                PF[:, :, 1, 1] += V11 * w_beta

    # 4. Allocate Output and Accumulate Loop
    VF = np.pad(np.zeros(T.shape), ((pad, pad), (pad, pad), (0, 0), (0, 0)))

    for x0 in range(T.shape[0]):
        for x1 in range(T.shape[1]):
            l0_mag = sorted_mags[x0, x1, 0]
            l1_mag = sorted_mags[x0, x1, 1]
            l0 = sorted_evals[x0, x1, 0]
            l1 = sorted_evals[x0, x1, 1]

            # Optimization: Skip completely empty tensors
            if l0_mag == 0 and l1_mag == 0:
                continue

            qx = E[x0, x1, 0, 1]
            qy = E[x0, x1, 1, 1]
            vf_slice = VF[x0:x0 + w, x1:x1 + w]

            # --- Stick Vote Execution ---
            scale_stick = (l1_mag - l0_mag) * np.sign(l1)
            
            if scale_stick != 0:
                qTd = qx * Dx + qy * Dy
                cos_2_theta = qTd**2
                sin_2_theta = 1.0 - cos_2_theta
                
                DECAY = d1 * np.power(sin_2_theta, power) + d2 * np.power(cos_2_theta, power)
                
                Rq_x = R11 * qx + R12 * qy
                Rq_y = R12 * qx + R22 * qy
                
                S00 = scale_stick * eta_val * DECAY * (Rq_x * Rq_x)
                S01 = scale_stick * eta_val * DECAY * (Rq_x * Rq_y)
                S11 = scale_stick * eta_val * DECAY * (Rq_y * Rq_y)

                # At the center (L=0): R=I, qTd=0, DECAY=d1=1, (Rq)(Rq)^T=qq^T.
                # Since d1[pad,pad]=1 above, the general expression already gives the
                # correct self-vote: scale_stick * eta_val * (qx*qx, qx*qy, qy*qy).

                vf_slice[:, :, 0, 0] += S00
                vf_slice[:, :, 0, 1] += S01
                vf_slice[:, :, 1, 0] += S01
                vf_slice[:, :, 1, 1] += S11

            # --- Plate Vote Execution ---
            if plate:
                scale_plate = l0
                
                if scale_plate != 0:
                    vf_slice[:, :, 0, 0] += scale_plate * PF[:, :, 0, 0]
                    vf_slice[:, :, 0, 1] += scale_plate * PF[:, :, 0, 1]
                    vf_slice[:, :, 1, 0] += scale_plate * PF[:, :, 0, 1]
                    vf_slice[:, :, 1, 1] += scale_plate * PF[:, :, 1, 1]

    return VF[pad:-pad, pad:-pad, :, :]