import numpy as np
import matplotlib.pyplot as plt
import tensors as ts

def polar_diff(gt_T, pred_T, ax, title="", use_mask=False):
    """
    Plot eigenvector angular difference (theta) vs eigenvalue magnitude difference (rho)
    on a polar axis. No normalization applied.

    Args:
        gt_T:     Ground truth tensor field, shape (H, W, 2, 2)
        pred_T:   Predicted tensor field, shape (H, W, 2, 2)
        ax:       Polar matplotlib axis
        title:    Plot title
        use_mask: If True, mask background pixels where gt_T has no dominant structure
    """
    alpha = ts.metrics._scale_factor(gt_T, pred_T)
    pred_T_scaled = alpha * pred_T
    gt_sym   = 0.5 * (gt_T   + np.swapaxes(gt_T,   -1, -2))
    pred_sym = 0.5 * (pred_T_scaled + np.swapaxes(pred_T_scaled, -1, -2))

    gt_evals,   gt_evecs   = np.linalg.eigh(gt_sym)
    pred_evals, pred_evecs = np.linalg.eigh(pred_sym)

    l1_gt   = np.abs(gt_evals[..., 1])
    v1_gt   = gt_evecs[..., :, 1]
    l1_pred = np.abs(pred_evals[..., 1])
    v1_pred = pred_evecs[..., :, 1]

    dot   = np.clip(np.abs(np.sum(v1_gt * v1_pred, axis=-1)), 0.0, 1.0)
    theta = np.arccos(dot)
    rho   = np.abs(l1_gt - l1_pred)

    theta_flat = theta.flatten()
    rho_flat   = rho.flatten()

    if use_mask:
        gap  = np.abs(gt_evals[..., 1]) - np.abs(gt_evals[..., 0])
        mask = ((l1_gt > 1e-8) & (gap > 1e-6)).flatten()
        theta_flat = theta_flat[mask]
        rho_flat   = rho_flat[mask]

    scatter = ax.scatter(theta_flat, rho_flat, c=rho_flat, cmap='magma', alpha=0.3, s=5)

    max_rho  = float(np.max(rho_flat)) if len(rho_flat) > 0 else 1.0
    rmax_val = max(max_rho * 1.05, 1e-9)

    ax.set_rmax(rmax_val)
    ax.set_rticks(np.linspace(rmax_val / 4, rmax_val, 4))
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_title(title, va='bottom', fontsize=14, pad=15)
    ax.set_xlabel("Orientation Difference (Degrees)", labelpad=10)

    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.7)
    cbar.set_label("Eigenvalue Difference (ρ)")


def sif_error(gt_T, pred_T, ax, title="", use_mask=False):
    """
    Plot the SIF error map between gt_T and pred_T on ax.

    Args:
        gt_T:   Ground truth tensor field, shape (H, W, 2, 2)
        pred_T: Predicted tensor field, shape (H, W, 2, 2)
        ax:     Matplotlib axis
        title:  Plot title
    """
    err = ts.metrics.sif(gt_T, pred_T)
    im = ax.imshow(err)
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_fields(fields_dict, gt_T, plot_fn, title="", use_mask=False, subplot_kw=None, figsize_per_plot=(4, 4)):
    """
    Call plot_fn for each field in fields_dict against gt_T and show them in a row.

    Args:
        fields_dict:      Dict mapping name -> tensor field, shape (H, W, 2, 2)
        gt_T:             Ground truth tensor field, shape (H, W, 2, 2)
        plot_fn:          Callable with signature plot_fn(gt_T, pred_T, ax, title, use_mask)
        title:            Optional figure-level title shown at the top
        use_mask:         Passed through to plot_fn
        subplot_kw:       Extra kwargs for plt.subplots (e.g. {'projection': 'polar'})
        figsize_per_plot: (width, height) per subplot panel
    """
    n  = len(fields_dict)
    kw = subplot_kw or {}
    w, h = figsize_per_plot

    fig, axes = plt.subplots(1, n, subplot_kw=kw, figsize=(w * n, h))
    if n == 1:
        axes = [axes]

    for ax, (name, field) in zip(axes, fields_dict.items()):
        plot_fn(gt_T, field, ax, title=name, use_mask=use_mask)

    if title:
        fig.suptitle(title, fontsize=14)

    plt.tight_layout()
    plt.show(block=False)
