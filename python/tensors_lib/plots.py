import numpy as np
import matplotlib.pyplot as plt
import tensors as ts


def _polar_diff_compute(gt_T, pred_T, use_mask=False):
    """Compute the theta and rho values for the polar_diff plot."""
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

    return theta_flat, rho_flat


def polar_rho(gt_T, pred_T, use_mask=False):
    """Return rho (eigenvalue magnitude difference) values used in polar_diff, for range computation."""
    _, rho_flat = _polar_diff_compute(gt_T, pred_T, use_mask)
    return rho_flat


def polar_diff(gt_T, pred_T, ax, title="", use_mask=False, clim=None):
    """
    Plot eigenvector angular difference (theta) vs eigenvalue magnitude difference (rho)
    on a polar axis. No normalization applied.

    Args:
        gt_T:     Ground truth tensor field, shape (H, W, 2, 2)
        pred_T:   Predicted tensor field, shape (H, W, 2, 2)
        ax:       Polar matplotlib axis
        title:    Plot title
        use_mask: If True, mask background pixels where gt_T has no dominant structure
        clim:     (vmin, vmax) tuple; vmax is used as rmax for the polar axis and colorbar
    """
    theta_flat, rho_flat = _polar_diff_compute(gt_T, pred_T, use_mask)

    if clim is not None:
        rmax_val = clim[1]
    else:
        max_rho  = float(np.max(rho_flat)) if len(rho_flat) > 0 else 1.0
        rmax_val = max(max_rho * 1.05, 1e-9)

    scatter = ax.scatter(theta_flat, rho_flat, c=rho_flat, alpha=0.5, cmap='magma', s=5,
                         vmin=0, vmax=rmax_val)

    ax.set_rmax(rmax_val)
    ax.set_rticks(np.linspace(rmax_val / 4, rmax_val, 4))
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.set_title(title, va='bottom', fontsize=14, pad=15)
    ax.set_xlabel("Orientation Difference (Degrees)", labelpad=10)

    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.7)
    cbar.set_label("Eigenvalue Difference (ρ)")


def sif_error(gt_T, pred_T, ax, title="", use_mask=False, clim=None):
    """
    Plot the SIF error map between gt_T and pred_T on ax.

    Args:
        gt_T:   Ground truth tensor field, shape (H, W, 2, 2)
        pred_T: Predicted tensor field, shape (H, W, 2, 2)
        ax:     Matplotlib axis
        title:  Plot title
        clim:   (vmin, vmax) tuple for the colorbar range
    """
    err = ts.metrics.sif(gt_T, pred_T)
    vmin, vmax = clim if clim is not None else (None, None)
    im = ax.imshow(err, cmap='magma', vmin=vmin, vmax=vmax)
    ax.set_title(title, fontsize=10)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_fields(fields_dict, gt_T, plot_fn, title="", use_mask=False, subplot_kw=None, figsize_per_plot=(4, 4), data_fn=None):
    """
    Call plot_fn for each field in fields_dict against gt_T and show them in a row.

    Args:
        fields_dict:      Dict mapping name -> tensor field, shape (H, W, 2, 2)
        gt_T:             Ground truth tensor field, shape (H, W, 2, 2)
        plot_fn:          Callable with signature plot_fn(gt_T, pred_T, ax, title, use_mask, clim)
        title:            Optional figure-level title shown at the top
        use_mask:         Passed through to plot_fn
        subplot_kw:       Extra kwargs for plt.subplots (e.g. {'projection': 'polar'})
        figsize_per_plot: (width, height) per subplot panel
        data_fn:          Optional callable (gt_T, pred_T) -> array; if provided, the global
                          (-max_abs, +max_abs) range across all fields is computed and passed
                          as clim to every plot_fn call so all subplots share the same scale
    """
    n  = len(fields_dict)
    kw = subplot_kw or {}
    w, h = figsize_per_plot

    fig, axes = plt.subplots(1, n, subplot_kw=kw, figsize=(w * n, h))
    if n == 1:
        axes = [axes]

    clim = None
    if data_fn is not None:
        all_vals = np.concatenate([np.ravel(data_fn(gt_T, field)) for field in fields_dict.values()])
        max_abs  = float(np.max(np.abs(all_vals)))
        clim     = (-max_abs, max_abs)

    for ax, (name, field) in zip(axes, fields_dict.items()):
        plot_fn(gt_T, field, ax, title=name, use_mask=use_mask, clim=clim)

    if title:
        fig.suptitle(title, fontsize=14)

    plt.tight_layout()
    plt.show(block=False)
