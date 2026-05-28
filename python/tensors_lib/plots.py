import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import tensors as ts



def _polar_diff_compute(gt_T, pred_T):
    """Compute theta, rho, and signed eccentricity difference for the polar_diff plot."""

    # scale predicted tensors to match the energy of the GT
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

    # signed: + means gt more elongated, - means pred more elongated
    ecc_diff = ts._eccentricity(gt_evals) - ts._eccentricity(pred_evals)

    theta_flat    = theta.flatten()
    rho_flat      = rho.flatten()
    ecc_diff_flat = ecc_diff.flatten()

    return theta_flat, rho_flat, ecc_diff_flat


def polar_rho(gt_T, pred_T):
    """Return rho (eigenvalue magnitude difference) values used in polar_diff, for range computation."""
    _, rho_flat, _ = _polar_diff_compute(gt_T, pred_T)
    return rho_flat


def polar_ecc_diff(gt_T, pred_T):
    """Return eccentricity difference values used in polar_diff, for color range computation."""
    _, _, ecc_diff_flat = _polar_diff_compute(gt_T, pred_T)
    return ecc_diff_flat


def polar_diff(gt_T, pred_T, ax, title="", clim=None, color_clim=None, show_colorbar=True):
    """
    Plot eigenvector angular difference (theta) vs eigenvalue magnitude difference (rho)
    on a polar axis, colored by the signed eccentricity difference (gt − pred).

    Args:
        gt_T:       Ground truth tensor field, shape (H, W, 2, 2)
        pred_T:     Predicted tensor field, shape (H, W, 2, 2)
        ax:         Polar matplotlib axis
        title:      Plot title
        
        clim:       (vmin, vmax) tuple; vmax is used as rmax for the polar axis
        color_clim: Symmetric colorbar half-range for eccentricity difference; if None,
                    computed from this subplot's data
    """
    theta_flat, rho_flat, ecc_diff_flat = _polar_diff_compute(gt_T, pred_T)

    if clim is not None:
        rmax_val = clim[1]
    else:
        max_rho  = float(np.max(rho_flat)) if len(rho_flat) > 0 else 1.0
        rmax_val = max(max_rho * 1.05, 1e-9)

    if color_clim is not None:
        ecc_lim = color_clim
    else:
        ecc_lim = float(np.max(np.abs(ecc_diff_flat))) if len(ecc_diff_flat) > 0 else 1.0
        ecc_lim = max(ecc_lim, 1e-9)

    scatter = ax.scatter(theta_flat, rho_flat, c=ecc_diff_flat, alpha=0.5, cmap='managua', s=8, #edgecolors='black',
                         vmin=-ecc_lim, vmax=ecc_lim)

    ax.set_rlim(0, rmax_val)
    ax.xaxis.grid(True, color='#b0b0b0')
    ax.set_rticks(np.linspace(rmax_val / 4, rmax_val, 4))
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
    ax.tick_params(labelsize=5)
    ax.set_thetamin(0)
    ax.set_thetamax(90)
    ax.spines['polar'].set_edgecolor('#b0b0b0')
    ax.spines['start'].set_edgecolor('#b0b0b0')
    ax.spines['end'].set_edgecolor('#b0b0b0')
    if title:
        ax.set_title(title, va='bottom', fontsize=10, pad=20)
    ax.set_xlabel("Largest Eigenvalue Diff", labelpad=14, fontsize=8)

    if show_colorbar:
        cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.7)
        cbar.set_label("Eccentricity Diff (gt − pred)", fontsize=8)
        cbar.ax.tick_params(labelsize=7)


def sif_error(gt_T, pred_T, ax, title="", clim=None, show_colorbar=True, **_):
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
    if show_colorbar:
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def plot_fields(fields_dict, gt_T, plot_fn, title="", subplot_kw=None, figsize_per_plot=(4, 4), data_fn=None, color_data_fn=None, params_dict=None, subtitle_dict=None, axes=None, show_subplot_titles=True):
    """
    Call plot_fn for each field in fields_dict against gt_T and show them in a row.

    Args:
        fields_dict:      Dict mapping name -> tensor field, shape (H, W, 2, 2)
        gt_T:             Ground truth tensor field, shape (H, W, 2, 2)
        plot_fn:          Callable with signature plot_fn(gt_T, pred_T, ax, title, clim, color_clim)
        title:            Optional figure-level title shown at the top
        subplot_kw:       Extra kwargs for plt.subplots (e.g. {'projection': 'polar'})
        figsize_per_plot: (width, height) per subplot panel; ignored when axes is provided
        data_fn:          Optional callable (gt_T, pred_T) -> array; if provided, the global
                          (-max_abs, +max_abs) range across all fields is computed and passed
                          as clim to every plot_fn call so all subplots share the same scale
        color_data_fn:    Optional callable (gt_T, pred_T) -> array; if provided, the global
                          symmetric color range is computed and passed as color_clim to every
                          plot_fn call so all colorbars share the same scale
        params_dict:      Optional dict mapping field name -> param string shown as a small
                          subtitle beneath each subplot
        axes:             Optional list of pre-created matplotlib axes, one per field.
                          When provided, no figure is created and plt.show() is not called,
                          allowing the caller to embed this row inside a larger figure.
    """
    n  = len(fields_dict)

    # create a new figure if axes are not provided
    if axes is not None:
        fig      = axes[0].get_figure()
        show_fig = False
    else:
        kw        = subplot_kw or {}
        w, h      = figsize_per_plot
        fig, axes = plt.subplots(1, n, subplot_kw=kw, figsize=(w * n, h))
        if n == 1:
            axes = [axes]
        show_fig  = True

    # 
    clim = None
    if data_fn is not None:
        all_vals = np.concatenate([np.ravel(data_fn(gt_T, field)) for field in fields_dict.values()])
        max_abs  = float(np.max(np.abs(all_vals)))
        clim     = (-max_abs, max_abs)

    color_clim = None
    if color_data_fn is not None:
        all_color = np.concatenate([np.ravel(color_data_fn(gt_T, field)) for field in fields_dict.values()])
        color_clim = max(float(np.max(np.abs(all_color))), 1e-9)

    for i, (ax, (name, field)) in enumerate(zip(axes, fields_dict.items())):
        is_last = (i == n - 1)
        subplot_title = name if show_subplot_titles else ""
        plot_fn(gt_T, field, ax, title=subplot_title, clim=clim, color_clim=color_clim, show_colorbar=is_last)
        if subtitle_dict and name in subtitle_dict:
            ax.set_title(ax.get_title() + "\n" + subtitle_dict[name], fontsize=10)
        if params_dict and name in params_dict and params_dict[name]:
            ax.text(0.5, -0.22, params_dict[name], transform=ax.transAxes,
                    ha='center', va='top', fontsize=10, color='dimgray')

    if title:
        if show_fig:
            fig.suptitle(title, fontsize=14)
        else:
            axes[0].set_ylabel(title, fontsize=11, labelpad=10)

    if show_fig:
        plt.tight_layout(rect=[0, 0.05, 1, 1])
        plt.show(block=False)


def tensor_orientation(gt_T, pred_T, ax, title="", clim=None, color_clim=None, **_):
    """
    Visualize the dominant eigenvector orientation of pred_T on ax.
    Hue encodes angle in [0, pi], saturation is scaled by eccentricity.
    """
    evals, evecs = np.linalg.eigh(pred_T)

    # angle of the dominant (largest) eigenvector
    theta = np.arctan2(evecs[..., 1, 1], evecs[..., 0, 1])

    # fold (-pi, 0) into (0, pi) to handle sign ambiguity of eigenvectors
    neg = theta < 0
    theta[neg] = np.pi - np.abs(theta[neg])
    theta /= np.pi                              # normalize to [0, 1] for HSV

    C = matplotlib.colormaps["hsv"](theta)         # RGBA in [0,1]
    zero_mask = evals[..., 1] < 1e-12              # largest eigenvalue ~ 0 -> zero tensor
    C[zero_mask] = [1, 1, 1, 1]                    # white

    #ecc = ts._eccentricity(evals)[..., np.newaxis]
    #C = ecc * C + (1 - ecc)                    # desaturate isotropic regions toward white

    im = ax.imshow(np.clip(C[..., :3], 0, 1))

    if title:
        ax.set_title(title, fontsize=10)
