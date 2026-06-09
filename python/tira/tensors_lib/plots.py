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

    scatter = ax.scatter(theta_flat, rho_flat, c=ecc_diff_flat, alpha=0.3, cmap='managua', s=5,
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


def tensor_orientation(gt_T, pred_T, ax, title="", clim=None, color_clim=None, minor=False, **_):
    """
    Visualize an eigenvector orientation of pred_T on ax, matching the tvote2
    interactive tool's Vector Display (default settings, Magnitude = Lighten).

    Hue encodes the selected eigenvector angle via the cyclic HSV colormap. The
    color is then blended toward white by two sequential "Lighten" steps -- by
    eccentricity, then by magnitude -- so the plot shows orientation only where
    the tensor is both stick-like AND strong. This is the tvote2 recipe:
    isotropic regions (low eccentricity) and low-energy background (low magnitude)
    both wash to white, leaving only the salient structure colored.

    minor=False -> dominant (largest-eigenvalue) eigenvector.
    minor=True  -> minor (smallest-eigenvalue) eigenvector.
    """
    # Degenerate vote fields (e.g. sigma collapsing to an empty grid) produce a
    # zero-size tensor array. Render blank instead of crashing on l1.max().
    if pred_T.size == 0 or pred_T.shape[0] == 0 or pred_T.shape[1] == 0:
        ax.imshow(np.ones((1, 1, 3)))
        if title:
            ax.set_title(title, fontsize=10)
        return

    evals, evecs = np.linalg.eigh(pred_T)         # ascending: [...,0]=minor, [...,1]=major

    # angle of the selected eigenvector
    col = 0 if minor else 1
    theta = np.arctan2(evecs[..., 1, col], evecs[..., 0, col])

    # fold (-pi, 0) into (0, pi) to handle sign ambiguity of eigenvectors
    neg = theta < 0
    theta[neg] = np.pi - np.abs(theta[neg])
    theta /= np.pi                              # normalize to [0, 1] for HSV

    C = matplotlib.colormaps["hsv"](theta)         # RGBA in [0,1]
    zero_mask = evals[..., 1] < 1e-12              # largest eigenvalue ~ 0 -> zero tensor
    C[zero_mask] = [1, 1, 1, 1]                    # white

    # tvote2 Lighten blends: color toward white by (1 - weight) for eccentricity
    # then magnitude, applied as two separate passes (not a single combined
    # weight). Magnitude uses |lambda1| normalized linearly by its global max --
    # tvote2's threshold -- not a percentile/gamma stretch, so genuinely weak
    # background washes out instead of being boosted into saturated hue noise.
    ecc = ts._eccentricity(evals)[..., np.newaxis]
    l1 = np.abs(evals[..., 1])
    l1_max = float(l1.max())
    mag_norm = np.clip(l1 / l1_max, 0.0, 1.0) if l1_max > 0 else np.zeros_like(l1)
    mag_norm = mag_norm[..., np.newaxis]

    rgb = C[..., :3]
    rgb = rgb * ecc + (1.0 - ecc) * 1.0            # eccentricity Lighten
    rgb = rgb * mag_norm + (1.0 - mag_norm) * 1.0  # magnitude Lighten
    C[..., :3] = rgb

    ax.imshow(np.clip(C[..., :3], 0, 1))

    if title:
        ax.set_title(title, fontsize=10)


def tensor_saliency(T, ax, title="", log=True, cmap="magma"):
    """
    Visualize tensor saliency = largest |eigenvalue| of T on ax, i.e. how
    strong the local response is. log=True applies log1p for dynamic range.
    """
    evals, _ = ts.eigmag(T)
    sal = np.abs(evals[..., 1])
    ax.imshow(np.log1p(sal) if log else sal, cmap=cmap)
    if title:
        ax.set_title(title, fontsize=10)


def tensor_eccentricity(T, ax, title="", cmap="viridis"):
    """
    Visualize tensor eccentricity of T on ax: how stick-like vs isotropic the
    local tensor is (1 = perfect line, 0 = blob).
    """
    evals, _ = ts.eigmag(T)
    ax.imshow(ts._eccentricity(evals), cmap=cmap, vmin=0.0, vmax=1.0)
    if title:
        ax.set_title(title, fontsize=10)
