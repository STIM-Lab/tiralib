#%%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

import tensors as ts
import tensors_lib.plots as tsplot
from tensors_lib.cuda_wrappers import tk_vote2_cuda, CUDA_TK_AVAILABLE

N = 20

# Ground truth cross field and a degraded version with 30% dropout.
G = ts.phantoms.cross((N, N), np.pi / 2, 0, 4)
I = ts.phantoms.dropout(G, 0.3)

# Initial parameters.
s1_0, s2_0, p_0 = 5.0, 0.0, 1

# CUDA gives a responsive real-time slider; CPU fallback is slower.
def run_vote(s1, s2, power):
    if CUDA_TK_AVAILABLE:
        return tk_vote2_cuda(I, sigma1=s1, sigma2=s2, power=power, plate=True, N=20)
    return ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=power, plate=True, N=20)


fig = plt.figure(figsize=(6, 8))
ax_ori = fig.add_axes([0.1, 0.40, 0.8, 0.52])

# Slider axes.
ax_s1 = fig.add_axes([0.25, 0.27, 0.6, 0.03])
ax_s2 = fig.add_axes([0.25, 0.22, 0.6, 0.03])
ax_p  = fig.add_axes([0.25, 0.17, 0.6, 0.03])

s_s1 = Slider(ax_s1, r"$\sigma_1$", 0.1, 30.0, valinit=s1_0, valstep=0.1)
s_s2 = Slider(ax_s2, r"$\sigma_2$", 0.0, 30.0, valinit=s2_0, valstep=0.1)
s_p  = Slider(ax_p,  "power", 1, 30, valinit=p_0, valstep=1)

# SIMF error readout under the sliders.
simf_text = fig.text(0.5, 0.07, "", ha="center", va="center", fontsize=12)

backend = "CUDA" if CUDA_TK_AVAILABLE else "CPU"
fig.suptitle(f"TK Voted Cross Field — Tensor Orientation ({backend})", fontsize=12)


def update(_=None):
    s1, s2, power = s_s1.val, s_s2.val, int(s_p.val)
    TK = run_vote(s1, s2, power)

    ax_ori.clear()
    tsplot.tensor_orientation(G, TK, ax_ori)
    ax_ori.set_xticks([])
    ax_ori.set_yticks([])

    simf = ts.metrics.simf(G, TK)
    simf_text.set_text(
        f"SIMF error (vs ground truth cross):  {simf:.4f}\n"
        rf"$\sigma_1$={s1:.2f},  $\sigma_2$={s2:.2f},  power={power}")
    fig.canvas.draw_idle()


s_s1.on_changed(update)
s_s2.on_changed(update)
s_p.on_changed(update)

update()
plt.show()
