import scipy.ndimage
from scipy.optimize import minimize_scalar, differential_evolution, dual_annealing
from . import metrics
import tensors as ts

# Minimizes SIMF(G, blur(I, σ)) over scalar σ using Brent's method on the bracketed interval.
# The unimodal assumption lets it find the minimum with O(log(1/ε)) evaluations via golden-section
# search narrowing the bracket until the tolerance is met.
def optimize_blur(G, I, bounds=(0.1, 50.0)):
    def objective(sg):
        return metrics.simf(G, scipy.ndimage.gaussian_filter(I, sigma=(sg, sg, 0, 0), mode="constant"))

    res = minimize_scalar(objective, bounds=bounds, method='bounded')
    B_opt = scipy.ndimage.gaussian_filter(I, sigma=(res.x, res.x, 0, 0), mode="constant")
    return res.x, B_opt, res.fun


# Minimizes SIMF(G, ATV(I, σ)) over scalar σ using the same bounded Brent's method as optimize_blur.
def optimize_atv(G, I, bounds=(0.5, 20.0)):
    def objective(sigma):
        return metrics.simf(G, ts.atv_vote2(I, sigma=sigma))

    res = minimize_scalar(objective, bounds=bounds, method='bounded')
    return res.x, ts.atv_vote2(I, sigma=res.x), res.fun


# Minimizes SIMF(G, TK(I, σ1, σ2, power, plate)) over three continuous and one binary parameter
# using differential evolution.
def optimize_tk(G, I, bounds=None, seed=42):
    if bounds is None:
        bounds = [(1.0, 30.0), (0.0, 20.0), (1, 30), (0, 1)]

    def objective(params):
        s1, s2, pw = float(params[0]), float(params[1]), max(1, int(round(params[2])))
        plate = bool(round(params[3]))
        return metrics.simf(G, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=plate, N=20))

    res = differential_evolution(objective, bounds=bounds, seed=seed,
                                 maxiter=150, popsize=5, tol=1e-3, polish=False)
    s1, s2, pw, plate = res.x[0], res.x[1], max(1, int(round(res.x[2]))), bool(round(res.x[3]))
    return s1, s2, pw, plate, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=plate, N=20), res.fun


# Minimizes the same TK objective using dual annealing — a hybrid of fast simulated annealing
# and classical annealing.
def optimize_tk_annealing(G, I, bounds=None, seed=42):
    if bounds is None:
        bounds = [(1.0, 30.0), (0.0, 20.0), (1, 30), (0, 1)]

    def objective(params):
        s1, s2, pw = float(params[0]), float(params[1]), max(1, int(round(params[2])))
        plate = bool(round(params[3]))
        return metrics.simf(G, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=plate, N=20))

    res = dual_annealing(objective, bounds=bounds, seed=seed, maxiter=150)
    s1, s2, pw, plate = res.x[0], res.x[1], max(1, int(round(res.x[2]))), bool(round(res.x[3]))
    return s1, s2, pw, plate, ts.tk_vote2(I, sigma1=s1, sigma2=s2, power=pw, plate=plate, N=20), res.fun
