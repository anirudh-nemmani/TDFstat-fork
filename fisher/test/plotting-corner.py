import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import corner

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.grid": True,
    "grid.alpha": 0.28,
    "axes.edgecolor": "0.15",
    "axes.linewidth": 1.3,
    "legend.frameon": False,
    "legend.fontsize": 16,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.minor.size": 5,
    "ytick.minor.size": 5,
    "xtick.minor.width": 1.0,
    "ytick.minor.width": 1.0,
})

rng = np.random.default_rng(9)


def fisher_samples(fisher_file, mean, n, rng):
    """Gaussian samples with covariance F^-1, drawn in normalised coordinates.

    Fisher matrix has a condition number of ~1e32, so decomposing F^-1 directly
    (as multivariate_normal does) loses the small variances (e.g. fdot).
    R = D^-1 F D^-1 with D = sqrt(diag F) has a condition number of ~1e10, which is safe.
    """
    F = np.loadtxt(fisher_file)
    d = np.sqrt(np.diag(F))
    R = F / np.outer(d, d)
    Rinv = np.linalg.inv(R)
    Rinv = 0.5 * (Rinv + Rinv.T)
    L = np.linalg.cholesky(Rinv)
    z = rng.standard_normal((n, len(d)))
    samples = mean + (z @ L.T) / d

    # sanity check: sample std must match the Fisher errors
    sigma = np.sqrt(np.diag(Rinv)) / d
    print("sample std  :", np.array2string(samples.std(axis=0), precision=3))
    print("Fisher sigma:", np.array2string(sigma, precision=3))
    return samples


mean_values = np.array([26.3, -8.5e-11, 3.86687548714555,
                        0.74835013347064, 2.984905359501588,
                        0.61465097870534, 0.11626181938654])

N = int(1e4)
samples = fisher_samples('./fisher_matrix_001_0017_analytical.txt', mean_values, N, rng)

plot_range = [(lo, hi) for lo, hi in zip(samples.min(axis=0), samples.max(axis=0))]

corner_plot_settings = {
    "bins": 100,
    "smooth": 0.5,
    "range": plot_range,
    "truths": mean_values,
    "truth_color": "k",
    "labels": [r'$f$ [Hz]', r'$\dot{f}$ [Hz/s]', r'$\alpha$ [rad]',
               r'$\delta$ [rad]', r'$\iota$ [rad]', r'$\psi$ [rad]', r'$\phi$ [rad]'],
    "title_fmt": ".2e",
    "title_kwargs": {"fontsize": 14},
    "label_kwargs": {"fontsize": 14},
    "quantiles": [0.16, 0.5, 0.84],
    "hist_kwargs": {"density": True, "alpha": 0.7},
    "fill_contours": True,
    "contourf_kwargs": {"alpha": [0.0, 0.2, 0.5, 0.8]},
    "plot_datapoints": False,
    "plot_density": False,
    "levels": [0.68, 0.95, 0.997],
    "use_math_text": True,
    "labelpad": 0,
}

fig = plt.figure(figsize=(18, 18))

corner.corner(samples, fig=fig, color="C1",
              show_titles=True, **corner_plot_settings)

handles = [Line2D([], [], color="C1", lw=3, label=r"Analytical Fisher")]
fig.legend(handles=handles, loc="upper right",
           bbox_to_anchor=(0.92, 0.92), fontsize=22)

def convert_sci_to_latex(s):
    def repl(match):
        coeff = match.group(1)
        exp = int(match.group(2))
        return rf"{coeff} \times 10^{{{exp}}}"
    return re.sub(r'(\d+\.\d+)e([+-]?\d+)', repl, s)

for ax in fig.axes:
    title = ax.get_title()
    if title:
        ax.set_title(convert_sci_to_latex(title), fontsize=12)

fig.subplots_adjust(wspace=0.15, hspace=0.15)
plt.show()
