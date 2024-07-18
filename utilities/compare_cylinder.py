import matplotlib.pyplot as plt
import numpy as np
import h5py
import treams
import treams.io
import matplotlib
import sys
from matplotlib.colors import Normalize, LogNorm
from validate import validate_hdf5_file
from validate import metric

if len(sys.argv) != 2:
    print("Usage: python script.py <path_to_hdf5_file_for_reference_cylinder>")
    sys.exit(1)


ref_path = (
    "../example_data/cylinder_tio2.tmat.h5"  # path to reference file computed with jcm
)
test_path = sys.argv[1]


def plotting(tmat, tmat2, wl, name1, name2):
    xs_1 = tmat.xs_ext_avg
    xs_2 = tmat2.xs_ext_avg
    fig, ax = plt.subplots(2, 3, figsize=(30, 20))
    len = 30
    tr = tmat[:len, :len].real
    ti = tmat[:len, :len].imag
    tr2 = tmat2[:len, :len].real
    ti2 = tmat2[:len, :len].imag
    mmax = np.max(
        [np.max(abs(tr2)), np.max(abs(tr)), np.max(abs(ti)), np.max(np.abs(ti2))]
    )
    norm = Normalize(vmin=-mmax, vmax=mmax)
    cmap = matplotlib.colormaps["coolwarm"]
    ax[0, 0].set_title(f"Real T {name1}")
    ax[1, 0].set_title(f"Imag T {name1}")
    ac1 = ax[0, 0].imshow(tr, norm=norm, cmap=cmap)  #
    ax[0, 0].axis("off")
    ac = ax[1, 0].imshow(ti, norm=norm, cmap=cmap)  # vmin=mmin, vmax=mmax,
    size = 0.45
    ax[1, 0].axis("off")
    cmap = matplotlib.colormaps["coolwarm"]
    ax[0, 1].set_title(f"Real T {name2}")
    ax[1, 1].set_title(f"Imag T {name2}")
    ac1 = ax[0, 1].imshow(tr2, norm=norm, cmap=cmap)  #

    ac = ax[1, 1].imshow(ti2, norm=norm, cmap=cmap)  # vmin=mmin, vmax=mmax,
    size = 0.45
    ax[0, 1].axis("off")
    ax[1, 1].axis("off")
    ax[0, 2].set_title("Real T Diff")
    ax[1, 2].set_title("Imag T Diff")
    ac1 = ax[0, 2].imshow(tr2 - tr, norm=norm, cmap=cmap)  #
    maxdiff = max([(tr2 - tr).max(), (ti2 - ti).max()])
    argmx = np.argmax(np.array([np.abs(tr2 - tr).max(), np.abs(ti2 - ti).max()]))
    if argmx:
        td = np.abs(ti2 - ti)
    else:
        td = np.abs(tr2 - tr)
    argdiff = np.argwhere(np.array(td) == np.max(np.array(td)))
    print(
        f"At {wl_max} nm the maximum absolute difference between corresponding T-matrix elements is {maxdiff} between elements at {argdiff}. Note that T-matrices have initial shapes {tmat.shape} and {tmat2.shape}"
    )
    ax[0, 2].axis("off")
    ac = ax[1, 2].imshow(ti2 - ti, norm=norm, cmap=cmap)  # vmin=mmin, vmax=mmax,
    size = 0.45
    ax[1, 2].axis("off")
    cb = plt.colorbar(ac, ax=ax.ravel().tolist(), fraction=0.047 * size)
    plt.suptitle(
        f"At {wl:.4f} nm: \n{name1}: $<\sigma_{{ext}}>$"
        + f"= {xs_1*1e-6:.4f} $\mu m^2$, {name2}: $<\sigma_{{ext}}>$"
        + f"= {xs_2*1e-6:.4f} $\mu m^2$"
    )
    figpath = "tmat_comparison.png"
    plt.savefig(figpath)
    print(f"Saved plot in {figpath}")


ref_tm = treams.io.load_hdf5(ref_path)
validate_hdf5_file(test_path)
ref_t00 = np.array([t[0, 0] for t in ref_tm])
ref_t11 = np.array([t[1, 1] for t in ref_tm])
ref_t10 = np.array([t[1, 0] for t in ref_tm])
ref_t01 = np.array([t[0, 1] for t in ref_tm])
ref_lambdas = np.array([2 * np.pi / tm.k0 for tm in ref_tm])
ref_xs = np.array([tm.xs_ext_avg for tm in ref_tm])
test_tm = treams.io.load_hdf5(test_path)
test_t00 = np.array([t[0, 0] for t in test_tm])
test_t11 = np.array([t[1, 1] for t in test_tm])
test_t10 = np.array([t[1, 0] for t in test_tm])
test_t01 = np.array([t[0, 1] for t in test_tm])
test_lambdas = np.array([2 * np.pi / tm.k0 for tm in test_tm])
test_xs = np.array([tm.xs_ext_avg for tm in test_tm])
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 10))
ax1.scatter(test_lambdas, test_xs, label="JCM")
ax1.plot(ref_lambdas, ref_xs, label="New method")
ax1.legend()
elements = ["$T_{00}$", "$T_{11}$", "$T_{10}$", "$T_{01}$"]
for i, t in enumerate(
    np.array(
        [
            (ref_t00, test_t00),
            (ref_t11, test_t11),
            (ref_t01, test_t01),
            (ref_t10, test_t10),
        ]
    )
):
    ax2.plot(ref_lambdas, t[0].real, label=f"JCM Real({elements[i]})")
    ax2.plot(ref_lambdas, t[0].imag, label=f"JCM Imag({elements[i]})")
    ax2.scatter(test_lambdas, t[1].real, label=f"New method: Real({elements[i]})")
    ax2.scatter(test_lambdas, t[1].imag, label=f"New method: Imag({elements[i]})")
ax2.set_ylabel("T-matrix elements")
ax2.set_xlabel("$\lambda$ (nm)")
plt.legend(fontsize=15)
ax1.set_ylabel("<$\sigma_{{ext}}$> (nm$^2$)")
ax1.set_xlabel("$\lambda$ (nm)")
figpath0 = "cross-section-tmatrix-2-2.png"
plt.savefig(figpath0)
print(f"Saved plots in {figpath0}")
common = np.isclose(test_lambdas[:, None], ref_lambdas)
if common.shape == (0,):
    print("The T-matrices were not computed at any coinciding wavelengths")
else:
    test_compare = test_xs[np.where(common)[1]]
    ref_compare = ref_xs[np.where(common)[0]]
    metric_ext = (ref_compare - test_compare) / ref_compare
    wl_compare = test_lambdas[np.where(common)[0]]
    wl_max = wl_compare[np.argmax(metric_ext)]
    print(
        f"Comparing T-matrices for reference TiO2 cylinder with radius 250 nm and height 300 nm: {len(wl_compare)} coinciding wavelengths were found. The relative difference between extinction cross-sections at these wavelengths has a maximum value of {metric_ext.max()}  at {wl_max:.4f} nm (value less than 0.01 is desired)."
    )
    metric2 = np.array([metric(ref_tm[i], test_tm[i]) for i in range(len(test_tm))])
    wl_max = wl_compare[np.argmax(metric2)]
    print(
        f"Using metric: half of the ratio of squared norm of the differences to the sum of squared norms of the T-matrices. Maximum value is {metric2.max()} (min: 0 (indentical), max: 1 (A = -B)) at {wl_max:.4f} nm."
    )
    tmat1 = ref_tm[np.argmax(metric2)]
    tmat2 = test_tm[np.argmax(metric2)]
    s1 = max([tmat1.shape[-2], tmat2.shape[-2]])
    s2 = max([tmat1.shape[-1], tmat2.shape[-1]])

    plotting(tmat1, tmat2, wl_max, "JCM", "New method")
