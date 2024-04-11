import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import treams
import treams.io

mpl.rcParams["lines.linewidth"] = 2.5
mpl.rcParams["font.size"] = 14

cyl = treams.io.load_hdf5("cylinder_tio2.h5")

positions = [[-300, 0, 0], [300, 0, 0]]
cyl_cluster = [
    treams.TMatrix.cluster([tm, tm], positions) for tm in cyl
]
cyl_cluster = [
    tm.interaction.solve().expand(
        treams.SphericalWaveBasis.default(4)
    )
    for tm in cyl_cluster
]

with h5py.File("two_cylinders.h5", "w") as fobj:
    treams.io.save_hdf5(fobj, cyl_cluster)

lmax = max(cyl[0].basis.l)
cyl_by_lmax = [
    [tm[treams.SphericalWaveBasis.default(i)] for tm in cyl]
    for i in range(1, lmax + 1)
]
cyl_cluster_by_lmax = [
    [tm[treams.SphericalWaveBasis.default(i)] for tm in cyl_cluster]
    for i in range(1, lmax + 1)
]

fig, axs = plt.subplots(1, 2, figsize=(14, 6))

freqs = np.array([tm.k0 * 299792.458 / (2 * np.pi) for tm in cyl])
lambdas = np.array([2 * np.pi / tm.k0 for tm in cyl])
linestyles = [":","--","-.", "-"]
for i in range(lmax):
    axs[0].plot(
        freqs,
        [tm.xs_ext_avg / 1e6 for tm in cyl_by_lmax[i]],
        linestyle=linestyles[i % lmax],
        zorder=-i,
    )

axs[0].set_xlabel("Frequency (THz)")
axs[0].set_ylabel("Cross-section (µm$^2$)")
axs[0].set_xlim([240, 400])

ax_top = axs[0].twiny()
ax_top.set_xlabel("Wavelength (nm)")
ax_top.set_xlim(lambdas[0], lambdas[-1])

axs[0].legend([f"$l_\\mathrm{{max}} = {i}$" for i in range(1, lmax + 1)])

for i in range(lmax):
    axs[1].plot(
        freqs,
        [tm.xs_ext_avg / 1e6 for tm in cyl_cluster_by_lmax[i]],
        linestyle=linestyles[i % lmax],
        zorder=-i,
    )

axs[1].set_xlabel("Frequency (THz)")
axs[1].set_ylabel("Cross-section (µm$^2$)")
axs[1].set_xlim([240, 400])

ax_top = axs[1].twiny()
ax_top.set_xlabel("Wavelength (nm)")
ax_top.set_xlim(lambdas[0], lambdas[-1])

fig.savefig("xs_cyl_tio2.png")