import matplotlib.animation as animation
import matplotlib.pyplot as plt
import meep as mp
import numpy as np
from matplotlib import ticker
from matplotlib.colors import CenteredNorm
from matplotlib.figure import Figure
from matplotlib.pyplot import Axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.typing import NDArray


def add_color_bar(fig: Figure, ax: Axes, im, ticks=None) -> None:
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    c = fig.colorbar(im, cax=cax, orientation="vertical", ticks=ticks)
    if ticks is None:
        c.ax.locator_params(nbins=3, tight=True)


def plot_im(im_data: NDArray, title="", centered=True):
    if mp.am_really_master():
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        add_im_to_plot(fig=fig, ax=ax, im_data=im_data, title=title, centered=centered)


def add_im_to_plot(
    fig: Figure, ax: Axes, im_data: NDArray, title="", centered=True, cbar_ticks=None
) -> None:
    ax.set_title(title)
    if centered and np.abs(im_data).sum() > 0:
        im = ax.imshow(im_data, cmap="RdBu", norm=CenteredNorm())
    else:
        im = ax.imshow(im_data, cmap="Purples")
    add_color_bar(fig, ax, im, cbar_ticks)
    ax.set_yticklabels([])
    ax.set_xticklabels([])


def plot_anim(data, path=None, title="", cmap="RdBu", interval=50) -> None:
    if mp.am_really_master():
        fig, ax = plt.subplots(1, 1, figsize=(6, 5))
        ax.set_title(title)
        im = ax.imshow(data[0], cmap=cmap)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        ax.set_aspect(1)
        plt.colorbar(im, cax=cax)
        ims = []
        for ii in range(data.shape[0]):
            im = ax.imshow(data[ii], cmap=cmap)
            ims.append([im])

        ani = animation.ArtistAnimation(fig, ims, interval=interval, repeat_delay=0)
        if path:
            ani.save(
                path,
                writer="imagemagick",
            )


def plot_complex_versus(
    data_1: NDArray,
    data_2: NDArray,
    title_1="data_1",
    title_2="data_2",
) -> None:
    def real(str: str):
        return f"$\Re (${str}$)$"

    def imag(str: str):
        return f"$\Im (${str}$)$"

    def abs(str: str):
        return f"Abs$(${str}$)$"

    if mp.am_really_master():
        fig, axes = plt.subplots(3, 3, figsize=(0.7 * 12, 0.7 * 9))
        add_im_to_plot(fig, axes[0, 0], np.real(data_1), real(title_1))
        add_im_to_plot(fig, axes[1, 0], np.real(data_2), real(title_2))
        add_im_to_plot(
            fig,
            axes[2, 0],
            np.real(data_1) - np.real(data_2),
            real(title_1) + " - " + real(title_2),
        )
        add_im_to_plot(fig, axes[0, 1], np.imag(data_1), imag(title_1))
        add_im_to_plot(fig, axes[1, 1], np.imag(data_2), imag(title_2))
        add_im_to_plot(
            fig,
            axes[2, 1],
            np.imag(data_1) - np.imag(data_2),
            imag(title_1) + " - " + imag(title_2),
        )
        add_im_to_plot(fig, axes[0, 2], np.abs(data_1), abs(title_1), centered=False)
        add_im_to_plot(fig, axes[1, 2], np.abs(data_2), abs(title_2), centered=False)
        add_im_to_plot(
            fig,
            axes[2, 2],
            np.abs(data_1) - np.abs(data_2),
            abs(title_1) + " - " + abs(title_2),
        )
        # non_zeros = np.abs(data_1) > 0.1
        # data_1 = np.asarray(data_1)[non_zeros]
        # data_2 = np.asarray(data_2)[non_zeros]
        # print("Relative differences of abs:")
        # print(np.max(np.abs(np.abs(data_1) - np.abs(data_2)) / np.abs(np.abs(data_1))))


def plot_complex(data: NDArray, title: str = "") -> None:
    if mp.am_really_master():
        fig, axes = plt.subplots(1, 3, figsize=(8, 4))
        add_im_to_plot(fig, axes[0], np.real(data), title + " real")
        add_im_to_plot(fig, axes[1], np.imag(data), title + " imag")
        add_im_to_plot(fig, axes[2], np.abs(data), title + " abs", centered=False)
