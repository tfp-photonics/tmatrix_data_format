import h5py
import numpy as np
import scipy
import matplotlib.pyplot as plt
import argparse
import os
import sys
import re


def metric(t, tnew):
    return (
        0.5
        * np.sum(np.abs(t - tnew) ** 2, axis=(-2, -1))
        / np.sum(np.abs(t) ** 2 + np.abs(tnew) ** 2, axis=(-2, -1))
    )


def visitfunc(f, x):
    if isinstance(f[x], h5py.Dataset):
        print(x, "Dataset:", f[x][...].shape)
    if isinstance(f[x], h5py.Group):
        print("Attr", f[x].attrs.items())
        pass


def czinf(tmat, js, ms, taus):
    """
    Check rotational symmetry (assume z is the axis of rotation)
    tmats, js, ms, taus
    """
    tol = 0
    ps = np.zeros(len(taus), dtype=int)
    ps[taus == taus[0]] = 1
    ps[taus != taus[0]] = 0
    m_out = ms + np.zeros_like(ms[:, None], dtype=int)
    m_in = ms[:, None] + np.zeros_like(ms)
    p_out = ps[None, :] + np.zeros_like(ps[:, None], dtype=int)
    p_in = ps[:, None] + np.zeros_like(ps, dtype=int)
    l_in = js[:, None] + np.zeros_like(js)
    l_out = js + np.zeros_like(js[:, None])

    tnew = np.zeros_like(tmat)
    for i in range(tmat.shape[0]):
        for j in range(tmat.shape[-1]):
            ind = (
                (m_in == -m_in[i, j])
                & (m_out == -m_out[i, j])
                & (l_in == l_in[i, j])
                & (l_out == l_out[i, j])
                & (p_in == p_in[i, j])
                & (p_out == p_out[i, j])
            )
            tnew[i, j] = tmat[
                (m_in == -m_in[i, j])
                & (m_out == -m_out[i, j])
                & (l_in == l_in[i, j])
                & (l_out == l_out[i, j])
                & (p_in == p_in[i, j])
                & (p_out == p_out[i, j])
            ]
    tnew = tnew * (-1.0) ** (p_in + p_out)
    tnew[m_in != m_out] = 0.0
    tol = metric(tnew, tmat)
    return tol


def passive(T):
    A = -T - np.conj(T).T - 2 * np.conj(T).T @ T
    Hc = np.conj(A).T
    E = np.linalg.eigvalsh(A)
    i = 1
    while np.abs(np.min(E) * 10**i) <= 1:
        i += 1
    bound = round(np.min(E), 1)
    return 10**-i


def mirrorxy(tmat, js, ms, taus):
    ts = np.zeros(len(taus), dtype=int)
    ts[taus == b"electric"] = -1
    ts[taus == b"magnetic"] = 1  # ?
    ms_out = ms + np.zeros_like(ms[:, None], dtype=int)
    ms_in = ms[:, None] + np.zeros_like(ms)
    t_out = ts[None, :] + np.zeros_like(ts[:, None], dtype=int)
    t_in = ts[:, None] + np.zeros_like(ts, dtype=int)
    js_in = js[:, None] + np.zeros_like(js)
    js_out = js + np.zeros_like(js[:, None])
    factor = (-1.0) ** (js_in + js_out + ms_in + ms_out) * t_in * t_out
    tnew = factor * tmat
    tol = metric(tnew, tmat)
    return tol


def pmirrorxy(tmat, js, ms, taus):
    l_in = js[:, None] + np.zeros_like(js)
    l_out = js + np.zeros_like(js[:, None])
    m_out = ms + np.zeros_like(ms[:, None], dtype=int)
    m_in = ms[:, None] + np.zeros_like(ms)
    factor = (-1.0) ** (l_in + l_out + m_in + m_out)
    return pmirrorxyxz(tmat, js, ms, taus, factor)


def pmirrorxz(tmat, js, ms, taus):
    m_out = ms + np.zeros_like(ms[:, None], dtype=int)
    m_in = ms[:, None] + np.zeros_like(ms)
    factor = (-1.0) ** (m_in + m_out)
    return pmirrorxyxz(tmat, js, ms, taus, factor)


def pmirrorxyxz(tmat, js, ms, taus, factor):
    tol = 0
    N_J = js.max()
    ts = np.zeros(len(taus), dtype=int)
    ts[taus == b"positive"] = -1
    ts[taus == b"negative"] = 1  # ??????????? IS IT SO?
    m_out = ms + np.zeros_like(ms[:, None], dtype=int)
    m_in = ms[:, None] + np.zeros_like(ms)
    p_out = ts[None, :] + np.zeros_like(ts[:, None], dtype=int)
    p_in = ts[:, None] + np.zeros_like(ts, dtype=int)
    l_in = js[:, None] + np.zeros_like(js)
    l_out = js + np.zeros_like(js[:, None])
    n_range = 2 * (2 * np.arange(1, N_J) + 1)
    Nrng = np.concatenate(([0], n_range, [len(js) + 1]))
    tnew_sum = 0
    for i, ni in enumerate(Nrng[:-1]):
        for j, nj in enumerate(Nrng[:-1]):
            tsubmat = tmat[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            p_inc = p_in[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            p_outc = p_out[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            m_inc = m_in[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            m_outc = m_out[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            t_new = np.zeros_like(tsubmat, dtype=complex)
            t_new[::2, ::2] = tsubmat[1::2, 1::2]
            t_new[1::2, ::2] = tsubmat[::2, 1::2]
            t_new[::2, 1::2] = tsubmat[1::2, ::2]
            t_new[1::2, 1::2] = tsubmat[::2, ::2]
            tol += np.sum(np.abs(tsubmat - t_new) ** 2)
            tnew_sum += np.sum(np.abs(tsubmat) ** 2)
    tol = tol / (np.sum(np.abs(tmat) ** 2) + tnew_sum)
    return tol


def pmirroryz(tmat, js, ms, lams):
    tol = 0
    N_J = js.max()
    MirrorMat_list = (
        []
    )  # Creates a list of mirror reflection transformation matrices for each angular momentum J / but without helicity change or tau-multiplication
    tnew_sum = 0
    for lam1 in [b"positive", b"negative"]:
        for J1 in range(1, N_J + 1):
            mindex1 = np.multiply(lams == lam1, js == J1)
            mindex1_mirror = np.multiply(np.logical_not(lams == lam1), js == J1)
            for lam2 in [b"positive", b"negative"]:
                for J2 in range(1, N_J + 1):
                    mindex2 = np.multiply(lams == lam2, js == J2)
                    mindex2_mirror = np.multiply(np.logical_not(lams == lam2), js == J2)

                    tsubmat_flipped_helicity = tmat[
                        np.outer(mindex1_mirror, mindex2_mirror)
                    ]
                    tsubmat_flipped_helicity = tsubmat_flipped_helicity.reshape(
                        2 * J1 + 1, 2 * J2 + 1
                    )
                    t_mirr_submat = np.matmul(
                        MirrorMat_list[J1 - 1],
                        np.matmul(tsubmat_flipped_helicity, MirrorMat_list[J2 - 1]),
                    )

                    tsubmat = tmat[np.outer(mindex1, mindex2)]
                    tsubmat = tsubmat.reshape(2 * J1 + 1, 2 * J2 + 1)

                    tol += np.sum(np.abs(tsubmat - t_mirr_submat) ** 2)
                    tnew_sum += np.sum(np.abs(t_mirr_submat) ** 2)
    tol = tol / (np.sum(np.abs(tmat) ** 2) + tnew_sum)
    return tol


def mirrorxz(tmats, js, ms, taus):

    tol = 0
    N_J = js.max()
    ts = np.zeros(len(taus), dtype=int)
    ts[taus == b"electric"] = -1
    ts[taus == b"magnetic"] = 1  # ??????????? IS IT SO?
    m_out = ms + np.zeros_like(ms[:, None], dtype=int)
    m_in = ms[:, None] + np.zeros_like(ms)
    p_out = ts[None, :] + np.zeros_like(ts[:, None], dtype=int)
    p_in = ts[:, None] + np.zeros_like(ts, dtype=int)
    l_in = js[:, None] + np.zeros_like(js)
    l_out = js + np.zeros_like(js[:, None])
    n_range = 2 * (2 * np.arange(1, N_J) + 1)
    Nrng = np.concatenate(([0], n_range, [len(js) + 1]))
    tnew_sum = 0
    for i, ni in enumerate(Nrng[:-1]):
        for j, nj in enumerate(Nrng[:-1]):
            tsubmat = tmats[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            p_inc = p_in[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            p_outc = p_out[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            m_inc = m_in[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            m_outc = m_out[ni : Nrng[i + 1], nj : Nrng[j + 1]]
            t_new = np.zeros_like(tsubmat, dtype=complex)

            t_new[::2, ::2] = tsubmat[::2, ::2][::-1, ::-1]
            t_new[1::2, ::2] = tsubmat[1::2, ::2][::-1, ::-1]
            t_new[::2, 1::2] = tsubmat[::2, 1::2][::-1, ::-1]
            t_new[1::2, 1::2] = tsubmat[1::2, 1::2][::-1, ::-1]
            t_new = t_new * p_inc * p_outc * (-1.0) ** (m_inc + m_outc)
            tol += np.sum(np.abs(tsubmat - t_new) ** 2)
            tnew_sum += np.sum(np.abs(t_new) ** 2)
    tol = tol / (np.sum(np.abs(tmats) ** 2) + tnew_sum)
    return tol


def mirroryz(tmats, js, ms, taus):
    N_J = js.max()
    import quaternionic
    import spherical

    wigner = spherical.Wigner(N_J + 1)  # WHY +!
    R = quaternionic.array.from_euler_angles(0, np.pi / 2, 0)
    D = wigner.D(
        R
    )  # Actually the package gives the compex conjugate of D, but here for real small wigner matrix it does not matter
    MirrorMat_list = (
        []
    )  # Creates a list of mirror reflection transformation matrices for each angular momentum J / but without helicity change or tau-multiplication
    for J in range(1, N_J + 1):
        MirrorMat = np.zeros((2 * J + 1, 2 * J + 1), dtype=complex)
        for Mpp in range(-J, J + 1):
            i_Mpp = Mpp + J
            for M in range(-J, J + 1):
                i_M = M + J
                for Mp in range(-J, J + 1):
                    MirrorMat[i_Mpp, i_M] += (
                        D[wigner.Dindex(J, Mpp, Mp)]
                        * D[wigner.Dindex(J, Mp, M)]
                        * (-1) ** (J + M)
                    )
        MirrorMat_list.append(MirrorMat)
    tol = 0
    tnew_sum = 0
    for Tau1 in [b"electric", b"magnetic"]:
        tau1_num = 1 if Tau1 == b"electric" else -1
        for J1 in range(1, N_J + 1):
            mindex1 = np.multiply(taus == Tau1, js == J1)
            mindex1_mirror = np.multiply(taus == Tau1, js == J1)
            for Tau2 in [b"electric", b"magnetic"]:
                tau2_num = 1 if Tau2 == b"electric" else -1
                for J2 in range(1, N_J + 1):
                    mindex2 = np.multiply(taus == Tau2, js == J2)
                    mindex2_mirror = np.multiply(taus == Tau2, js == J2)
                    tsubmat_flipped_helicity = tmats[
                        np.outer(mindex1_mirror, mindex2_mirror)
                    ]

                    tsubmat_flipped_helicity = tsubmat_flipped_helicity.reshape(
                        2 * J1 + 1, 2 * J2 + 1
                    )
                    t_mirr_submat = (
                        tau1_num
                        * tau2_num
                        * np.matmul(
                            MirrorMat_list[J1 - 1],
                            np.matmul(tsubmat_flipped_helicity, MirrorMat_list[J2 - 1]),
                        )
                    )
                    tsubmat = tmats[np.outer(mindex1, mindex2)]
                    tsubmat = tsubmat.reshape(2 * J1 + 1, 2 * J2 + 1)

                    tol += np.sum(np.abs(tsubmat - t_mirr_submat) ** 2)
                    tnew_sum += np.sum(np.abs(t_mirr_submat) ** 2)
    tol = tol / (np.sum(np.abs(tmats) ** 2) + tnew_sum)
    return tol


def reciprocal_check(tmat, l, m, p):
    p = np.where(p == p[0], 1, 0)
    m_out = m + np.zeros_like(m[:, None])
    m_in = m[:, None] + np.zeros_like(m)
    p_out = p + np.zeros_like(p[:, None])
    p_in = p[:, None] + np.zeros_like(p)
    l_in = l[:, None] + np.zeros_like(l)
    l_out = l + np.zeros_like(l[:, None])
    tmat_reorder = np.zeros_like(tmat)
    for i in range(tmat.shape[0]):
        for j in range(
            tmat.shape[0]
        ):  # This is the one that's true for global mats, so the correct one
            bools = (
                (m_in == -m_out[i, j])
                & (m_out == -m_in[i, j])
                & (p_in == p_out[i, j])
                & (p_out == p_in[i, j])
                & (l_in == l_out[i, j])
                & (l_out == l_in[i, j])
            )
            tmat_reorder[np.argwhere(bools)[0][0], np.argwhere(bools)[0][1]] = (
                -1.0
            ) ** (-m_in[i, j] - m_out[i, j]) * tmat[i, j]
    tol = metric(tmat_reorder, tmat)
    return tol


def validate_hdf5_file(filepath):
    gr = []
    dt = []
    Groups = ["computation", "scatterer", "embedding", "modes"]
    Datasets = np.array(
        [
            "tmatrix",
            "vacuum_wavelength",
            "vacuum_wavenumber",
            "angular_vacuum_wavenumber",
            "frequency",
            "angular_frequency",
            "mesh",
        ]
    )
    Equiv = np.array(
        [
            "vacuum_wavelength",
            "vacuum_wavenumber",
            "angular_vacuum_wavenumber",
            "frequency",
            "angular_frequency",
        ]
    )
    print("Validate file", filepath)
    with h5py.File(filepath, "r") as f:
        # CHECK NAMES
        for key in f.keys():
            try:
                f[key]
            except:
                print(f" Entry {key} does not open")
                continue
            if isinstance(f[key], h5py.Group):
                gr.append(key)
            elif isinstance(f[key], h5py.Dataset):
                dt.append(key)
                if key in Equiv:
                    if "unit" not in f[key].attrs.keys():
                        raise ValueError(
                            "Unit of frequency (wavelength etc) has to be given"
                        )
                    if not isinstance(f[key].attrs["unit"], str):
                        raise TypeError(
                            "Unit of frequency (wavelength etc) has to be a string"
                        )
                    xunit = f[key].attrs["unit"]
                    if dt == "frequency":
                        xunit = FREQUENCIES[xunit]
                    if dt == "angular_frequency":
                        xunit = FREQUENCIES[xunit]
                    if dt == "vacuum_wavelength":
                        xunit = LENGTHS[xunit]
                    if dt == "vacuum_wavenumber":
                        xunit = INVLENGTHS[xunit]
                    if dt == "angular_vacuum_wavenumber":
                        xunit = INVLENGTHS[xunit]

        if bool(set(Equiv) & set(dt)):
            dt.extend(Equiv)

        gr_left = np.setdiff1d(Groups, gr, assume_unique=False)
        dt_left = np.setdiff1d(Datasets, dt, assume_unique=False)
        if len(gr_left):
            # raise NameError(f'The following mandatory groups are not present:', *gr_left)
            print(f"The following mandatory groups are not present:", *gr_left)
        if len(dt_left):
            # raise NameError(f'The following mandatory datasets are not present:', *dt_left)
            print(f"The following mandatory datasets are not present:", *dt_left)

        # CHECK ORDER IS THE SAME WITH T-MATRIX DIMENSION
        sizel = np.array(f["modes"]["l"]).max()
        if "index" not in f["modes"].keys():
            num = 1
        else:
            num = np.max(f["modes/index"][...])
        tmatsize = 2 * sizel * (sizel + 2) * num
        tmatsize1 = f["tmatrix"][...].shape[-1]
        if not (tmatsize == tmatsize1):
            raise ValueError("The modes should include only integer types.")

        # CHECK IF M AND L ARE IN CORRECT ORDER AND TYPE:
        m = f["modes"]["m"][...]
        l = f["modes"]["l"][...]
        p = f["modes"]["polarization"][...]

        if not (
            np.issubdtype(l.dtype, np.integer) & np.issubdtype(m.dtype, np.integer)
        ):
            raise TypeError("The orders l and m  should include only integer types.")

        if (b"electric" in p) & (b"magnetic" in p):
            # if not (p[:int(len(p)/2)] == b'electric').all() &  (p[int(len(p)/2):] == b'magnetic').all():
            if not (p[::2] == b"electric").all() & (p[1::2] == b"magnetic").all():
                raise ValueError(
                    f"Polarization ordering does not follow the accepted convention: 'electric' ('tm'), 'magnetic' ('te') alternating sequence"
                )

        elif (b"positive" in p) & ("negative" in p):
            if not (p[::2] == b"positive").all() & (p[1::2] == b"negative").all():
                raise ValueError(
                    "Polarization ordering does not follow the accepted convention: positive , negative alternating sequence"
                )
        else:
            raise ValueError(
                "The accepted entries are: 'electric' and 'magnetic' or 'positive' 'negative'."
            )

        mm = [li for lm in np.unique(l) for li in range(-lm, lm + 1)]
        mm = np.repeat(mm, 2)

        if not (m == mm).all():
            print("The ordering in m does not match the correct one:", mm)

        # CHECK MESH PRESENCE
        if not any(
            re.compile(r"mesh").search(key) for key in f["computation"].keys()
        ) | any(
            re.compile(r"mesh").search(key) for key in f["scatterer/geometry/"].keys()
        ):
            print("Mesh is not available")

        # CHECK KEYWORDS:
        if "keywords" in f.attrs:
            if f["tmatrix"][...].ndim == 3:
                tmat = f["tmatrix"][0]
            else:
                tmat = f["tmatrix"][...]
            if "czinfinity" in f.attrs["keywords"]:
                tol = czinf(tmat, l, m, p)
                print(
                    f"The object is claimed rotationally symmetric with respect to z-axis, this holds with a tolerance of {tol} for the T-matrix, which is half of the ratio of squared norm of the differences to squared norm of the sum of original and transformed matrix. If calculated in 2D, the tolerance should be indistinguishable from zero."
                )

            if "passive" in f.attrs["keywords"]:
                tol = passive(tmat)
                print(
                    f"The object is claimed passive, a corresponding  equation is Hermitian semi-definite  with a tolerance of {tol}"
                )

            if "lossless" in f.attrs["keywords"]:
                tnew = -(np.conj(tmat).T + 2 * np.conj(tmat).T @ tmat)
                deviation = metric(
                    tnew, tmat
                )  # np.max(np.abs(( -T -np.conj(T).T -2* np.conj(T).T @ T ) ))
                print(
                    f"The object is claimed lossless, the energy conservation condition holds with a maximum deviation of {deviation}"
                )

            if "reciprocal" in f.attrs["keywords"]:
                tol = reciprocal_check(tmat, l, m, p)
                print(
                    f"The object is claimed reciprocal,  this holds with a tolerance of {tol},  which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                )

            for z in f.attrs["keywords"]:
                if p[0] in [b"electric", b"magnetic"]:
                    if ("xy" in z) | ("XY" in z) | ("yx" in z) | ("YX" in z):
                        tol = mirrorxy(f["tmatrix"][...][0], l, m, p)
                        print(
                            f"The object is claimed mirror symmetric with respect to xy plane with a tolerance of {tol},  which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                        )

                    if ("yz" in z) | ("YZ" in z) | ("zy" in z) | ("ZY" in z):
                        tol = mirroryz(f["tmatrix"][...][0], l, m, p)
                        print(
                            f"The object is claimed mirror symmetric with respect to yz plane, this holds with a tolerance of {tol}, which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                        )
                    if (("x" in z) & ("z" in z)) | (("Z" in z) & ("X" in z)):
                        tol = mirrorxz(f["tmatrix"][...][0], l, m, p)
                        print(
                            f"The object is claimed mirror symmetric with respect to xz plane with a tolerance of {tol}, which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                        )
                elif p[0] in [b"positive", b"negative"]:
                    if ("xy" in z) | ("XY" in z) | ("yx" in z) | ("YX" in z):
                        tol = pmirrorxy(f["tmatrix"][...][0], l, m, p)
                        print(
                            f"The object is claimed mirror symmetric with respect to xy plane with a tolerance of {tol} which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                        )
                    if ("yz" in z) | ("YZ" in z) | ("zy" in z) | ("ZY" in z):
                        tol = pmirroryz(f["tmatrix"][...][0], l, m, p)
                        print(
                            f"The object is claimed mirror symmetric with respect to yz plane with a tolerance of {tol}, which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                        )
                    if (("x" in z) & ("z" in z)) | (("Z" in z) & ("X" in z)):
                        tol = pmirrorxz(f["tmatrix"][...][0], l, m, p)
                        print(
                            f"The object is claimed mirror symmetric with respect to xz plane with a tolerance of {tol}, which is half of the ratio of squared norm of the differences to the sum of squared norms of original and transformed matrix"
                        )
        return False


def validate_hdf5_directory(directory):
    print(f"Validate files in directory {directory}")
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".h5") or file.endswith(".hdf5"):
                filepath = os.path.join(root, file)
                validate_hdf5_file(filepath)


def main():
    parser = argparse.ArgumentParser(description="Validate HDF5 files.")
    parser.add_argument(
        "path", help="Path to an HDF5 file or directory containing HDF5 files."
    )
    args = parser.parse_args()
    path = args.path
    if os.path.isfile(path):
        validate_hdf5_file(path)
    elif os.path.isdir(path):
        validate_hdf5_directory(path)
    else:
        print("Invalid path. Please provide a valid HDF5 file or directory.")


LENGTHS = {
    "ym": 1e-24,
    "zm": 1e-21,
    "am": 1e-18,
    "fm": 1e-15,
    "pm": 1e-12,
    "nm": 1e-9,
    "um": 1e-6,
    "µm": 1e-6,
    "mm": 1e-3,
    "cm": 1e-2,
    "dm": 1e-1,
    "m": 1,
    "dam": 1e1,
    "hm": 1e2,
    "km": 1e3,
    "Mm": 1e6,
    "Gm": 1e9,
    "Tm": 1e12,
    "Pm": 1e15,
    "Em": 1e18,
    "Zm": 1e21,
    "Ym": 1e24,
}

INVLENGTHS = {
    r"ym^{-1}": 1e24,
    r"zm^{-1}": 1e21,
    r"am^{-1}": 1e18,
    r"fm^{-1}": 1e15,
    r"pm^{-1}": 1e12,
    r"nm^{-1}": 1e9,
    r"um^{-1}": 1e6,
    r"µm^{-1}": 1e6,
    r"mm^{-1}": 1e3,
    r"cm^{-1}": 1e2,
    r"dm^{-1}": 1e1,
    r"m^{-1}": 1,
    r"dam^{-1}": 1e-1,
    r"hm^{-1}": 1e-2,
    r"km^{-1}": 1e-3,
    r"Mm^{-1}": 1e-6,
    r"Gm^{-1}": 1e-9,
    r"Tm^{-1}": 1e-12,
    r"Pm^{-1}": 1e-15,
    r"Em^{-1}": 1e-18,
    r"Zm^{-1}": 1e-21,
    r"Ym^{-1}": 1e-24,
}

FREQUENCIES = {
    "yHz": 1e-24,
    "zHz": 1e-21,
    "aHz": 1e-18,
    "fHz": 1e-15,
    "pHz": 1e-12,
    "nHz": 1e-9,
    "uHz": 1e-6,
    "µHz": 1e-6,
    "mHz": 1e-3,
    "cHz": 1e-2,
    "dHz": 1e-1,
    "Hz": 1,
    "daHz": 1e1,
    "hHz": 1e2,
    "kHz": 1e3,
    "MHz": 1e6,
    "GHz": 1e9,
    "THz": 1e12,
    "PHz": 1e15,
    "EHz": 1e18,
    "ZHz": 1e21,
    "YHz": 1e24,
    r"ys^{-1}": 1e24,
    r"zs^{-1}": 1e21,
    r"as^{-1}": 1e18,
    r"fs^{-1}": 1e15,
    r"ps^{-1}": 1e12,
    r"ns^{-1}": 1e9,
    r"us^{-1}": 1e6,
    r"µs^{-1}": 1e6,
    r"ms^{-1}": 1e3,
    r"cs^{-1}": 1e2,
    r"ds^{-1}": 1e1,
    r"s^{-1}": 1,
    r"das^{-1}": 1e-1,
    r"hs^{-1}": 1e-2,
    r"ks^{-1}": 1e-3,
    r"Ms^{-1}": 1e-6,
    r"Gs^{-1}": 1e-9,
    r"Ts^{-1}": 1e-12,
    r"Ps^{-1}": 1e-15,
    r"Es^{-1}": 1e-18,
    r"Zs^{-1}": 1e-21,
    r"Ys^{-1}": 1e-24,
}

if __name__ == "__main__":
    main()
