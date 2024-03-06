import jax.numpy as jnp
import numpy as np
import scipy.optimize as sc
import treams
import treams.special
from numpy.typing import ArrayLike, NDArray

from config.config import Config
from geometry.rotation_data import RotationData
from sim.sim_data import SimData
from sources.source_data import SourceData
from utility.vector_geometry import cartesian_to_spherical


def get_e_amp_from_e_in(
    c: Config, source_data: SourceData, sim_res: SimData
) -> tuple[NDArray, NDArray]:

    k_abs = c.wave_numbers
    pos = np.reshape(sim_res.pos, (-1, 3))
    nfreq = sim_res.e.shape[0]
    e_abs_amp = np.zeros(nfreq)
    e_phase = np.zeros(nfreq)

    for ii in range(nfreq):
        k_pw = k_abs[ii] * source_data.k_src / np.linalg.norm(source_data.k_src)
        field = sim_res.e[ii]
        field = np.reshape(field, (-1, 3))
        field = np.dot(field, source_data.pol_src)
        field_phases = np.angle(field)
        field = np.real(field)

        a_est = np.max(np.abs(field))

        phi_est = np.mod(
            np.dot(pos, k_pw) - field_phases,
            2 * np.pi,
        )
        at_border = False
        border_threshold = 0.05 * phi_est.shape[0]
        if (
            np.sum(phi_est < 1) > border_threshold
            and np.sum(phi_est > 5) > border_threshold
        ):
            phi_est = np.mod(phi_est + np.pi, 2 * np.pi)
            at_border = True
        phi_est = np.mod(np.average(phi_est), 2 * np.pi)
        if at_border:
            phi_est -= np.pi
            phi_est = np.mod(phi_est, 2 * np.pi)

        def pw_func(x, amp, phi):
            return amp * np.exp(1j * (np.dot(x, k_pw) - phi))

        def pw_func_real(x, amp, phi):
            return np.real(pw_func(x, amp, phi))

        popt, pcov = sc.curve_fit(
            pw_func_real,
            xdata=pos,
            ydata=field,
            p0=[a_est, phi_est],
            bounds=(
                [a_est, 0.8 * phi_est],
                [1.1 * a_est, 1.2 * phi_est],
            ),
        )
        e_abs_amp[ii], e_phase[ii] = popt

    e_amp_phased = e_abs_amp * np.exp(-1j * e_phase)
    # shape = (frequency, xyz)
    e_amp = e_amp_phased[:, None] * source_data.pol_src[None, :]
    # plot_incident_field(sim_res=sim_res)
    # plot_fitted_incident_field(
    #     sim_res=sim_res,
    #     source_data=source_data,
    #     e_abs_amp=e_abs_amp,
    #     e_phase=e_phase,
    # )
    return e_amp


def get_inc_coefs(
    c: Config,
    source_data: SourceData,
    rotation_data: RotationData,
    e_amp: NDArray,
) -> NDArray:
    # incident field complex vector amplitude
    inc_coefs = np.zeros((c.n_freq, c.min_sim_amount), dtype=np.complex128)
    e_amp_rot = (rotation_data.rotation_matrix @ e_amp.T).T
    wave_numbers = c.wave_numbers
    for ii in range(c.n_freq):
        k_abs = wave_numbers[ii]
        k_pw_rot = rotation_data.rotation_matrix @ (
            k_abs * source_data.k_src / np.linalg.norm(source_data.k_src)
        )
        pw = treams.plane_wave(
            k_pw_rot,
            [*e_amp_rot[ii]],
            k0=k_abs,
            material=treams.Material(),
            poltype="parity",
        )
        inc_coefs[ii] = np.array(pw.expand(treams.SphericalWaveBasis.default(c.l_max)))
    print('in inc coe', inc_coefs.shape)
    return inc_coefs


def get_sca_coef(
    c: Config,
    l: int,
    m: int,
    sim_res: SimData,
    rotation_data: RotationData,
) -> NDArray:
    def rotate_ndarray_last_axis(arr: ArrayLike):
        return jnp.tensordot(arr, rotation_data.rotation_matrix, axes=([-1], [-1]))

    sim_res_rot = SimData(
        e=rotate_ndarray_last_axis(sim_res.e),
        b=rotate_ndarray_last_axis(sim_res.b),
        pos=rotate_ndarray_last_axis(sim_res.pos),
        normals=rotate_ndarray_last_axis(sim_res.normals),
        weights=sim_res.weights,
    )
    rtp = cartesian_to_spherical(sim_res_rot.pos[None, :, :, :, :])
    r, theta, phi = (
        rtp[..., 0],
        rtp[..., 1],
        rtp[..., 2],
    )
    k = c.wave_numbers
    kr = r * k[:, None, None, None]
    omega = c.omegas
    # rotE = -dB/dt = i*omega*B
    rot_E_sc_rot = 1j * omega[:, None, None, None, None] * sim_res_rot.b
    N_1_star = jnp.conj(
        treams.special.vsph2car(treams.special.vsw_rN(l, m, kr, theta, phi), rtp)
    )
    M_1_star = jnp.conj(
        treams.special.vsph2car(treams.special.vsw_rM(l, m, kr, theta, phi), rtp)
    )

    def get_sc_coef_via_integral(vswf_a, vswf_b):
        integrand = jnp.cross(rot_E_sc_rot, vswf_a) - k[
            :, None, None, None, None
        ] * jnp.cross(vswf_b, sim_res_rot.e)

        integral = jnp.sum(
            integrand
            * sim_res_rot.normals[None, :, None, None, :]
            * sim_res_rot.weights[None, :, :, :, None],
            (1, 2, 3, 4),
        )
        return 1j * k * integral

    a = get_sc_coef_via_integral(N_1_star, M_1_star)
    b = get_sc_coef_via_integral(M_1_star, N_1_star)
    return jnp.stack([a, b])


def get_sca_coefs(c: Config, rotation_data: RotationData, sim_res: SimData) -> NDArray:
    sca_coefs = jnp.array([])
    for l in range(1, c.l_max + 1):
        for m in range(-l, l + 1):
            ab = get_sca_coef(c, l, m, sim_res, rotation_data)
            sca_coefs = jnp.concatenate((sca_coefs, ab)) if sca_coefs.size else ab
    return sca_coefs.T
