import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp
import meep as mp
import meep.adjoint as mpa
import numpy as np
from jax import value_and_grad

from config.config import Config
from geometry.geometry_factory import (
    get_epsgrid_from_geometric_objects,
    get_geometry_from_eps_grid,
    get_matgrid_from_epsgrid,
)
from geometry.rotate import rotate_eps_grid
from geometry.rotation_data import RotationData
from optimization.sim_jax_wrapper import SimJaxWrapper
from sources.source_factory import get_plane_wave_source


def test_reshape():
    import jax.numpy as jnp
    from jax import grad

    from sim.sim import _reshape_and_complete_monitor_values

    class MockSim:
        def get_array_metadata(self, center=None, size=None, return_pw=True):
            return np.zeros((5, 5, 3)), np.zeros((5, 5))

    class MockMonitor:
        sim = MockSim()
        volume = mp.Volume()

    mockMonitors = []
    for _ in range(36):
        mockMonitors.append(MockMonitor())
    monitor_vals = jnp.ones((36, 3, 4, 4))

    def loss(monitor_vals, mockMonitors):
        x = _reshape_and_complete_monitor_values(monitor_vals, mockMonitors)
        return jnp.sum(x[0])

    g = grad(loss)(monitor_vals, mockMonitors)
    print(g)


def modified_opt_problem():
    def get_wrapped_meep():
        geometry = [
            mp.Sphere(radius=0.25, material=mp.Medium(epsilon=1.2)),
        ]

        c = Config(
            path="/scratch/local/pscherer/tmat_meep_out/diff_test",
            object_size=np.array([0.5, 0.5, 0.5]),
            resolution=24,
            sim_amount_mult=1,
            load_simulations=False,
            l_max=1,
            start_geometry=geometry,
        )

        cell = mp.Vector3(*c.with_pml_size)
        pml = [mp.PML(c.dpml)]
        sim = mp.Simulation(
            cell_size=cell,
            resolution=c.resolution,
            boundary_layers=pml,
            k_point=mp.Vector3(),
            sources=get_plane_wave_source(c=c).meep_sources,
        )
        eps_grid = jnp.asarray(c.start_geometry)

        l = c.no_pml_size[0]
        monitor1 = mpa.FourierFields(
            sim=sim,
            component=mp.Ex,
            volume=mp.Volume(center=mp.Vector3(0, 0, l), size=mp.Vector3(l, l, 0)),
        )
        monitor2 = mpa.FourierFields(
            sim=sim,
            component=mp.Ex,
            volume=mp.Volume(center=mp.Vector3(0, 0, -l), size=mp.Vector3(l, l, 0)),
        )

        wrapped_meep = SimJaxWrapper(
            simulation=sim,
            design_volume=mp.Volume(size=mp.Vector3(*c.no_pml_size)),
            monitors=[monitor1, monitor2],
            frequencies=c.frequencies,
        )
        return wrapped_meep, eps_grid.shape

    wrapped_meep, eps_grid_shape = get_wrapped_meep()
    np.random.seed(123)
    eps_grid = jnp.asarray(np.random.random(eps_grid_shape)) + 1.0

    def fom(eps_grid):
        fields = wrapped_meep.simulate(eps_grid)
        f = jnp.abs(fields.sum())
        return f

    f, grads = jax.value_and_grad(fom)(eps_grid)

    wrapped_meep2, eps_grid_shape = get_wrapped_meep()
    direction = np.zeros_like(eps_grid)
    print(eps_grid.shape)
    indices = (12, 12, 12)
    direction[indices] = 1
    step = 0.001
    eps_grid2 = eps_grid + step * jnp.asarray(direction)
    fields2 = wrapped_meep2.simulate(eps_grid2)
    f2 = jnp.abs(fields2.sum())
    df = f2 - f
    deps = eps_grid2[indices] - eps_grid[indices]
    grads_num = df * (1.0 / deps)
    print(grads[indices], grads_num)


def meep_jax_wrapper():
    geometry = [
        mp.Sphere(radius=0.5, material=mp.Medium(epsilon=1.2)),
    ]

    c = Config(
        path="/scratch/local/pscherer/tmat_meep_out/diff_test",
        object_size=np.array([1.0, 1.0, 1.0]),
        resolution=10,
        sim_amount_mult=1,
        load_simulations=False,
        l_max=1,
        start_geometry=geometry,
    )

    cell = mp.Vector3(*c.with_pml_size)
    pml = [mp.PML(c.dpml)]
    sim = mp.Simulation(
        cell_size=cell,
        resolution=c.resolution,
        boundary_layers=pml,
        k_point=mp.Vector3(),
    )
    eps_grid = jnp.array(
        get_epsgrid_from_geometric_objects(2 * c.resolution, c.object_size, geometry)
    )
    matgrid = get_matgrid_from_epsgrid(eps_grid)
    geo = get_geometry_from_eps_grid(eps_grid, c.object_size)
    design_region_matgrid = mp.MaterialGrid(
        mp.Vector3(*eps_grid.shape),
        medium1=mp.Medium(epsilon=c.opt_eps_min),
        medium2=mp.Medium(epsilon=c.opt_eps_max),
    )
    design_region = mpa.DesignRegion(
        design_parameters=design_region_matgrid,
        size=mp.Vector3(*c.no_pml_size),
        center=mp.Vector3(),
    )

    l = c.no_pml_size[0]
    monitor1 = mpa.FourierFields(
        sim=sim,
        component=mp.Ex,
        volume=mp.Volume(center=mp.Vector3(0, 0, l), size=mp.Vector3(l, l, 0)),
    )
    monitor2 = mpa.FourierFields(
        sim=sim,
        component=mp.Ex,
        volume=mp.Volume(center=mp.Vector3(0, 0, -l), size=mp.Vector3(l, l, 0)),
    )
    monitor_test = mpa.EigenmodeCoefficient(
        sim=sim,
        volume=mp.Volume(center=mp.Vector3(0, 0, l), size=mp.Vector3(l, l, 0)),
        mode=1,
    )
    test_eigenmode_moitors = False
    wrapped_meep = mpa.MeepJaxWrapper(
        simulation=sim,
        sources=get_plane_wave_source(c=c).meep_sources,
        until_after_sources=mp.stop_when_dft_decayed(tol=1e-8),
        design_regions=[design_region],
        monitors=[monitor_test] if test_eigenmode_moitors else [monitor1, monitor2],
        frequencies=c.frequencies,
    )

    print(eps_grid)
    fields = wrapped_meep([eps_grid])
    print(fields)


def meep_tutorial():
    import meep as mp
    import meep.adjoint as mpa
    import nlopt
    import numpy as np
    from autograd import grad
    from autograd import numpy as npa
    from autograd import tensor_jacobian_product
    from matplotlib import pyplot as plt
    from matplotlib.patches import Circle
    from scipy import signal, special

    mp.quiet(quietval=True)
    Si = mp.Medium(index=3.4)
    SiO2 = mp.Medium(index=1.44)

    waveguide_width = 0.5
    design_region_width = 2.5
    design_region_height = 2.5

    waveguide_length = 0.5

    pml_size = 1.0

    resolution = 40

    frequencies = 1 / np.linspace(1.5, 1.6, 3)
    minimum_length = 0.09  # minimum length scale (microns)
    eta_i = 0.5  # blueprint (or intermediate) design field thresholding point (between 0 and 1)
    eta_e = 0.55  # erosion design field thresholding point (between 0 and 1)
    eta_d = 1 - eta_e  # dilation design field thresholding point (between 0 and 1)
    filter_radius = mpa.get_conic_radius_from_eta_e(minimum_length, eta_e)
    print(filter_radius)
    design_region_resolution = int(1 * resolution)
    Sx = 2 * pml_size + 2 * waveguide_length + design_region_width + 2
    Sy = 2 * pml_size + design_region_height + 2
    cell_size = mp.Vector3(Sx, Sy)

    pml_layers = [mp.PML(pml_size)]

    fcen = 1 / 1.55
    width = 0.2
    fwidth = width * fcen
    source_center = [-Sx / 2 + pml_size + waveguide_length / 3, 0, 0]
    source_size = mp.Vector3(0, Sy, 0)
    kpoint = mp.Vector3(1, 0, 0)
    src = mp.GaussianSource(frequency=fcen, fwidth=fwidth)
    source = [
        mp.EigenModeSource(
            src,
            eig_band=1,
            direction=mp.NO_DIRECTION,
            eig_kpoint=kpoint,
            size=source_size,
            center=source_center,
        )
    ]

    Nx = int(design_region_resolution * design_region_width) + 1
    Ny = int(design_region_resolution * design_region_height) + 1

    design_variables = mp.MaterialGrid(mp.Vector3(Nx, Ny), SiO2, Si)
    design_region = mpa.DesignRegion(
        design_variables,
        volume=mp.Volume(
            center=mp.Vector3(),
            size=mp.Vector3(design_region_width, design_region_height, 0),
        ),
    )

    x_g = np.linspace(-design_region_width / 2, design_region_width / 2, Nx)
    y_g = np.linspace(-design_region_height / 2, design_region_height / 2, Ny)
    X_g, Y_g = np.meshgrid(x_g, y_g, sparse=True, indexing="ij")

    left_wg_mask = (X_g == -design_region_width / 2) & (
        np.abs(Y_g) <= waveguide_width / 2
    )
    top_wg_mask = (Y_g == design_region_width / 2) & (
        np.abs(X_g) <= waveguide_width / 2
    )
    Si_mask = left_wg_mask | top_wg_mask

    border_mask = (
        (X_g == -design_region_width / 2)
        | (X_g == design_region_width / 2)
        | (Y_g == -design_region_height / 2)
        | (Y_g == design_region_height / 2)
    )
    SiO2_mask = border_mask.copy()
    SiO2_mask[Si_mask] = False

    def mapping(x, eta, beta):
        x = npa.where(Si_mask.flatten(), 1, npa.where(SiO2_mask.flatten(), 0, x))
        # filter
        filtered_field = mpa.conic_filter(
            x,
            filter_radius,
            design_region_width,
            design_region_height,
            design_region_resolution,
        )

        # projection
        projected_field = mpa.tanh_projection(filtered_field, beta, eta)

        # interpolate to actual materials
        return projected_field.flatten()

    geometry = [
        mp.Block(
            center=mp.Vector3(x=-Sx / 4), material=Si, size=mp.Vector3(Sx / 2, 0.5, 0)
        ),  # horizontal waveguide
        mp.Block(
            center=mp.Vector3(y=Sy / 4), material=Si, size=mp.Vector3(0.5, Sy / 2, 0)
        ),  # vertical waveguide
        mp.Block(
            center=design_region.center,
            size=design_region.size,
            material=design_variables,
        ),  # design region
    ]

    sim = mp.Simulation(
        cell_size=cell_size,
        boundary_layers=pml_layers,
        geometry=geometry,
        sources=source,
        default_material=SiO2,
        resolution=resolution,
    )

    Ez_top = mpa.FourierFields(
        sim,
        mp.Volume(
            center=mp.Vector3(0, Sy / 2 - pml_size - 0.1, 0),
            size=mp.Vector3(x=waveguide_width),
        ),
        mp.Ez,
    )

    ob_list = [Ez_top]

    def J(top):
        power = npa.abs(top[1, 7]) ** 2
        return power

    opt = mpa.OptimizationProblem(
        simulation=sim,
        objective_functions=J,
        objective_arguments=ob_list,
        design_regions=[design_region],
        frequencies=frequencies,
        decay_by=1e-6,
    )
    rho_vector = np.random.rand(Nx * Ny)
    opt.update_design([mapping(rho_vector, eta_i, 1e3)])
    opt.plot2D(True)
    plt.show()
    evaluation_history = []
    cur_iter = [0]

    def f(v, gradient, cur_beta):
        print("Current iteration: {}".format(cur_iter[0] + 1))

        f0, dJ_du = opt([mapping(v, eta_i, cur_beta)])  # compute objective and gradient

        if gradient.size > 0:
            gradient[:] = tensor_jacobian_product(mapping, 0)(
                v, eta_i, cur_beta, np.sum(dJ_du, axis=1)
            )  # backprop

        evaluation_history.append(np.real(f0))

        plt.figure()
        ax = plt.gca()
        opt.plot2D(
            False,
            ax=ax,
            plot_sources_flag=False,
            plot_monitors_flag=False,
            plot_boundaries_flag=False,
        )
        circ = Circle((2, 2), minimum_length / 2)
        ax.add_patch(circ)
        ax.axis("off")
        # plt.savefig('media/bend_{:03d}.png'.format(cur_iter[0]),dpi=300)
        plt.show()

        cur_iter[0] = cur_iter[0] + 1

        return np.real(f0)

    algorithm = nlopt.LD_MMA
    n = Nx * Ny  # number of parameters

    # Initial guess
    x = np.ones((n,)) * 0.5
    x[Si_mask.flatten()] = 1  # set the edges of waveguides to silicon
    x[SiO2_mask.flatten()] = 0  # set the other edges to SiO2

    # lower and upper bounds
    lb = np.zeros((Nx * Ny,))
    lb[Si_mask.flatten()] = 1
    ub = np.ones((Nx * Ny,))
    ub[SiO2_mask.flatten()] = 0

    cur_beta = 4
    beta_scale = 2
    num_betas = 6
    update_factor = 15
    for iters in range(num_betas):
        solver = nlopt.opt(algorithm, n)
        solver.set_lower_bounds(lb)
        solver.set_upper_bounds(ub)
        solver.set_max_objective(lambda a, g: f(a, g, cur_beta))
        solver.set_maxeval(update_factor)
        solver.set_xtol_rel(1e-4)
        x[:] = solver.optimize(x)
        cur_beta = cur_beta * beta_scale


def test_new_rot():
    n_pixels = 10
    eps = np.ones((n_pixels, n_pixels, n_pixels))
    n = 5
    eps[:n, :n, :n] = 3
    rot = RotationData(theta=0.0, phi=np.pi / 2, alpha=0.0)

    def fom_func(eps):
        return jnp.sum(rotate_eps_grid(eps, rot)[:n, n:, :n])

    fom, grad = value_and_grad(fom_func)(eps)
    abs_grad = np.asarray(jnp.abs(grad))
    new_eps = rotate_eps_grid(eps=eps, rotation_data=rot)
    import matplotlib.pyplot as plt

    def plot_voxel(eps, threshold=1.1):
        ax = plt.figure().add_subplot(projection="3d")
        ax.voxels(filled=eps > threshold)

    plot_voxel(eps)
    plt.savefig("eps.png")
    plot_voxel(new_eps)
    plt.savefig("new_eps.png")
    plot_voxel(abs_grad, threshold=0.1)
    plt.savefig("abs_grad.png")
