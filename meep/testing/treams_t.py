import numpy as np
import treams


def get_treams_t_two_spheres(l_max=1, f=1) -> treams.PhysicsArray:
    objects = [
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.4, [treams.Material(1.1), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.3, [treams.Material(1.2), treams.Material()]
        ),
    ]
    positions = [[-0.4, -0.4, -0.4], [0.5, 0.5, 0]]
    t_treams = treams.TMatrix.cluster(objects, positions).interaction.solve()
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(l_max))
    return t_treams


def get_treams_t_three_spheres(l_max=1, f=1) -> treams.PhysicsArray:
    objects = [
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.4, [treams.Material(1.1), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.3, [treams.Material(1.2), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.45, [treams.Material(1.3), treams.Material()]
        ),
    ]
    positions = [[-0.4, -0.4, -0.4], [0.5, 0.5, -0.1], [-0.2, 0.3, 0.3]]
    t_treams = treams.TMatrix.cluster(objects, positions).interaction.solve()
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(l_max))
    return t_treams


def get_treams_t_four_spheres(l_max=1, f=1) -> treams.PhysicsArray:
    objects = [
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.4, [treams.Material(1.1), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.3, [treams.Material(1.2), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.45, [treams.Material(1.3), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.35, [treams.Material(1.4), treams.Material()]
        ),
    ]
    positions = [
        [-0.4, -0.4, -0.4],
        [0.5, 0.5, -0.1],
        [-0.2, 0.3, 0.3],
        [0.4, -0.2, -0.3],
    ]
    t_treams = treams.TMatrix.cluster(objects, positions).interaction.solve()
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(l_max))
    return t_treams


def get_treams_t_core_shell(l_max=1, f=1) -> treams.PhysicsArray:
    t_treams = treams.TMatrix.sphere(
        l_max,
        f * 2 * np.pi,
        np.array([0.25, 0.5]),
        [treams.Material(1.1), treams.Material(1.2), treams.Material()],
    )
    return t_treams


def get_treams_t_core_shell_translated(l_max=1, f=1) -> treams.PhysicsArray:
    t_treams = treams.TMatrix.sphere(
        l_max,
        f * 2 * np.pi,
        np.array([0.25, 0.5]),
        [treams.Material(1.1), treams.Material(1.2), treams.Material()],
    )
    t_treams = t_treams.translate([0.5, 0.0, 0.0])
    return t_treams


def get_treams_t_six_spheres(l_max=1, f=1) -> treams.PhysicsArray:
    objects = [
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.25, [treams.Material(1.15), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.25, [treams.Material(1.15), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.25, [treams.Material(1.15), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.25, [treams.Material(1.15), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.25, [treams.Material(1.15), treams.Material()]
        ),
        treams.TMatrix.sphere(
            l_max, f * 2 * np.pi, 0.25, [treams.Material(1.15), treams.Material()]
        ),
    ]
    positions = [
        [-0.5, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        [0.0, -0.5, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, -0.5],
        [0.0, 0.0, 0.5],
    ]
    t_treams = treams.TMatrix.cluster(objects, positions).interaction.solve()
    t_treams = t_treams.expand(treams.SphericalWaveBasis.default(l_max))
    return t_treams


def get_treams_t_sphere(l_max=1, f=1) -> treams.PhysicsArray:
    t_treams = treams.TMatrix.sphere(
        l_max, f * 2 * np.pi, 0.5, [treams.Material(1.15), treams.Material()]
    )
    return t_treams


def get_treams_t_sphere_embedded(l_max=1, f=1) -> treams.PhysicsArray:
    t_treams = treams.TMatrix.sphere(
        l_max, f * 2 * np.pi, 0.5, [treams.Material(1.2), treams.Material(1.1)]
    )
    return t_treams


def get_t_treams(l_max=1, f=1) -> treams.PhysicsArray:
    t_treams = get_treams_t_two_spheres(l_max=l_max, f=f)
    return t_treams.changepoltype()
