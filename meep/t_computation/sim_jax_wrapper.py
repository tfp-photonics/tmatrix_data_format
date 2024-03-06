from typing import Callable, List, Optional, Union

import jax.numpy as jnp
import meep as mp
import numpy as np
from jax import custom_vjp
from meep.adjoint import DesignRegion, ObjectiveQuantity, utils
from numpy.typing import NDArray

from geometry.geometry_factory import (
    get_geometry_from_eps_grid,
    get_matgrid_from_epsgrid,
)


class SimJaxWrapper:
    def __init__(
        self,
        simulation: mp.Simulation,
        design_volume: mp.Volume,
        monitors: List[ObjectiveQuantity],
        frequencies: Optional[Union[float, List[float]]],
        tol=1e-8,
        finite_difference_step: Optional[float] = utils.FD_DEFAULT,
    ):
        self._sim = simulation
        self._monitors = monitors
        self._design_volume = design_volume
        self._frequencies = frequencies
        self._finite_difference_step = finite_difference_step
        self._tol = tol
        self._forward_sources = self._sim.sources
        self._gradient = []
        self._simulate_func = self._get_jax_simulation_function()
        # The optimizer has three allowable states : "INIT", "FWD", and "ADJ".
        #    INIT - The optimizer is initialized and ready to run a forward simulation
        #    FWD  - The optimizer has already run a forward simulation
        #    ADJ  - The optimizer has already run an adjoint simulation (but not yet calculated the gradient)
        self._current_state = "INIT"

    def simulate(self, eps_grid) -> NDArray:
        return self._simulate_func(eps_grid)

    def _get_jax_simulation_function(self) -> Callable[[NDArray], NDArray]:
        @custom_vjp
        def simulate(eps_grid):
            self._install_eps_grid(eps_grid)
            if self._current_state == "INIT":
                print("Starting forward run...")
                self._forward_run()
            else:
                raise ValueError(
                    f"Incorrect solver state detected: {self._current_state}"
                )
            return jnp.asarray(self.results_list)

        def simulate_fwd(eps_grid):
            return simulate(eps_grid), None

        def simulate_bwd(res, grads):
            if self._current_state == "FWD":
                print("Starting adjoint run...")
                self._adjoint_run(grads)
                print("Calculating gradient...")
                self._calculate_gradient()
            else:
                raise ValueError(
                    f"Incorrect solver state detected: {self._current_state}"
                )
            return (jnp.asarray(self._gradient),)

        simulate.defvjp(simulate_fwd, simulate_bwd)
        return simulate

    def _install_eps_grid(self, eps_grid: NDArray) -> None:
        self._eps_grid_shape = eps_grid.shape
        self._sim.geometry = get_geometry_from_eps_grid(
            eps=eps_grid,
            center=self._design_volume.center,
            size=self._design_volume.size,
        )
        self._design_region = DesignRegion(
            design_parameters=get_matgrid_from_epsgrid(eps_grid),
            center=self._design_volume.center,
            size=self._design_volume.size,
        )

    def _prepare_forward_run(self) -> None:
        # prepare forward run
        self._sim.reset_meep()

        # add forward sources
        self._sim.change_sources(self._forward_sources)

        # register user specified monitors
        self.forward_monitors = [
            m.register_monitors(self._frequencies) for m in self._monitors
        ]

        # register design region
        self.forward_design_region_monitors = utils.install_design_region_monitors(
            self._sim, [self._design_region], self._frequencies
        )[0]

    def _forward_run(self) -> None:
        # set up monitors
        self._prepare_forward_run()

        # Forward run
        self._sim.run(until_after_sources=mp.stop_when_dft_decayed(self._tol))

        # record objective quantities from user specified monitors
        self.results_list = [m() for m in self._monitors]

        # update solver's current state
        self._current_state = "FWD"

    def _prepare_adjoint_run(self, dJ_total) -> None:
        # Compute adjoint sources
        self.adjoint_sources = []
        for mi, m in enumerate(self._monitors):
            dJ = dJ_total[mi]
            # get gradient of objective w.r.t. monitor
            if np.any(dJ):
                self.adjoint_sources += m.place_adjoint_source(np.array(dJ))

    def _adjoint_run(self, grads) -> None:
        # set up adjoint sources and monitors
        self._prepare_adjoint_run(grads)

        # flip the m number
        if utils._check_if_cylindrical(self._sim):
            self._sim.change_m(-self._sim.m)

        # flip the k point
        if self._sim.k_point:
            self._sim.change_k_point(-1 * self._sim.k_point)

        # Reset the fields
        self._sim.restart_fields()
        self._sim.clear_dft_monitors()

        # Update the sources
        self._sim.change_sources(self.adjoint_sources)

        # register design dft fields
        self.adjoint_design_region_monitors = utils.install_design_region_monitors(
            self._sim,
            [self._design_region],
            self._frequencies,
        )[0]

        self._sim._evaluate_dft_objects()

        # Adjoint run
        self._sim.run(until_after_sources=mp.stop_when_dft_decayed(self._tol))

        # reset the m number
        if utils._check_if_cylindrical(self._sim):
            self._sim.change_m(-self._sim.m)

        # reset the k point
        if self._sim.k_point:
            self._sim.change_k_point(-1 * self._sim.k_point)

        # update optimizer's state
        self._current_state = "ADJ"

    def _calculate_gradient(self) -> None:
        # calculate gradient
        self._gradient = self._design_region.get_gradient(
            self._sim,
            self.adjoint_design_region_monitors,
            self.forward_design_region_monitors,
            self._frequencies,
            self._finite_difference_step,
        )

        for i in range(3):
            # note that dft_fields::remove calls delete on its chunks, and the
            # destructor ~dft_chunk automatically removes it from the fields object
            self.forward_design_region_monitors[i].remove()

        # Cleanup list of lists
        self._gradient = np.reshape(
            np.sum(self._gradient, axis=-1), self._eps_grid_shape
        )

        # Return optimizer's state to initialization
        self._current_state = "INIT"
