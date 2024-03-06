import jax


from config.config_factory import  get_test_config
from config.goal_t import get_goal_T
from t_computation.scat import calculate_T 
from t_computation.optimizer import optimize
from testing.t_plots import plot_me_vs_treams
from testing.treams_t import get_t_treams
from store_format import persist_T
from decomposition.t_data import TData
import meep as mp

jax.config.update("jax_enable_x64", True)


def opt() -> None:
    c = get_config()
    eps = optimize(c, goal_T=get_goal_T(c))


def T() -> None:
    c = get_config()
    plot_me_vs_treams(t, get_t_treams(l_max=c.l_max, f=c.f_cen), c.path)


def save() -> None:

    c = get_test_config()
    t = calculate_T(c)
    # t = TData.from_hdf(c.path, 0).t 
    persist_T(c, t, "one_sphere", "A test save file")



if __name__ == "__main__":
    save()
