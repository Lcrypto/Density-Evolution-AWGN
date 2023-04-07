import numpy as np
import de_utils as de_u
import matplotlib.pyplot as plt
import ga_discrete as ga_d
from multiprocessing import Pool, shared_memory
from multiprocessing.sharedctypes import Value

from config import cfg
cfg_de = cfg.get('density_evolution')
cfg_disc = cfg_de.get('ga_discrete')
run_id = cfg.get("run_id")
n_processes = cfg.get("n_processes")


def evaluate(j):
    global pop_name, fit_name, pop_shape, fit_shape, i_idx

    shm_pop = shared_memory.SharedMemory(name=pop_name)
    shm_fit = shared_memory.SharedMemory(name=fit_name)
    pop = np.ndarray(pop_shape, dtype=np.int64, buffer=shm_pop.buf)
    fit = np.ndarray(fit_shape, dtype=np.float64, buffer=shm_fit.buf)
    i = i_idx.value

    if fit[i, j] == -np.inf:
        cn_degrees = np.sum(pop[j, :, :], 0)
        vn_degrees = np.sum(pop[j, :, :], 1)

        condition1 = np.all(cn_degrees >= 2) and np.all(vn_degrees >= 2)
        condition2 = np.all(cn_degrees <= 25)

        if condition1 and condition2:
            rho_node = np.bincount(vn_degrees)[1:]
            rho_edge = np.arange(1, len(rho_node)+1) * rho_node / \
                np.sum(np.arange(1, len(rho_node)+1) * rho_node)
            lam_node = np.bincount(cn_degrees)[1:]
            lam_edge = np.arange(1, len(lam_node)+1) * lam_node / \
                np.sum(np.arange(1, len(lam_node)+1) * lam_node)

            min = cfg_de.get("min_rber")
            max = cfg_de.get("max_rber")
            fitness = de_u.bisection_search(min, max, rho_edge, lam_edge)
            fit[i, j] = fitness*100

        else:
            fit[i, j] = -100
    shm_pop.close()
    shm_fit.close()


def pool_initializer():
    global pop_name, fit_name, pop_shape, i_idx


def ga_discrete_parallel():
    global pop_name, fit_name, pop_shape, fit_shape, i_idx
    load_population = cfg_de.get("load_population")
    de_u.set_cfg(cfg)

    n_pop = cfg_de.get("Np")
    n_generations = cfg_de.get("generations")
    n_vn = cfg_disc.get("n_vn")
    n_cn = cfg_disc.get("n_cn")
    p_vertical = cfg_disc.get("p_vertical")
    p_horizontal = cfg_disc.get("p_horizontal")
    p_mutation = cfg_disc.get("p_mutation")
    print_terminal = cfg_de.get("print_terminal")

    de_u.save_params()

    if load_population:
        fname = "data/" + run_id + ".npz"
        data = np.load(fname)
        population = data["population"]
        fitness = data["fitness"]
        i_start = int(data["generation"][0]) + 1
        dim_0 = np.size(fitness, axis=0)
        if n_generations != np.size(fitness, axis=0):
            fitness_new = np.zeros((n_generations, n_pop))
            fitness_new[:dim_0] = fitness
            fitness = fitness_new
    else:
        population = ga_d.init_population(n_pop, n_cn, n_vn)
        fitness = np.full((n_generations, n_pop), -np.inf)
        i_start = 0

    best_idx = 0
    best_rber = 0
    code_rate = (n_vn-n_cn)/n_vn

    header = f"""
    Running ga_disrete.py
    ===================================================================
    Optimizing {n_cn}x{n_vn} protograph.
    Code rate: {code_rate:.2f}.
    Number of individuals: {n_pop}.
    Number of generations: {n_generations}.
    Mutation probability: {p_mutation:.2f}.
    Vertical probability: {p_vertical:.2f}.
    Horizontal probability: {p_horizontal:.2f}.
    -------------------------------------------------------------------
    """
    if print_terminal:
        print(header)
    else:
        de_u.log(header, 'w')

    # Initialize shared_memory in multiprocessing
    shm_pop = shared_memory.SharedMemory(create=True, size=population.nbytes)
    shm_fit = shared_memory.SharedMemory(create=True, size=fitness.nbytes)

    pop = np.ndarray(population.shape, dtype=np.int64, buffer=shm_pop.buf)
    fit = np.ndarray(fitness.shape, dtype=np.float64, buffer=shm_fit.buf)

    pop[:] = population[:]
    fit[:] = fitness[:]

    pop_name = shm_pop.name
    fit_name = shm_fit.name

    pop_shape = pop.shape
    fit_shape = fit.shape

    i_idx = Value('i', i_start)
    i = i_start
    # Use Pool for parallel computing
    pool = Pool(processes=n_processes, initializer=pool_initializer)
    pop_idx = range(n_pop)

    try:
        # Initial evaluation
        if i_start == 0:
            pool.map(evaluate, pop_idx)
            best_idx = np.argmax(fit[i_start, :])
            best_rber = np.max(fit[i_start, :])

            status = f"{i} generations completed. RBER: best: {np.max(fit[i_start,:]):.2f}, min: {np.min(fit[i_start,:]):.2f}, mean: {np.mean(fit[i_start,:]):.2f}, variance: {np.var(fit[i_start,:]):.2f}.                                "
            if print_terminal:
                print(status)
            else:
                de_u.log(status, 'a')

            i_start += 1
            i_idx.value += 1

        for i in range(i_start, n_generations):

            # Select with elitism
            pop_new = np.zeros(pop.shape, int)
            pop_new[0] = pop[np.argmax(fit[i-1, :])]
            fit[i, 0] = fit[i-1, np.argmax(fit[i-1, :])]
            for j in range(1, n_pop):
                pop_new[j], fit[i, j] = ga_d.tournament(
                    pop, fit[i-1, :])
            pop[:] = pop_new[:]
            #fit = fit_new

            # Crossover
            for j in range(1, n_pop-1, 2):
                if np.random.rand() < p_vertical:
                    pop[j], pop[j +
                                1] = ga_d.vertical_crossover(pop[j], pop[j+1])
                    fit[i-1, j], fit[i-1, j+1] = -np.inf, -np.inf
                if np.random.rand() < p_horizontal:
                    pop[j], pop[j +
                                1] = ga_d.horizontal_crossover(pop[j], pop[j+1])
                    fit[i, j], fit[i, j+1] = -np.inf, -np.inf

            # Mutation
            for j in range(1, n_pop):
                if np.random.rand():
                    pop[j] = ga_d.mutation(pop[j])
                    fit[i, j] = -np.inf

            # Evaluate new individuals
            pool.map(evaluate, pop_idx)

            i_idx.value += 1

            # Keep track of best individual
            best_idx = np.argmax(fit[i, :])
            best_rber = np.max(fit[i, :])

            status = f"{i} generations completed. RBER: best: {np.max(fit[i,:]):.2f}, min: {np.min(fit[i,:]):.2f}, mean: {np.mean(fit[i,:]):.2f}, variance: {np.var(fit[i,:]):.2f}.                                "
            if print_terminal:
                print(status)
            else:
                de_u.log(status, 'a')

            if i % int(cfg_de.get("save_interval")) == 0:
                de_u.save_population(pop, fit, i, best_idx,
                                     best_rber, "discrete")

        status = f"""
    -------------------------------------------------------------------
    Finished!
    ===================================================================
            """
        if print_terminal:
            print(status)
        else:
            de_u.log(status, 'a')

    finally:
        de_u.save_population(pop, fit, i, best_idx, best_rber, "discrete")
        if i < n_generations-1:
            status = f"""
    -------------------------------------------------------------------
    Optimization interrupted.
    ===================================================================
                """
            if print_terminal:
                print(status)
            else:
                de_u.log(status, 'a')

        shm_pop.close()
        shm_fit.close()
        shm_pop.unlink()
        shm_fit.unlink()
