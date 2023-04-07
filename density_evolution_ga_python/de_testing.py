#%%
import de_utils
import numpy as np
import matplotlib.pyplot as plt 

#%%

condition1 = False
condition2 = True
n_cn = 9
n_vn = 73

while not condition1 and condition2:

    proto = np.zeros((n_cn, n_vn), int)
    for j in range(n_vn):
        k = np.random.choice(n_cn, 2, replace=False)
        proto[k, j] = 1    

    cn_degrees = np.sum(proto, 0)
    vn_degrees = np.sum(proto, 1)

    condition1 = np.all(cn_degrees >= 2) and np.all(vn_degrees >= 2)
    condition2 = np.all(cn_degrees <= 25)

#%% 
if condition1 and condition2:
    rho_node = np.bincount(vn_degrees)[1:]
    rho_edge = np.arange(1, len(rho_node)+1) * rho_node / \
        np.sum(np.arange(1, len(rho_node)+1) * rho_node)
    lam_node = np.bincount(cn_degrees)[1:]
    lam_edge = np.arange(1, len(lam_node)+1) * lam_node / \
        np.sum(np.arange(1, len(lam_node)+1) * lam_node)
        
#%%
def eval(rber, rho_edge, lam_edge):
    plot = False

    f_grid, g_grid, pdf = de_utils.init_pdf(rber, 4096, 60)
    cdf = de_utils.to_cdf(pdf)

    if plot:
        plt.plot(f_grid, cdf)
        plt.show()

    result = de_utils.symmetric_density_evolution(
        cdf, f_grid, g_grid, rho_edge, lam_edge, tol=rber*0.05,  plot=False, prnt=True)

    return result


def bisection_search(min, max, rho_edge, lam_edge, tol=0.0001):
    while max - min > tol:
        x = (min + max)/2
        result = eval(x, rho_edge, lam_edge)

        if result == 0:
            min = x
        elif result > 0:
            max = x
        else:
            raise Exception

    return (min + max)/2

bisection_search(0.005, 0.05, rho_edge, lam_edge)

#%%
rber = 0.015
f_grid, g_grid, pdf = de_utils.init_pdf(rber, n_grid=4096, llr_max=30)
cdf = de_utils.to_cdf(pdf)

de_utils.symmetric_density_evolution(cdf, f_grid, g_grid, rho_edge, lam_edge, n_iter=50, tol=1e-4, plot=True, prnt=True)
# %%
rber = 0.035
f_grid, g_grid, pdf = de_utils.init_pdf(rber, n_grid=1024*4, llr_max=20)
cdf = de_utils.to_cdf(pdf)

de_utils.symmetric_density_evolution(cdf, f_grid, g_grid, rho_edge, lam_edge, n_iter=50, tol=1e-6, plot=True, prnt=True)

# %%
