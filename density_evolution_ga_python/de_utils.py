"""
This file contains methods for a numeric discretized implementation of density evolution for a general distribution.

Sources:
https://arxiv.org/pdf/cs/0509014.pdf [1]
https://ieeexplore.ieee.org/abstract/document/910578?casa_token=27RtOEW11OcAAAAA:OX_XwVvok1YSH5PgL930MjWe1a1MySsUdVvu_CCVeIWA3X2m9AwPTN9MR6BS7gDfqQuq_U8o0A [2]

"""

from threading import Thread
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import scipy.signal as sp
import pandas as pd
import os
import sys
sys.path.insert(
    0, '/home/fredrikblomgren/CAS/MasterThesis/ldpc-investigation-master-thesis/python-files')


global data, cfg_de, cfg_cont, cfg_disc, run_id
data = None

def set_cfg(cfg):
    global cfg_de, cfg_cont, cfg_disc, run_id
    cfg_de = cfg.get('density_evolution')
    cfg_cont = cfg_de.get('ga_continuous')
    cfg_disc = cfg_de.get('ga_discrete')
    run_id = cfg.get('run_id')

def to_cdf(pdf):
    """Returns the discrete pdf from a given cdf"""
    return np.cumsum(pdf)


def to_pdf(cdf):
    """Returns the discrete pdf from a given cdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]


def convolution_pad(x, final_size):
    """Returns the given array padded to at least double the final_length to avoid circular convolution, rounded to the nearest power of two"""
    nearest_power = np.ceil(np.log2(final_size))
    padding = int(2**nearest_power - x.size)
    x = np.pad(x, (0, padding))

    return x


def gamma(F, F_grid, G_grid):
    """
    Given a discrete stochastic variable f defined by a cdf with probabilities F for each value on F_grid,
    calculates the cdf of the transformed variable
        g = gamma(f)
    with gamma defined as in [1].
    Because G has to be an equidistant grid for numeric purposes, this creates a small quantization error.
    """
    zero_index = F.size//2
    F_step = abs(F_grid[1] - F_grid[0])

    G0_indices = np.floor(-np.log(np.tanh(G_grid/2)) / F_step)
    G0_indices = np.clip(G0_indices, 0, zero_index - 1).astype(int)
    G0 = 1 - F.take(G0_indices + zero_index)

    G1_indices = np.floor(np.log(np.tanh(G_grid/2)) / F_step)
    G1_indices = np.clip(G1_indices, -zero_index, -1).astype(int)
    G1 = F.take(G1_indices + zero_index)

    return np.stack((G0, G1))


def gamma_inv(G, F_grid, G_grid):
    """
    Given a discrete stochastic variable g defined by a 2 dimensional cdf with probabilities G for each value on (GF(2) x G_grid),
    calculates the cdf of the transformed variable
        f = gamma^-1(g)
    with gamma^-1 defined as the inverse of gamma in [2].
    """
    zero_index = F_grid.size//2
    G_step = abs(G_grid[1] - G_grid[0])

    F_neg_indices = np.floor(-np.log(np.tanh(-F_grid[:zero_index]/2)) / G_step)
    F_neg_indices = np.clip(F_neg_indices, 0, G[1, :].size-1).astype(int)
    F_neg = G[1, :].take(F_neg_indices)

    f_0 = G[1, -1]

    F_pos_indices = np.floor(
        -np.log(np.tanh(F_grid[zero_index+1:]/2)) / G_step)
    F_pos_indices = np.clip(F_pos_indices, 0, G[0, :].size-1).astype(int)
    F_pos = 1 - G[0, :].take(F_pos_indices)

    return np.hstack((F_neg, f_0, F_pos))


# def rho(x, coeffs):
#     """
#     """
#     final_size = (x[0, :].size - 1)*len(coeffs) + 1
#     dx0, dx1 = to_pdf(x[0, :]), to_pdf(x[1, :])
#     x0, x1 = convolution_pad(
#         x[0, :], final_size), convolution_pad(x[1, :], final_size)
#     dx0, dx1 = convolution_pad(
#         dx0, final_size), convolution_pad(dx1, final_size)

#     x0_ft, x1_ft = np.fft.fft(x0), np.fft.fft(x1)
#     dx0_ft, dx1_ft = np.fft.fft(dx0), np.fft.fft(dx1)
#     y0_ft, y1_ft = np.zeros(x0.shape, complex), np.zeros(x1.shape, complex)

#     for coeff in coeffs[1:]:
#         x0_ft, x1_ft = x0_ft*dx0_ft + x1_ft*dx1_ft, \
#             x0_ft*dx1_ft + x1_ft*dx0_ft
#         y0_ft, y1_ft = y0_ft + coeff * x0_ft, y1_ft + coeff * x1_ft

#     y = np.stack((np.abs(np.fft.ifft(y0_ft)[:final_size]), np.abs(
#         np.fft.ifft(y1_ft)[:final_size])))
#     y[0, np.argmax(y[0, :]):] = np.max(y[0, :])
#     y[1, np.argmax(y[1, :]):] = np.max(y[1, :])

#     return y


# def lambd(x, coeffs):
#     """
#     """
#     final_size = (x.size - 1)*len(coeffs) + 1
#     dx = to_pdf(x)

#     x = convolution_pad(x, final_size)
#     dx = convolution_pad(dx, final_size)

#     x_ft = np.fft.fft(x)
#     dx_ft = np.fft.fft(dx)
#     y_ft = np.zeros(x.size, complex)

#     for coeff in coeffs[1:]:
#         x_ft = x_ft * dx_ft
#         y_ft += coeff*x_ft

#     y = np.abs(np.fft.ifft(y_ft)[:final_size])
#     y[np.argmax(y):] = 1
#     return y


# def conv(x, x0):
#     """
#     """
#     dx = to_pdf(x)
#     final_size = x.size + x0.size - 1

#     x0 = convolution_pad(x0, final_size)
#     dx = convolution_pad(dx, final_size)

#     y = np.abs(np.fft.ifft(np.fft.fft(dx)*np.fft.fft(x0))[:final_size])
#     y[np.argmax(y):] = 1
#     return y


def conv(x, x0):
    """
    Naive implementation of conv using sp.convolve, kept for reference.
    """
    x0 = np.pad(x0, (x0.size, x0.size), constant_values=(0, 1))
    x = np.pad(x, (x.size, x.size), constant_values=(0, 1))
    dx = to_pdf(x)

    y = sp.convolve(x0, dx)

    current_size = y.size
    y = y[:np.argmax(y)+1]
    y = np.pad(y, (0, current_size-y.size), constant_values=y.max())

    y = y[y.size//4: -y.size//4]
    return y


def lambd(x, coeffs):
    """
    Naive implementation of lambd using sp.convolve
    """
    x = np.pad(x, (x.size//2, x.size//2), constant_values=(0, 1))
    dx = to_pdf(x)
    final_size = x.size * len(coeffs)

    y = np.zeros(final_size)
    for coeff in coeffs[1:]:
        x = sp.convolve(x, dx)
        current_size = x.size

        x = x[:np.argmax(x)+1]
        padding1 = int(np.ceil((-current_size+final_size)/2))
        padding2 = int(np.floor((-current_size+final_size)/2))
        padding3 = int((current_size-x.size))
        x = np.pad(x, (padding1, padding2 + padding3), constant_values=(0, 1))

        y += coeff*x
        x = x[padding1:-padding2]

    y = y[y.size//4:-y.size//4]
    return y


def rho(x, coeffs):
    """

    """
    dx = np.stack((to_pdf(x[0, :]), to_pdf(x[1, :])))
    final_size = x.size//2 * len(coeffs)

    x0 = x[0, :]
    x1 = x[1, :]
    y = np.zeros((2, final_size))
    for coeff in coeffs[1:]:
        x0, x1 = sp.convolve(x0, dx[0, :]) + sp.convolve(x1, dx[1, :]),\
            sp.convolve(x1, dx[0, :]) + sp.convolve(x0, dx[1, :])
        current_size = x.size//2

        x0 = x0[:np.argmax(x0)+1]
        x1 = x1[:np.argmax(x1)+1]
        x0 = np.pad(x0, (0, final_size-x0.size), constant_values=x0.max())
        x1 = np.pad(x1, (0, final_size-x1.size), constant_values=x1.max())
        y[0, :] += coeff * x0
        y[1, :] += coeff * x1

        x0 = x0[:current_size]
        x1 = x1[:current_size]

    return y


def init_pdf(rber, n_grid=512, llr_max=15):
    """Returns density values of a DMC pdf and its grid in F and G from a look up table."""
    mu1 = -1
    mu2 = 1

    # Load database if not loaded
    global data
    file_path = cfg_de.get("dmc_file")
    if data is None:
        data = np.loadtxt(file_path, delimiter=" ", skiprows=1)

    # Interpolate values on rber
    rber_step = data[1, 0] - data[0, 0]
    rber_index = (rber - data[0, 0])/rber_step
    index0 = int(np.floor(rber_index))
    index1 = int(np.ceil(rber_index))
    values = data[index0] + \
        (data[index1] - data[index0]) * (rber_index - index0)

    # Set values
    sigma = values[1]
    sigma1 = values[2]
    sigma2 = values[3]

    if len(values) == 8:
        thresholds = np.array([values[5]])
        llrs = values[6:]
    elif len(values) == 12:
        thresholds = values[5:8]
        llrs = values[8:]
    else:
        raise Exception("Invalid dmc-file")

    # Calculate probabilities
    p0 = norm.cdf(np.append(thresholds, np.inf), mu2, sigma2)
    p0 = np.ediff1d(p0, to_begin=p0[0])
    #p1 = 1 - norm.cdf(np.array(thresholds), mu1, sigma1)
    #p1 = np.ediff1d(p1, to_begin=p1[0])
    p1 = norm.cdf(np.append(thresholds,np.inf),mu1,sigma1)
    p1 = np.ediff1d(p1, to_begin=p1[0])

    # Symmetrize probabilities
    # Source: https://books.google.se/books?hl=en&lr=&id=ZJrZPObOe60C&oi=fnd&pg=PR13&dq=modern+coding+theory+urbanke&ots=WogyOo3ipu&sig=qKPObOLJNFbt8_4h9NKPd-bzwo4&redir_esc=y#v=onepage&q=modern%20coding%20theory%20urbanke&f=false
    p1 = np.flip(p1)
    p_sym = (p0+p1)/2

    step = 2*llr_max / n_grid
    bins = np.round(-llrs / step) + n_grid//2
    bins = np.clip(bins, 0, n_grid-1)

    x1 = np.linspace(-llr_max, llr_max, n_grid, endpoint=False)
    y = np.zeros(x1.shape)
    np.put(y, bins.astype(int), p_sym)

    f_step = abs(x1[1] - x1[0])
    max_val = -np.log(np.tanh(f_step))
    x2 = np.linspace(0, max_val, n_grid)

    return x1, x2, y


def compute_cbp(pdf, f_grid,is_zero):
    y = np.exp(-f_grid/2)*pdf
    ber = np.sum(pdf[:int(f_grid.size/2)])
    cbp_avg = np.sum(y)
    is_valid = (2*ber <= cbp_avg and cbp_avg <= 2*np.sqrt(ber*(1-ber)))
    message = f"CBP should be bounded by 2*ber<= CBP <= 2*sqrt(ber*(1-ber)). CBP = {cbp_avg}. BER = {ber}."
    #assert is_valid, message
    # if not is_valid:
    #    message = f"CBP should be bounded by 2*ber<= CBP <= 2*sqrt(ber*(1-ber)). CBP = {cbp_avg}. BER = {ber}."
    #    print(message)
    #plt.plot(f_grid, y)
    # plt.show()
    return cbp_avg


def poly_eval(poly, x):
    # sum = 0
    # for i, coeff in enumerate(poly):
    #    sum += coeff*x**i
    #
    # return sum
    return np.sum(poly*x**(np.arange(poly.size)))


def compute_eps(lam, rho, r):
    lam_2 = lam[1]
    rho_prime = rho[1:]*np.arange(1, rho.size)

    rho_prime_1 = poly_eval(rho_prime, 1)
    f1 = lam_2*rho_prime_1*r

    if f1 < 1:
        if np.abs(lam_2-1) > 1e-1:
            f2 = poly_eval(lam, rho_prime_1)*r
            eps = (1-f1)/(f2-f1)
            if eps > r or eps < 0:
                eps = 1
        else:
            eps = -1
    else:
        eps = -1
    return eps


def compute_threshold(pdf, f_grid, lam, rho):
    # Based on "C. Stability" section of:
    # https://arxiv.org/pdf/cs/0509014.pdf
    r = compute_cbp(pdf, f_grid,True)

    eps = compute_eps(lam, rho, r)

    return eps


def symmetric_density_evolution(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, tol,  plot=False, prnt=False):
    global cfg_de
    n_iter = cfg_de.get("de_iter")
    if plot:
        fig, axes = plt.subplots(3, 2)

    # assert np.sum(rho_coeffs[1:]) == 1, "Invalid rho polynom"
    # assert np.sum(lambda_coeffs[1:]) == 1, "Invalid lambda polynom"

    pl = cdf
    p0 = cdf
    pl_old = 0
    i = 0
    diff = np.inf
    error = np.inf
    while True:
        # Compute pdf to determine convergence to zero
        pdf_l = to_pdf(pl)
        pdf_l = np.where(pdf_l > 1e-10, pdf_l, 0) # Numerical thing
        prob_l = np.sum(pdf_l)

        is_valid = (np.abs(1-prob_l)<1e-5)
        message = f"Total probability has to sum to one with tolerance 1e-5. Total probability = {prob_l:.6f}."
        assert is_valid, message

        # Compute correct f_grid
        if i == 0 or i == 1:
            coeff = pl.size/f_grid.size
            start = coeff*f_grid[0]
            step = (f_grid[-1]-f_grid[0])/(f_grid.size-1)
            f_grid_l = np.arange(start, -start, step)

        # Check convergence
        cbp = compute_cbp(pdf_l, f_grid_l,False)
        is_zero = (cbp < tol)
        max_iter = (i == n_iter)
        is_converged = False

        if is_zero:
            error = 0
            if prnt:
                print("Is zero")
        if is_converged and prnt:
            print("Converged")
        if max_iter and prnt:
            print("MAX")
        if is_zero or is_converged or max_iter:
            if error < 0:
                me = "Tjabba"
            if plot:
                plt.show()

            return error

        is_zero = (cbp < tol)
        
        # Perform one iteration of density evolution
        x1 = gamma(pl, f_grid, g_grid)
        x2 = rho(x1, rho_coeffs)
        x3 = gamma_inv(x2, f_grid, g_grid)
        x4 = lambd(x3, lambda_coeffs)
        pl = conv(x4, p0)
        pl[np.argwhere(pl<1e-16)] = 0 # Numerical thing

        diff = sum((pl_old - pl)**2)
        pl_old = pl

        zero_index = pl.size//2
        error = pl[zero_index]

        if plot:
            axes[0, 0].plot(x1[0, :])
            axes[0, 0].plot(x1[1, :])
            axes[0, 1].plot(x2[0, :])
            axes[0, 1].plot(x2[1, :])
            axes[1, 0].plot(x3)
            axes[1, 1].plot(x4)
            axes[2, 0].plot(pl)
            axes[2, 1].scatter(i, error)

        is_converged = False #(diff < float(cfg_de.get("is_converged_tol")))
        is_zero = (error < tol)
        max_iter = (i == n_iter)

        if is_zero:
            error = 0
            if prnt:
                print("Is zero", error)
        if is_converged and prnt:
            print("Converged", error)
        if max_iter and prnt:
            print("MAX", error)
        if is_zero or is_converged or max_iter:
            if plot:
                plt.show()

            return error

        i += 1


def eval(rber, rho_edge, lam_edge):
    plot = False
    n_grid = cfg_de.get("n_grid")

    global X1, X2
    f_grid, g_grid, pdf = init_pdf(rber, X1, X2)
    cdf = to_cdf(pdf)
    if plot:
        plt.plot(f_grid, cdf)
        plt.show()

    eps = compute_threshold(pdf, f_grid, lam_edge, rho_edge)
    if eps == -1:  # will not converge to zero
        result = 1
    elif eps == 1:  # will converge to zero
        result = 0
    else:
        result = symmetric_density_evolution(
            cdf, f_grid, g_grid, rho_edge, lam_edge, tol=eps,  plot=False)

    return result


def bisection_search(min, max, rho_edge, lam_edge):
    global cfg_de
    tol=float(cfg_de.get("bs_tol"))
    
    while max - min > tol:
        x = (min + max)/2
        result = eval(x, rho_edge, lam_edge)

        if result == 0:
            min = x
        elif result > 0:
            max = x
        else:
            error_message = f"The probability from density evolution is negative ({result})."
            raise Exception(error_message)

    return (min + max)/2


def convert_edge_to_node(p_edge):
    p_node = p_edge/(np.arange(1, len(p_edge)+1) *
                     np.sum(1/np.arange(1, len(p_edge)+1)*p_edge))
    return p_node


def set_continuous_params(start_pt, const_mat):
    global x_0, C_c
    x_0 = start_pt
    C_c = const_mat


def best_individual_discrete(best_ind):
    cn_degrees = np.sum(best_ind, 0)
    vn_degrees = np.sum(best_ind, 1)
    rho_node = np.bincount(vn_degrees)[1:]
    lam_node = np.bincount(cn_degrees)[1:]
    rho_node = rho_node/np.sum(rho_node)
    lam_node = lam_node/np.sum(lam_node)
    return lam_node, rho_node


def best_individual_continous(best_ind):
    global x_0, C_c
    x = x_0 + np.matmul(C_c, best_ind)
    dc = cfg_cont.get("dc")
    dv = cfg_cont.get("dv")
    rho_edge = x[:dc]
    lam_edge = x[dc:]
    lam_node = convert_edge_to_node(lam_edge)
    rho_node = convert_edge_to_node(rho_edge)
    return lam_node, rho_node


def save_params():
    algorithm = cfg_de.get("algorithm")
    data = {

        "Np": [cfg_de.get("Np")],
        "generations": [cfg_de.get("generations")],
        "n_grid": [cfg_de.get("n_grid")],
        "algorithm": [algorithm]
    }

    if algorithm == "ga_continuous" or algorithm == "ga_continuous_parallel":
        data_c = {
            "R": [cfg_cont.get("R")],
            "F": [cfg_cont.get("F")],
            "Cr": [cfg_cont.get("Cr")],
            "dv": [cfg_cont.get("dv")],
            "dc": [cfg_cont.get("dc")]
        }
        data.update(data_c)

    if algorithm == "ga_discrete" or algorithm == "ga_discrete_parallel":
        data_d = {
            "n_vn": [cfg_disc.get("n_vn")],
            "n_cn": [cfg_disc.get("n_cn")],
            "p_vertical": [cfg_disc.get("p_vertical")],
            "p_horizontal": [cfg_disc.get("p_horizontal")],
            "p_mutation": [cfg_disc.get("p_mutation")]
        }
        data.update(data_d)

    df = pd.DataFrame.from_dict(data)

    os.makedirs('data/', exist_ok=True)
    fname = "data/" + run_id + "_de_params.csv"
    df.to_csv(fname)


def save_population(population, fitness, generation, best_idx, best_rber, algorithm):
    if algorithm == "discrete":
        lam_node, rho_node = best_individual_discrete(
            population[best_idx, :, :])
    elif algorithm == "continuous":
        lam_node, rho_node = best_individual_continous(population[:, best_idx])

    fname = "data/" + run_id + ".npz"
    with open(fname, 'wb') as f:
        np.savez(f, population=population, fitness=fitness, generation=np.array([generation]), best_idx=np.array(
            [best_idx]), best_rber=np.array([best_rber]),  lam_node=lam_node, rho_node=rho_node)


def log(txt, options):
    fname = "data/log_" + run_id + ".txt"
    with open(fname, options) as f:
        f.write(txt+"\n")

# %%
