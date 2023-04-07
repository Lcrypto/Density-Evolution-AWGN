"""
This file implements the full asymmetric density evolution algorithm for optimisng on the discritized A-BIAWGN channel.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf
https://arxiv.org/pdf/2001.01249.pdf

"""
import ga_discrete as ga_d
import ga_continuous as ga_c
#import ga_discrete_parallel as ga_d_p
#import ga_continuous_parallel as ga_c_p
import matplotlib.pyplot as plt
import numpy as np
import time
from config import cfg


np.seterr(divide='ignore')
cfg_de = cfg.get('density_evolution')


def main():
    algorithm = cfg_de.get("algorithm")
    if algorithm == "ga_continuous":
        ga_c.ga_continous()
    elif algorithm == "ga_discrete":
        ga_d.ga_discrete()
    elif algorithm == "ga_continuous_parallel":
        ga_c_p.ga_continous_parallel()
    elif algorithm == "ga_discrete_parallel":
        ga_d_p.ga_discrete_parallel()
    else:
        raise Exception("No valid algorithm chosen.")


if __name__ == "__main__":
    main()
