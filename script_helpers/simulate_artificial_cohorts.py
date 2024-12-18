import os
import logging

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from scipy.stats import truncnorm, norm


DATA_FOLDER = 'data'
N_STEPS_PER_AGE = 35
YOUNG_AGE_THRES = 46
OLD_AGE_THRES = 80
N_SAMPLES_PER_AGE = 5


@np.vectorize
def p_c(effsize, gamma):
    """
    calculate the probability of methylation state for particular effect size of CpG (effsize_c) and 
    global probability of a DNAm change (gamma)
    """
    return 1 - np.exp(-gamma * np.abs(effsize))


@np.vectorize
def change_state_status(effsize_c, gamma):
    """
    emits the state of methylation, based on defined effect size of CpG (effsize_c) and 
    global probability of a DNAm change (gamma)
    """
    p = p_c(effsize_c, gamma)
    return np.random.choice([0, 1], replace=True, p=[1-p, p])


def r_c(sigma):
    """
    emits random deviation from truncated (one-sided) 
    normal distribution with m=0, sigma=sigma
    """
    return truncnorm.rvs(a=0, b=np.inf, loc=0, scale=sigma)


def x_c(beta_c):
    """
    converts beta-values of DNAm to norm quantiles
    """
    return norm.ppf(beta_c)


def x_c_next(x_c_, sigma, effsize_c):
    """
    calculates DNAm level in inversed quntiles 
    based on size of sigma and sign of effect size
    """
    return x_c_ + np.sign(effsize_c) * r_c(sigma)


@np.vectorize
def beta_c_next(effsize_c, gamma, beta_c, sigma):
    """
    calculates DNAm level for next step of simulation
    """
    ch_st = change_state_status(effsize_c, gamma)
    if ch_st == 1:
        x_c_val = x_c(beta_c)
        x_c_next_val = x_c_next(x_c_val, sigma, effsize_c)
        return norm.cdf(x_c_next_val)
    else:
        return beta_c


def simulate_dnam(age_diff):
    """
    simulare dynamics of DNAm change for given age
    """
    es_c = effect_sizes.copy()
    b_c_init = starting_dnam_frac.copy()
    if age_diff == 0:  # for age equal to initial state just small differences
        n_steps = np.random.randint(low=1, high=N_STEPS_PER_AGE)
    else:
        n_steps = age_diff * N_STEPS_PER_AGE

    logging.debug(f'Simulation for {age_diff=}')

    # first step depends on avgDNAm(young)
    beta_c_cur = beta_c_next(effsize_c=es_c, gamma=gamma, beta_c=b_c_init, sigma=sigma)

    # next steps depend on previous states
    for _ in range(n_steps - 1):
        beta_c_cur = beta_c_next(effsize_c=es_c, gamma=gamma, beta_c=beta_c_cur, sigma=sigma)

    return np.append(beta_c_cur, age_diff)


if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(asctime)s:%(message)s', level=logging.DEBUG)

    dataset = 'CD14'
    dnam_fracs_path = os.path.join(DATA_FOLDER, f'{dataset}_DNAm_frac.pkl')
    dnam_age_path = os.path.join(DATA_FOLDER, f'{dataset}_sample.pkl')

    n_threads = 20

    # define optimal gamma and sigma
    gamma = 9
    sigma = 0.0004


    dnam_fracs = pd.read_pickle(dnam_fracs_path)
    sample_age = pd.read_pickle(dnam_age_path).rename(columns={'Sample_title': 'sample'}).set_index('sample')

    dnam_fracs_T = dnam_fracs.T
    dnam_df = dnam_fracs_T.join(sample_age)


    young_dnam_df = dnam_df.query('age < @YOUNG_AGE_THRES')
    old_dnam_df = dnam_df.query('age > @OLD_AGE_THRES')

    # calculate average DNAm fractions
    young_dnam_avg = young_dnam_df.mean()[:-1]
    old_dnam_avg = old_dnam_df.mean()[:-1]

    # save the average DNAm fraction of youngest as the starting DNAm fraction 
    starting_dnam_frac = young_dnam_avg.copy()

    # calculate an effect size as a difference between young and old
    effect_sizes = old_dnam_avg - young_dnam_avg

    # form a package for estimation
    for cohort in ['train', 'test', 'validation']:
        age_differences = np.arange(0, 39, 1)

        # run the function in parallel
        with ProcessPoolExecutor(max_workers=n_threads) as executor:
            results = list(executor.map(simulate_dnam, np.tile(age_differences, N_SAMPLES_PER_AGE)))
        # write results into a pickle file
        pd.DataFrame(results).to_pickle(os.path.join(DATA_FOLDER, f'simulated_cohort_{cohort}.pkl'))
        
        logging.info(f'Done simulation of {cohort}')
