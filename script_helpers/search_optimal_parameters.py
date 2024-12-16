import os
import logging

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from scipy.stats import truncnorm, norm
from sklearn.metrics import mean_absolute_error


DATA_FOLDER = 'data'
RANDOMSTATE = 42
N_STEPS = 37 * 35


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
    return np.random.choice([0, 1], 1, True, [1-p, p])[0]


def r_c(sigma):
    return truncnorm.rvs(a=0, b=np.inf, loc=0, scale=sigma)


def x_c(beta_c):
    return norm.ppf(beta_c)


def x_c_next(x_c_, sigma, effsize_c):
    return x_c_ + np.sign(effsize_c) * r_c(sigma)


@np.vectorize
def beta_c_next(effsize_c, gamma, beta_c, sigma):
    ch_st = change_state_status(effsize_c, gamma)
    if ch_st == 1:
        x_c_val = x_c(beta_c)
        x_c_next_val = x_c_next(x_c_val, sigma, effsize_c)
        return norm.cdf(x_c_next_val)
    else:
        return beta_c


def calculate_score_params(params):
    es_c = effect_sizes.copy()
    b_c_init = starnting_dnam_frac.copy()
    b_c_finit = old_dnam_avg.copy()
    gamma, sigma = params

    logging.debug(f'Simulation for {gamma=} {sigma=}')

    # first step depends on avgDNAm(young)
    beta_c_cur = beta_c_next(effsize_c=es_c, gamma=gamma, beta_c=b_c_init, sigma=sigma)

    # next steps depend on previous states
    for _ in range(N_STEPS - 1):
        beta_c_cur = beta_c_next(effsize_c=es_c, gamma=gamma, beta_c=beta_c_cur, sigma=sigma)

    # calculate mae for obtained DNAms
    pred_es = beta_c_cur - b_c_init
    mae = mean_absolute_error(es_c, pred_es)
    return gamma, sigma, mae
    

if __name__ == '__main__':
    logging.basicConfig(format='%(levelname)s:%(asctime)s:%(message)s', level=logging.DEBUG)

    dataset = 'CD14'
    dnam_fracs_path = os.path.join(DATA_FOLDER, f'{dataset}_DNAm_frac.pkl')
    dnam_age_path = os.path.join(DATA_FOLDER, f'{dataset}_sample.pkl')
    n_threads = 20

    dnam_fracs = pd.read_pickle(dnam_fracs_path)
    sample_age = pd.read_pickle(dnam_age_path).rename(columns={'Sample_title': 'sample'}).set_index('sample')

    dnam_fracs_T = dnam_fracs.T
    dnam_df = dnam_fracs_T.join(sample_age)


    young_dnam_df = dnam_df.query('age < 46')
    old_dnam_df = dnam_df.query('age > 80')

    # calculate average DNAm fractions
    young_dnam_avg = young_dnam_df.mean()[:-1]
    old_dnam_avg = old_dnam_df.mean()[:-1]

    # save the average DNAm fraction of youngest as the starting DNAm fraction 
    starnting_dnam_frac = young_dnam_avg.copy()

    # calculate an effect size as a difference between young and old
    effect_sizes = old_dnam_avg - young_dnam_avg

    # form a package for estimation
    gammas = np.arange(1, 11, 0.5)
    sigmas = np.arange(0.0001, 0.0030, 0.00015)

    package = zip(np.repeat(gammas, sigmas.shape[0]), np.tile(sigmas, gammas.shape[0]))

    # run the function in parallel
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        results = list(executor.map(calculate_score_params, package))
    
    with open(os.path.join(DATA_FOLDER, 'mae_res.tsv'), 'wt') as out_f:
        for res in results:
            out_f.write('\t'.join([str(i) for i in res]))
            out_f.write('\n')
