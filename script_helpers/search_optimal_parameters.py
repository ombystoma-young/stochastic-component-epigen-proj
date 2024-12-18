import os
import logging

import numpy as np
import pandas as pd

from concurrent.futures import ProcessPoolExecutor
from sklearn.metrics import mean_absolute_error

from simulate_artificial_cohorts import beta_c_next

DATA_FOLDER = 'data'
N_STEPS = 37 * 35
YOUNG_AGE_THRES = 46
OLD_AGE_THRES = 80


def calculate_score_params(params):
    """
    simulare dynamics of DNAm change for given parameters
    """
    es_c = effect_sizes.copy()
    b_c_init = starting_dnam_frac.copy()
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

    # form a package of parameters for estimation
    gammas = np.arange(1, 11, 0.5)
    sigmas = np.arange(0.0001, 0.0030, 0.00015)
    package = zip(np.repeat(gammas, sigmas.shape[0]), np.tile(sigmas, gammas.shape[0]))

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

    # run the function in parallel
    with ProcessPoolExecutor(max_workers=n_threads) as executor:
        results = list(executor.map(calculate_score_params, package))

    with open(os.path.join(DATA_FOLDER, 'mae_res.tsv'), 'wt') as out_f:
        for res in results:
            out_f.write('\t'.join([str(i) for i in res]))
            out_f.write('\n')
