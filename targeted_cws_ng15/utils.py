#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 16:04:17 2024

@author: bjornlarsen
"""

import json, os
import la_forge.core as co
import numpy as np
from enterprise.signals import parameter

def get_prior_log_weights(chain, prior1=parameter.LinearExpPrior,
                          prior0=parameter.UniformPrior,
                          prior1_args=[7,11], prior0_args=[7,11]):
    '''
    chain: MCMC samples over the parameter you are reweighting (e.g., log10_mc)
    prior1: Prior function corresponding to your new model.
        default: parameter.LinearExpPrior for an "upper limit run"
    prior0: Prior function corresponding to your old model. This should match
        what was used to compute the chain.
        default: parameter.UniformPrior for a "detection run"
    prior1_args: arguments for the new prior function. This should match
        what was used to compute the chain.
    prior0_args: arguments for the old prior function. This should match
        what was used to compute the chain.
    '''
    lnp1 = np.log(prior1(chain, *prior1_args))
    lnp0 = np.log(prior0(chain, *prior0_args))
    return lnp1 - lnp0

def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights. From stackexchange
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def get_initial_sample(pta, dataset, outdir_label, project_path,
                       pt_chains=True, try_alt_dir=False):
    """
    Get initial sample as MLV from a previous run

    Parameters
    ----------
    psrname : string, name of pulsar
    model_label : string, name of model from which to take a sample

    Returns
    -------
    Dictionary of param names and sample value
    """
    chaindir = f'{project_path}/data/chains/{dataset}'
    if try_alt_dir:
        chaindir = f'{chaindir}_old_Tspans/{outdir_label}'
    else:
        chaindir = f'{chaindir}/{outdir_label}'
    try:
        with open(chaindir+'/model_params.json' , 'r') as fin:
            model_params = json.load(fin)
    except:
        with open(chaindir+'/0/model_params.json' , 'r') as fin:
            model_params = json.load(fin)
    try:
        corepath = f'{chaindir}/core.h5'
        c = co.Core(corepath=corepath, params=model_params,
                    pt_chains=pt_chains)
    except:
        if pt_chains:
            chain_file = 'chain_1.0.txt'
        else:
            chain_file = 'chain_1.txt'
        if not os.path.isfile(f'{chaindir}/{chain_file}'):
            chaindir += '/1'
        c = co.Core(chaindir=chaindir, params=model_params,
                    pt_chains=pt_chains)
    param_dict = c.get_map_dict()
    x0 = []
    rand_sample_pnames = []
    for i, pname in enumerate(pta.param_names):
        if pname in list(param_dict.keys()):
            x0.append(param_dict[pname])
        else:
            #print(f'no sample for {pname} in {model_label} chain')
            rand_sample_pnames.append(pname)
            x0.append(pta.params[i].sample())
        if 'log10_rho' in pname:
            # pta.params is not the same size as pta.param_names if using FS
            i -= 1
    if len(rand_sample_pnames) > 0:
        print(f'no initial samples for params {rand_sample_pnames}')
    return np.hstack(x0)