#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 17:02:10 2022

@author: bjornlarsen
"""

import numpy as np
import argparse, json, os, glob, pickle
import h5pulsar
from mpi4py import MPI

from targeted_cws_ng15.models import cw_model_2, gwb_only_model
#from DR3_noise_modeling.utils import get_n_samples, get_initial_sample
#from DR3_noise_modeling.sampler import setup_sampler
from enterprise_extensions.sampler import save_runtime_info
from targeted_cws_ng15.jump_proposal import JumpProposal
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
import targeted_cws_ng15.Dists_Parameters as dists
from targeted_cws_ng15.empirical_distr_new import EmpiricalDistribution2D_v2
from targeted_cws_ng15.utils import get_initial_sample
from la_forge import core as co
import multiprocessing as mp
import time


# reweight a GWB or CW run to use HD instead of CRN ORF


PARSER = argparse.ArgumentParser()

# label path to directory where la_forge core should be loaded and results output
PARSER.add_argument('--coredir', action='store',
                    type=str, default=None)

# label source for loading priors
PARSER.add_argument('-sn', '--source_name', action='store',
                    type=str, dest='source_name', default=None,
                    nargs='?', const=None)

# label dataset to run analysis on
PARSER.add_argument('-d', '--dataset', action='store',
                    type=str, dest='dataset', default=None)

# input pulsars to drop from dataset
PARSER.add_argument('--droppsr', action='store',
                    type=str, dest='droppsr', default=None,
                    nargs='?', const=None)

# flag if detection run (log-uniform mass priors)
# priority over upper limit
PARSER.add_argument('--detection', action='store_true',
                    dest='detection', default=False,
                    help='Bool type')

# flag if upper limit run (uniform mass priors)
PARSER.add_argument('--upper_limit', action='store_true',
                    dest='upper_limit', default=False,
                    help='Bool type')

# flag to vary CRN (vary gamma and log10_A w/ log-uniform prior)
PARSER.add_argument('--vary_crn', action='store_true',
                    dest='vary_crn', default=False,
                    help='Bool type')

# flag to vary fGW
PARSER.add_argument('--vary_fgw', action='store_true',
                    dest='vary_fgw', default=False,
                    help='Bool type')

# flag to set prior used for varied fgw (default log-normal)
PARSER.add_argument('--log10_fgw_prior', action='store',
                    type=str, dest='log10_fgw_prior', default='normal')

# flag to vary cos_gwtheta and gwphi (uniform priors)
PARSER.add_argument('--all_sky', action='store_true',
                    dest='all_sky', default=False,
                    help='Bool type')

# flag to turn off the CW signal
PARSER.add_argument('--no_cw', action='store_true',
                    dest='no_cw', default=False,
                    help='Bool type')

# flag to turn off the GWB signal
PARSER.add_argument('--no_gwb', action='store_true',
                    dest='no_gwb', default=False,
                    help='Bool type')

# flag to run a GWB broken power law
PARSER.add_argument('--bpl', action='store_true',
                    dest='bpl', default=False,
                    help='Bool type')

# flag to run a GWB free spectrum
PARSER.add_argument('--freespec', action='store_true',
                    dest='freespec', default=False,
                    help='Bool type')

# flag to fix all custom chromatic noise params
PARSER.add_argument('--fixedpoint', action='store_true',
                    dest='fixedpoint', default=False,
                    help='Bool type')

# label noise dictionary
PARSER.add_argument('-nd', '--noisedict', action='store',
                    type=str, dest='noisedict', default=None)

# label psr distance pickle
PARSER.add_argument('-pd', '--psrdists', action='store',
                    type=str, dest='psrdists', default=None)

PARSER.add_argument('-pp', '--project_path', action='store',
                    type=str, dest='project_path', default=os.getcwd())

# specific ra and dec to use - overrides priors file
# format is hmsdms: '##h##m##.###s +##d##m##.###s' (+ or -)
PARSER.add_argument('--ra_dec', action='store',
                    type=str, dest='ra_dec', default=None,
                    nargs='?', const=None)

# target number of samples to reweight
PARSER.add_argument("-N", dest="Niter", default=100_000,
                    type=int, action='store')

# number of simultaneous processes
PARSER.add_argument("--Ncpus", default=1,
                    type=int, action='store')

# flag to restart the likelihood reweighting run
PARSER.add_argument('--restart', action='store_true',
                    default=False, help='Bool type')

args = PARSER.parse_args()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

print(f'kicking off likelihood reweighting (Ncpus = {args.Ncpus})')
if args.dataset == None:
    raise NameError('Missing dataset')
elif args.coredir == None:
    raise NameError('Missing chains for reweighting')
else:
    print(f'Dataset: {args.dataset}')
    print(f'Corepath: {args.coredir}/core.h5')
if args.source_name == None or args.source_name == 'no_cw':
    if args.no_cw or args.source_name == 'no_cw':
        print('No CW signal - modeling GWB and pulsar noise only')
    else:
        raise NameError('Missing target')
else:
    print(f'Target: {args.source_name}')

# setup directory for chains
# head directory (dataset info)
outdir = f'{args.coredir}/HD_reweighting/'
if not os.path.isdir(outdir):
    os.mkdir(outdir)
print(f'Output to: {outdir}')

# load samples
core_crn = co.Core(corepath=f'{args.coredir}/core.h5', burn=0)

# determine thinning
Ns = len(core_crn('lnpost'))
thin = Ns//args.Niter
print(f'Aiming for {args.Niter} reweighted samples')
print(f'Starting from {Ns} reweighted samples')
print(f'Thinning by {thin}')

if os.path.isfile(f'{outdir}/pta.pkl'):
    # save time making the PTA
    with open(f'{outdir}/pta.pkl','rb') as f:
        pta = cloudpickle.load(f)
else:
    # Load PSR objects, depending on dataset
    psrs = []
    # og dataset or dropout psr dataset
    if args.dataset == 'ng15_v1p1' or 'ng15_v1p1_drop' in args.dataset:
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'v1p1_de440_pint_bipm2019.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # custom J1713 noise dataset
    elif args.dataset == 'ng15_v1p1_customJ1713':
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'v1p1_de440_pint_bipm2019_customJ1713.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # custom 6 psr dataset
    elif args.dataset == 'ng15_v1p1_custom6':
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'v1p1_de440_pint_bipm2019_custom6.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # sim dataset with GWB
    elif args.dataset == 'ng15_sim_GWB':
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'psrs_GWB.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # sim dataset with GWB + CW
    elif (args.dataset == 'ng15_sim_GWB_HS1630' or
          args.dataset == 'ng15_sim_GWB_NGC3115' or
          args.dataset == 'ng15_sim_kneeGWB_NGC3115'):
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'psrs_GWB_CW.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # sim dataset with CW
    elif args.dataset == 'ng15_sim_NGC3115':
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'psrs_CW.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # sim dataset with noise only
    elif args.dataset == 'ng15_sim_noiseonly':
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'psrs_noise.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # load from hdf5 - these are old so probably shouldn't be used
    else:
        print('loading hdf5 pulsars dataset (are you sure?)')
        datastr = f'{args.project_path}/data/ePSRs/{args.dataset}_hdf5/*.hdf5'
        for hdf5_file in glob.glob(datastr):
            psrs.append(h5pulsar.FilePulsar(hdf5_file))
    if args.droppsr is not None:
        idx_drop = [i for i, p in enumerate(psrs) if p.name == args.droppsr][0]
        print(f'Dropping PSR {psrs[idx_drop].name}')
        psrs.pop(idx_drop)
    print(f'Loaded {len(psrs)} pulsars from {args.dataset}')

    # get kwargs
    if not args.no_cw:
        fname = f'{args.project_path}/priors/{args.source_name}_priors.json'
        with open(fname, 'r') as fin:
            priors = json.load(fin)
        # overwrite ra and dec on file
        if isinstance(args.ra_dec, str):
            print(f'Overwriting RA and DEC on file to {args.ra_dec}')
            ra = args.ra_dec[:args.ra_dec.index('_')-1]
            dec = args.ra_dec[args.ra_dec.index('_')+1:]
            priors['RA'] = ra
            priors['DEC'] = dec

    # CRN params
    if args.vary_crn:
        print(f'Vary CRN params')
        log10_A = None
        gamma = None
    elif args.no_gwb:
        print('No CRN or GWB model')
        log10_A = None
        gamma = None
    else:
        print(f'Using fixed CRN params. ORF = {args.orf}')
        log10_A = np.log10(6.4e-15)
        gamma = 3.2

    # Upper limit?
    if not args.no_cw:
        if args.detection:
            print('Using log-uniform priors on chirp mass (detection)')
            log10_mc_prior = 'uniform'
        elif args.upper_limit:
            print('Using uniform priors on chirp mass (upper limit)')
            log10_mc_prior = 'linearexp'
        else:
            print('Using Gaussian priors for chirp mass (astro)')
            log10_mc_prior = 'normal'

    # set up PTA object
    noisedict_path = f'{args.project_path}/noise_dicts/{args.noisedict}.json'
    psr_distance_path=f'{args.project_path}/psr_distances/{args.psrdists}.pkl'
    if args.no_cw:
        if args.freespec:
            print('run GWB FS model (no CW, HD ORF)')
            gw_psd = 'spectrum'
            gw_components = 30
        elif args.bpl:
            print('run GWB BPL model (no CW, HD ORF)')
            gw_psd = 'broken_powerlaw'
            gw_components = 30
        else:
            print('run GWB PL model (no CW, HD ORF)')
            gw_psd = 'powerlaw'
            gw_components = 14
        hd_pta = gwb_only_model(psrs, noisedict_path, orf='hd', gw_psd=gw_psd,
                                gw_components=gw_components,
                                log10_A_val=log10_A, gamma_val=gamma,
                                fixedpoint=args.fixedpoint)
    else:
        print('run CW model (GWB w/ HD ORF)')
        hd_pta = cw_model_2(psrs, priors, noisedict_path=noisedict_path,
                            psr_distance_path=psr_distance_path, orf='hd',
                            log10_A_val=log10_A, gamma_val=gamma,
                            log10_mc_prior=log10_mc_prior, log10_fgw_prior=args.log10_fgw_prior,
                            vary_fgw=args.vary_fgw, all_sky=args.all_sky,
                            fixedpoint=args.fixedpoint, nogwb=args.no_gwb)

    with open(f'{outdir}/pta_summary.txt','w') as f:
        f.write(hd_pta.summary())

    with open(f'{outdir}/pta.pkl','wb') as f:
        cloudpickle.dump(hd_pta, f)
    
# time initial likelihood
#params = {p.name:p.sample() for p in hd_pta.params}
#t0 = time.time()
#lnlike = hd_pta.get_lnlikelihood(params)
#print(f'lnlike = {lnlike}')
#print(f'time = {time.time() - t0} s')

# ----------------------------
# SAMPLER SETUP
# ----------------------------

# resuming run
if os.path.isfile(f'{outdir}/ln_weights.txt') and not args.restart:
    # load in idxs + values so far
    chain_idxs = np.loadtxt(f'{outdir}/chain_idxs.txt') # which values of the chain to reference
    ln_weights = np.loadtxt(f'{outdir}/ln_weights.txt') # computed ln_weights for each sample in the chain
    new_array_idxs = np.arange(len(ln_weights))#reference to keep track of where in the chain_idxs and ln_weights array we are at
    # print stuff
    if not np.any(np.isnan(ln_weights)):
        print('You have already finished processing this run! (no more nans in output ln_weights file)')
        raise ValueError('You have already finished processing this run! (no more nans in output ln_weights file)')
    print(f'resuming run ({len(chain_idxs[~np.isnan(ln_weights)])}/{len(chain_idxs)} processed)')
    
# starting new run
else:
    # setup arrays of chain idxs + ln weight values + ln weight idxs
    chain_idxs = np.arange(0,Ns,thin) # which values of the chain to reference
    ln_weights = np.zeros(len(chain_idxs))*np.nan # computed ln_weights for each sample in the chain
    new_array_idxs = np.arange(len(chain_idxs))#reference to keep track of where in the chain_idxs and ln_weights array we are at
    # save idxs to process
    np.savetxt(f'{outdir}/chain_idxs.txt', chain_idxs)
    # print stuff
    print(f'Starting new run (0/{len(chain_idxs)} processed)')
    if os.path.isfile(f'{outdir}/ln_weights.txt'):
        print("Overwriting the existing file!!!")
    
# setup mask to filter out processed values (only process samples which are still nan in the ln_weights array)
mask = np.isnan(ln_weights)

# setup[ parallel processing
print('setting up parallel processing')
manager = mp.Manager()
lock = manager.Lock() # used to ensure that writing to the shared array is thread-safe

# Turn ln_weights into a shared array all workers can write to
ln_weights_shared = mp.Array('d', len(ln_weights))  # 'd' for double precision
for i in range(len(ln_weights)):
    ln_weights_shared[i] = ln_weights[i]

# Create a pool of workers
pool = mp.Pool(args.Ncpus)

# Create try/except block to ensure resources are properly released after making pool
try:
    # Split chain_idxs and new_array_idxs into chunks
    # TODO: Improve the code readability by only making new_array_idx_chunks, which can be used to get chain_idx values
    Npoints = len(ln_weights[mask])
    chunk_size = int(np.ceil(Npoints/args.Ncpus))
    chain_idx_chunks = [chain_idxs[mask][i:i + chunk_size] for i in range(0, Npoints, chunk_size)]
    new_array_idx_chunks = [new_array_idxs[mask][i:i + chunk_size] for i in range(0, Npoints, chunk_size)]

    def process_chunk(chain_idxs, new_array_idxs, ln_weights_shared, lock):
        print('iteration!')
        for chain_idx, new_array_idx in zip(chain_idxs, new_array_idxs):
            # get sample
            params = {p:core_crn(p)[chain_idx] for p in core_crn.params}
            # rename GWB param for HD model
            params['gwb_log10_A'] = params['crn_log10_A']
            # get ln_weight for given sample
            ln_weight = hd_pta.get_lnlikelihood(params) - core_crn('lnlike')[chain_idx]
            with lock: # ensures array updating occurs smoothly across processes
                ln_weights_shared[new_array_idx] = ln_weight

    # Start processing chunks asynchronously
    for chain_idx_chunk, new_array_idx_chunk in zip(chain_idx_chunks, new_array_idx_chunks):
        print('applying async')
        pool.apply_async(process_chunk, args=(chain_idx_chunk, new_array_idx_chunk,
                                              ln_weights_shared, lock))

    start_time = time.time()

    #while np.any(np.isnan(ln_weights_shared)):

    # wait a minute before each update
    time.sleep(20)

    # compute stats
    Nvals = np.count_nonzero(~np.isnan(ln_weights_shared))
    if Nvals > 0:
        mean_wgt = np.nanmean(np.exp(ln_weights_shared))
        sigma_wgt = np.nanstd(np.exp(ln_weights_shared))
        neff = Nvals/(1 + (sigma_wgt/mean_wgt)**2)
        eff = neff/Nvals
        sigma_B = sigma_wgt/np.sqrt(neff)

        # print progress
        postfix_str = (f'{Nvals}/{Npoints} | B = {np.round(mean_wgt,decimals=3)} +- {np.round(sigma_B,decimals=3)} | '
                       f'neff = {np.round(neff,decimals=3)} | efficiency = {np.round(eff,decimals=3)}')
        print(postfix_str)

        # update array
        np.savetxt(f'{outdir}/ln_weights.txt', ln_weights_shared)
    else:
        print('No non-nan values in array yet')
        
    time.sleep(20)

    # Save final results
    pool.close()
    pool.join()

    print('FINISHED!!')
    Nvals = np.count_nonzero(~np.isnan(ln_weights_shared))
    mean_wgt = np.nanmean(np.exp(ln_weights_shared))
    sigma_wgt = np.nanstd(np.exp(ln_weights_shared))
    neff = Nvals/(1 + (sigma_wgt/mean_wgt)**2)
    eff = neff/Nvals
    sigma_B = sigma_wgt/np.sqrt(neff)
    postfix_str = (f'{Nvals}/{Npoints} | B = {np.round(mean_wgt,decimals=3)} +- {np.round(sigma_B,decimals=3)} | '
                   f'neff = {np.round(neff,decimals=3)} | efficiency = {np.round(eff,decimals=3)}')
    print(postfix_str)
    np.savetxt(f'{outdir}/ln_weights.txt', ln_weights_shared)
    print(f'Total time: {time.time() - start_time}')
except Exception as e:
    print(f"An error occurred: {e}")
    pool.terminate()
    pool.join()
    raise e
