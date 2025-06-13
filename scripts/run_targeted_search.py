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
import pandas as pd

from targeted_cws_ng15.models import cw_model_2, gwb_only_model
#from DR3_noise_modeling.utils import get_n_samples, get_initial_sample
#from DR3_noise_modeling.sampler import setup_sampler
from enterprise_extensions.sampler import save_runtime_info
from targeted_cws_ng15.jump_proposal import JumpProposal
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
import targeted_cws_ng15.Dists_Parameters as dists
from targeted_cws_ng15.empirical_distr_new import EmpiricalDistribution2D_v2
from targeted_cws_ng15.utils import get_initial_sample


if __name__ == '__main__':
    
    PARSER = argparse.ArgumentParser()
    
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
    
    # input pulsar to run dropout analysis
    PARSER.add_argument('--dropout', action='store',
                        type=str, dest='dropout', default=None,
                        nargs='?', const=None)
    
    # label orf for the CRN. Currently only supports 'hd' or None
    PARSER.add_argument('-orf', '--orf', action='store',
                        type=str, dest='orf', default=None,
                        nargs='?', const=None)
    
    # flag to turn off the CW signal
    PARSER.add_argument('--no_cw', action='store_true',
                        dest='no_cw', default=False,
                        help='Bool type')
    
    # flag to turn off the GWB signal
    PARSER.add_argument('--no_gwb', action='store_true',
                        dest='no_gwb', default=False,
                        help='Bool type')
    
    # frequencies to include in the GWB
    PARSER.add_argument('--gw_components', action='store',
                        type=int, dest='gw_components', default=None,
                        nargs='?', const=None)
    
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
    
    PARSER.add_argument('-op', '--outdir_path', action='store',
                        type=str, dest='outdir_path', default=os.getcwd())
    
    PARSER.add_argument('--pkldir_path', action='store',
                        type=str, dest='pkldir_path', default=os.getcwd())
    
    # specific ra and dec to use - overrides priors file
    # format is hmsdms: '##h##m##.###s +##d##m##.###s' (+ or -)
    PARSER.add_argument('--ra_dec', action='store',
                        type=str, dest='ra_dec', default=None,
                        nargs='?', const=None)

    PARSER.add_argument('-hm', '--human', action='store',
                        type=str, dest='human', default='blarsen')

    PARSER.add_argument('-ed', '--emp_dist_name', action='store',
                        type=str, dest='emp_dist_name', default=None,
                        nargs='?', const=None)
    
    PARSER.add_argument("-T", dest="T", default=10,
                        type=int, action='store')
    
    PARSER.add_argument("-N", dest="Niter", default=500_000,
                        type=int, action='store')
    
    # this is useful if diagnosing PT issues
    PARSER.add_argument('--write_hot_chains', action='store_true',
                        dest='write_hot_chains', default=False,
                        help='Bool type')
    
    args = PARSER.parse_args()
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    print('kicking off parallel tempering')
    if args.dataset == None:
        raise NameError('Missing dataset')
    else:
        print(f'Dataset: {args.dataset}')
    if args.source_name == None or args.source_name == 'no_cw':
        if args.no_cw or args.source_name  == 'no_cw':
            print('No CW signal - modeling GWB and pulsar noise only')
        else:
            raise NameError('Missing target')
    else:
        print(f'Target: {args.source_name}')
        
    # setup directory for chains
    # head directory (dataset info)
    outdir = f'{args.outdir_path}/data/chains/{args.dataset}'
    if args.droppsr is not None:
        outdir += '_drop'+args.droppsr+'/'
    else:
        outdir += '/'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    # sub directory (target+run info)
    if args.no_cw:
        outdir += 'no_cw'
    else:
        outdir += f'{args.source_name}'
    if args.detection:
        outdir += '_det'
    elif args.upper_limit:
        outdir += '_UL'
    if args.vary_fgw:
        outdir += '_varyfgw'
    if args.all_sky:
        outdir += '_allsky'
    if args.dropout is not None:
        outdir += f'_dropout{args.dropout}'
    elif isinstance(args.ra_dec, str):
        outdir += f'_{args.ra_dec}'
    if args.vary_crn and args.orf == 'hd':
        outdir += '_varygwb'
    elif args.vary_crn:
        outdir += '_varycrn'
    elif args.no_gwb:
        outdir += '_nogwb'
    if args.freespec:
        outdir += '_FS'
    elif args.bpl:
        outdir += '_BPL'
    #outdir += str(args.gw_components)
    if args.fixedpoint:
        outdir += '_fixedpoint'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # check if using job array
    if "SLURM_ARRAY_TASK_ID" in list(os.environ.keys()):
        task_id = os.environ["SLURM_ARRAY_TASK_ID"]
        outdir += f'/{task_id}'
        if not os.path.isdir(outdir):
            try:
                os.mkdir(outdir)
            except Exception as e:
                print('Exception when setting up directories:', e)
                print('Directory probably exists already!')
    print(f'Output to: {outdir}')
    
    # check if using ecorr
    if 'mdc' in args.dataset:
        print('Using MDC data: no ecorr and use TNEQUAD')
        inc_ecorr = False
        tnequad = True
    else:
        inc_ecorr = True
        tnequad = False
            
    # check number of samples
    try:
        with open(f'{outdir}/chain_1.0.txt', 'r') as f:
            n_samples = len(f.readlines())
    except:
        n_samples = 0
    if n_samples < 100:
        print(f'{n_samples} samples so far, setting resume = False!')
        resume = False
    else:
        print(f'{n_samples} samples so far')
        resume = True
    
    # Load PSR objects, depending on dataset
    psrs = []
    # og dataset or dropout psr dataset
    if args.dataset == 'ng15_v1p1' or 'ng15_v1p1_drop' in args.dataset:
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}/'
        pfile += 'v1p1_de440_pint_bipm2019.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    # 12p5yr dataset
    elif args.dataset == 'ng12p5':
        pfile = f'{args.project_path}/../12p5yr_custom_noise_models/psr_pkls/'
        pfile += 'channelized_12yr_v3_partim_py3.pkl'
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
    # sim dataset with Rohan and Gondor
    elif args.dataset == 'ng15_sim_GWB_Rondor':
        pfile = f'{args.pkldir_path}/data/ePSRs/{args.dataset}/enterprise_psrs.pkl'
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
    else:
        pfile = f'{args.project_path}/data/ePSRs/{args.dataset}.pkl'
        with open(pfile, 'rb') as f:
            psrs = pickle.load(f)
    if args.droppsr is not None:
        idx_drop = [i for i, p in enumerate(psrs) if p.name == args.droppsr][0]
        print(f'Dropping PSR {psrs[idx_drop].name}')
        psrs.pop(idx_drop)
    print(f'Loaded {len(psrs)} pulsars from {args.dataset}')
    
    # get kwargs
    if not args.no_cw:
        df = pd.read_csv(f'{args.project_path}/priors/all_targets_info.csv')
        df = df.set_index('Target')
        priors = df.loc[[args.source_name]].iloc[0]
        # overwrite ra and dec on file
        if isinstance(args.ra_dec, str):
            print(f'Overwriting RA and DEC on file to {args.ra_dec}')
            ra = args.ra_dec[:args.ra_dec.index('_')-1]
            dec = args.ra_dec[args.ra_dec.index('_')+1:]
            priors['RA'] = ra
            priors['DEC'] = dec
        
    # CRN params
    if args.vary_crn:
        print(f'Vary CRN params. ORF = {args.orf}')
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
            print('run GWB FS model (no CW)')
            gw_psd = 'spectrum'
            gw_components = 30
        elif args.bpl:
            print('run GWB BPL model (no CW)')
            gw_psd = 'broken_powerlaw'
            gw_components = 30
        else:
            print('run GWB PL model (no CW)')
            gw_psd = 'powerlaw'
            gw_components = 14
        if args.gw_components:
            gw_components = args.gw_components
        pta = gwb_only_model(psrs, noisedict_path, orf=args.orf, gw_psd=gw_psd,
                             gw_components=gw_components, tnequad=tnequad,
                             log10_A_val=log10_A, gamma_val=gamma,
                             fixedpoint=args.fixedpoint, inc_ecorr=inc_ecorr)
    else:
        print('run CW model')
        if args.gw_components:
            gw_components = args.gw_components
        else:
            gw_components = 14
        pta = cw_model_2(psrs, priors, noisedict_path=noisedict_path,
                         psr_distance_path=psr_distance_path, orf=args.orf,
                         gw_components=args.gw_components, tnequad=tnequad,
                         log10_A_val=log10_A, gamma_val=gamma,
                         log10_mc_prior=log10_mc_prior, log10_fgw_prior=args.log10_fgw_prior,
                         vary_fgw=args.vary_fgw, all_sky=args.all_sky,
                         fixedpoint=args.fixedpoint, nogwb=args.no_gwb,
                         dropout=args.dropout, inc_ecorr=inc_ecorr)
    
    # initial sample
    if args.emp_dist_name and not args.emp_dist_name == '15yr_emp_distr':
        try:
            print('Attempting to get initial sample using '
                  f'{args.emp_dist_name}...')
            outdir_label = args.emp_dist_name.replace('_emp_dist', '')
            x0 = get_initial_sample(pta, args.dataset, outdir_label,
                                    args.outdir_path)
        except:
            try:
                print('Try the old directory for an initial sample')
                x0 = get_initial_sample(pta, args.dataset, outdir_label,
                                        args.outdir_path, try_alt_dir=True)
            except Exception as e:
                print('Could not get initial sample! '
                      'Sampling instead from the prior')
                print(f'Exception: {e}')
                x0 = np.hstack([p.sample() for p in pta.params])
    else:
        x0 = np.hstack([p.sample() for p in pta.params])
    ndim = len(x0)
    
    # initial covariance matrix
    '''if args.emp_dist_name and not args.emp_dist_name == '15yr_emp_distr':
        try:
            print('Attempting to get initial covariance using '
                  f'{args.emp_dist_name}...')
            outdirlbl = args.emp_dist_name.replace('_emp_dist', '')
            cov_path = f'{args.outdir_path}/data/chains/{outdirlbl}/0/cov.npy'
            cov = np.load(cov_path)
        except Exception as e:
            print('Could not get initial cov matrix! '
                  'Setting proportional to identity')
            print(f'Exception: {e}')
            cov = np.diag(np.ones(ndim) * 0.1**2)'''
    cov = np.diag(np.ones(ndim) * 0.1**2)
    
    # get groups    
    groups = dists.get_parameter_groups_CAW_target(pta)
    with open(f'{outdir}/groups.txt', 'w') as fi:
        for group in groups:
            line = np.array(pta.param_names)[np.array(group)]
            fi.write("[" + " ".join(line) + "]\n")
        
    # setup sampler
    sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov,
                     groups=groups, outDir=outdir, resume=resume)
    save_runtime_info(pta, outdir=outdir, human=args.human)
    
    # load psr distances
    with open(psr_distance_path, 'rb') as fp:
        dists_file = pickle.load(fp)
    psr_dist = {}
    for p in psrs:
        psr_dist[p.name] = np.array(np.array(dists_file[p.name][:2]))
    
    # get empirical distributions
    if args.emp_dist_name == '15yr_emp_distr':
        emp_dist_file = f'{args.project_path}/empirical_dists/'
        emp_dist_file += f'{args.emp_dist_name}.json'
        with open(emp_dist_file, 'r') as fi:
            emp_dists = json.load(fi)
        emp_distr = []
        for key in emp_dists.keys():
            emp_distr.append(EmpiricalDistribution2D_v2(emp_dists[key]))
    elif args.emp_dist_name:
        emp_dist_file = f'{args.project_path}/empirical_dists/'
        emp_dist_file += f'{args.emp_dist_name}.pkl'
        with open(emp_dist_file, 'rb') as fi:
            emp_distr = pickle.load(fi)
    else:
        emp_dist_file = None
        emp_distr = None
        
    # set up jump proposals
    jp = JumpProposal(pta, empirical_distr=emp_distr)
    # noise prior draws
    print('Adding red noise prior draws...')
    sampler.addProposalToCycle(jp.draw_from_red_prior, 30)
    if 'dm_gp' in jp.snames and not args.fixedpoint:
        print('Adding dm gp prior draws...')
        sampler.addProposalToCycle(jp.draw_from_dm_gp_prior, 10)
    if 'chrom_gp' in jp.snames and not args.fixedpoint:
        print('Adding chrom gp prior draws...')
        sampler.addProposalToCycle(jp.draw_from_chrom_gp_prior, 10)
    if 'exp_1' in jp.snames and not args.fixedpoint:
        print('Adding exp dip prior draws...')
        sampler.addProposalToCycle(jp.draw_from_expdip_prior, 10)
    if 'sw_r2' in jp.snames and not args.fixedpoint:
        print('Adding sw_r2 prior draws...')
        sampler.addProposalToCycle(jp.draw_from_sw_prior, 10)
    if emp_dist_file:
        print('Adding prior draws from empirical distribution...')
        sampler.addProposalToCycle(jp.draw_from_empirical_distr, 30)
    # gwb orior draw
    if args.vary_crn and not args.freespec:
        print('Adding GWB power law prior draws...')
        sampler.addProposalToCycle(jp.draw_from_gwb_log_uniform_distribution,5)
    elif args.vary_crn and args.freespec:
        print('Adding GWB free spectrum prior draws...')
        sampler.addProposalToCycle(jp.draw_from_gw_rho_prior, 25)
    if not args.no_cw:
        jpCW = dists.JumpProposalCW(pta, fgw=10**priors['log10_freq'],
                                    psr_dist=psr_dist,
                                    empirical_distr=emp_distr)
        #pick a cw param & jump
        print('Adding CW prior draws...')
        sampler.addProposalToCycle(jp.draw_from_cw_prior, 20)
        #draw from uniform h
        #sampler.addProposalToCycle(jp.draw_from_cw_log_uniform_distribution,
        #                           10)
        #draw from Mc
        print('Adding Chirp mass prior draws...')
        sampler.addProposalToCycle(jp.draw_from_par_prior(['log10_mc']), 5)
        if args.vary_fgw:
            print('Adding fGW prior draws...')
            sampler.addProposalToCycle(jp.draw_from_par_prior(['log10_fgw']),
                                       10)
        #if 'phase0' in pta.param_names and 'psi' in pta.param_names:
        print('Adding phase/psi reverse jumps...')
        sampler.addProposalToCycle(jpCW.phase_psi_reverse_jump, 1)
        # Pulsar term
        pdist_pars = [p for p in pta.param_names if 'p_dist' in p]
        pphase_pars = [p for p in pta.param_names if 'p_phase' in p]
        #draw from p_dists
        print('Adding p_dist + p_phase prior draws...')
        sampler.addProposalToCycle(jpCW.draw_from_many_par_prior(pdist_pars,
                                                                 'p_dist'), 30)
        #draw from p_phases
        sampler.addProposalToCycle(jpCW.draw_from_many_par_prior(pphase_pars,
                                                                 'p_phase'),30)
        #draw from k_drop
        kdrop_pars = [p for p in pta.param_names if 'k_drop' in p]
        if len(kdrop_pars) > 0:
            print('Adding k_drop prior draws...')
            sampler.addProposalToCycle(jpCW.draw_from_many_par_prior(kdrop_pars, 'k_drop'),30)
        pdist_pars = [p for p in pta.param_names if 'p_dist' in p]
        # sky loc prior draws and custom draws
        if args.all_sky:
            print('Adding sky location prior draws...')
            sampler.addProposalToCycle(jp.draw_from_par_prior(['cos_gwtheta']),
                                       10)
            sampler.addProposalToCycle(jp.draw_from_par_prior(['gwphi']),
                                       10)
        print('Adding CW cyclic par auxiliary jump...')
        sampler.addAuxilaryJump(jpCW.fix_cyclic_pars)
    print('Adding generic prior draws...')
    sampler.addProposalToCycle(jp.draw_from_prior, 10)
    
    # saving model params in chain directory
    with open(outdir+'/model_params.json', 'w') as fout:
        json.dump(pta.param_names, fout, sort_keys=True,
                  indent=4, separators=(',', ': '))
    
    sampler.sample(x0, args.Niter, burn=3_000, thin=args.T, Tskip=10000,
                   SCAMweight=40, AMweight=25, DEweight=25,
                   writeHotChains=args.write_hot_chains)
