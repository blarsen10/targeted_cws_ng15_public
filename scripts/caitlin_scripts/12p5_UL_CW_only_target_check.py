#!/usr/bin/env python

from __future__ import division
import numpy as np
import json
#import glob
#import matplotlib.pyplot as plt
from astropy import units as u
#import astropy.constants as c
from astropy.coordinates import SkyCoord
#import scipy as sp
import os
import pickle as pickle

from enterprise.signals import parameter
#from enterprise.pulsar import Pulsar
from enterprise.signals import selections
from enterprise.signals import signal_base
from enterprise.signals import white_signals
from enterprise.signals import gp_signals
from enterprise.signals import deterministic_signals
#import enterprise.constants as const
from enterprise.signals import utils


import new_delays_nd as nd

#from enterprise_extensions.models import t_process
#from enterprise_extensions.models import InvGamma, InvGammaPrior, InvGammaSampler

#from enterprise_extensions import model_utils

#from enterprise_extensions.model_utils import HyperModel


from enterprise_extensions.deterministic import CWSignal#, cw_delay
from jump_proposal import JumpProposal


from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc
#from corner import corner, quantile

#from bound_norm_py import BoundedNormal

# import targeted_functions as fns

import argparse
import Dists_Parameters as dists


parser = argparse.ArgumentParser(description = "CW run at frequency: ")
# 
# parser.add_argument("-freq_index", required = True, type = int, help = "Index of the frequency in the file you'd like to search")
# ##run a batch of all values in your freq file!
# parser.add_argument("-batch_number", required = True, type = str, help = "Index of the frequency in the file you'd like to search")
##this is the file number you signed up for!
parser.add_argument("-N_iter", required = True, type = float, help = "number of iterations", default = 1e6)
##this is the file number you signed up for!


args = parser.parse_args()

# num = args.batch_number
# freq_ind = args.freq_index
N_iter = args.N_iter


#use a reasonable number for your machine's walltime! I use 1e6 for a week's walltime, then restart.
# N = int(N_iter)


# #update this to your system!
# freq_file = '/scratch/caw0057/nano12p5_gwb/cw_12p5/detection_scripts/12p5_cw_freqs_'+num +'.txt'
# 
# freqs = np.loadtxt(freq_file)
# 
# freq = freqs[freq_ind]

#update this to your system!
chaindir = '/scratch/caw0057/nano12p5_gwb/cw_12p5/target_scripts_e_rn_emp_nd/chains/'
outdir = chaindir+'/3c66b_UL_CW_only/'

def run_ul(N):
    #update this to your system!
    data = 'channelized_12yr_v3_partim_DE438.pkl'
    datadir = '/scratch/caw0057/nano12p5/'
    
    #update this to your system!
    initdir = '/ocean/projects/phy210013p/cawitt/nano12p5/'
    empirical_distr = initdir + 'new_12p5yr_rn_emp_distr.pkl'
    # empirical_distr = initdir + '12yr_emp_dist_RNonly_py3.pkl'
    
    
    with open(datadir+data, 'rb') as fp:
        psrs = pickle.load(fp)
    
    
    
    #update this to your system!
    noise = 'channelized_12p5yr_v3_full_noisedict.json'
    with open(datadir+noise, 'r') as fp:
        noisedict = json.load(fp)
        
        
    c = SkyCoord('02h23m11.4112s', '+42d59m31.384s', frame='icrs')
    freq = 60.4e-9
    log_f = np.log10(freq)
    log_dist = np.log10(85.8)
    
        
    #everything below here should be fine!
    #######################################################################
    
    
    tmin = [p.toas.min() for p in psrs]
    tmax = [p.toas.max() for p in psrs]
    Tspan = np.max(tmax) - np.min(tmin)
    #for key in list(noisedict.keys()):
    #    if any(no in key for no in nos):
    #        del(noisedict[key])
    
    
    efac = parameter.Constant()
    equad = parameter.Constant()
    ecorr = parameter.Constant()
    
    # define selection by observing backend
    selection = selections.Selection(selections.by_backend)
    
    # define white noise signals
    ef = white_signals.MeasurementNoise(efac=efac, selection=selection)
    eq = white_signals.EquadNoise(log10_equad=equad, selection=selection)
    ec = white_signals.EcorrKernelNoise(log10_ecorr=ecorr, selection=selection)
    
    log10_A = parameter.Uniform(-20, -11)
    gamma = parameter.Uniform(0, 7)
    
    # define powerlaw PSD and red noise signal
    pl = utils.powerlaw(log10_A=log10_A, gamma=gamma)
    rn = gp_signals.FourierBasisGP(pl, components=30)
    
    
    
    cos_gwtheta = parameter.Constant(np.cos((np.pi)/2-c.dec.rad))('cos_gwtheta')
    gwphi = parameter.Constant(c.ra.rad)('gwphi')
    
    
    
    
    log10_fgw = parameter.Constant(log_f)('log10_fgw')
    
    
    if freq>=191.3e-9:
        m = (1./(6**(3./2)*np.pi*freq*u.Hz))*(1./4)**(3./5)*(c.c**3/c.G)
        m_max = np.log10(m.to(u.Msun).value)
    else:
        m_max = 10
    
    log10_mc = parameter.LinearExp(7,m_max)('log10_mc')
    
    phase0 = parameter.Uniform(0, 2*np.pi)('phase0')
    psi = parameter.Uniform(0, np.pi)('psi')
    cos_inc = parameter.Uniform(-1, 1)('cos_inc')
    
    ##sarah's change
    p_phase = parameter.Uniform(0, 2*np.pi)
    #p_dist = parameter.Normal(0, 1)
    p_dist_PX = dists.Dist_PX_Parameter()
    p_dist_DM = dists.Dist_DM_Parameter()
    
    # log10_h = parameter.Uniform(-18, -11)('log10_h')
    log10_dL = parameter.Constant(log_dist)('log10_dL')
    tref = max(tmax)
    
    
    wf_DM = nd.cw_delay_new(cos_gwtheta=cos_gwtheta, gwphi=gwphi, cos_inc=cos_inc,
                      log10_mc=log10_mc, log10_fgw=log10_fgw,
                      log10_dist=log10_dL,
                      phase0=phase0, psi=psi,
                      psrTerm=True, p_dist=p_dist_DM, p_phase=p_phase,
                      evolve=True, check=False,
                      tref=tref)
    cw_DM = CWSignal(wf_DM, ecc=False, psrTerm=True, name = 'cw')
    
    
    wf_PX = nd.cw_delay_new(cos_gwtheta=cos_gwtheta, gwphi=gwphi, cos_inc=cos_inc,
                      log10_mc=log10_mc, log10_fgw=log10_fgw,
                      log10_dist=log10_dL,
                      phase0=phase0, psi=psi,
                      psrTerm=True, p_dist=p_dist_PX, p_phase=p_phase,
                      evolve=True, check=False,
                      tref=tref)
    cw_PX = CWSignal(wf_PX, ecc=False, psrTerm=True, name = 'cw')
    
    
    # log10_Agw = parameter.Constant(-16.27)('gwb_log10_A')
    # gamma_gw = parameter.Constant(6.6)('gwb_gamma')
    log10_Agw =  parameter.Constant(-15.80)('gwb_log10_A')
    gamma_gw = parameter.Constant(6.08)('gwb_gamma')
    cpl = utils.powerlaw(log10_A=log10_Agw, gamma=gamma_gw)
    crn = gp_signals.FourierBasisGP(cpl, components=5, Tspan=Tspan,
                                            name='gw')
    
    tm = gp_signals.TimingModel()
    
    
    
    eph = deterministic_signals.PhysicalEphemerisSignal(use_epoch_toas=True, model = 'setIII')
    
    
    
    s_DM = ef + eq + ec + rn + cw_DM + eph + tm #+ crn
    
    s_PX = ef + eq + ec + rn + cw_PX + eph + tm #+ crn
    
    file = 'pulsar_distances_12p5.json'
    with open(file, 'r') as fp:
        dists_file = json.load(fp)
    
    psrs_DM = [psr for psr in psrs if 'DM' in dists_file[psr.name]]
    psrs_PX = [psr for psr in psrs if 'PX' in dists_file[psr.name]]
    
    
    #models = [s(psr) for psr in psrs]
    models = [s_DM(psr) for psr in psrs_DM]
    models.extend([s_PX(psr) for psr in psrs_PX])
    
    pta = signal_base.PTA(models)
    pta.set_default_params(noisedict)
    
    
    xs = {par.name: par.sample() for par in pta.params}
    print(pta.get_lnlikelihood(xs));
    print(pta.get_lnprior(xs));
    #print(pta.param_names)
    
    
    
    
    
    
    ndim= len(pta.param_names)
    
    # initial jump covariance matrix
    cov = np.diag(np.ones(ndim) * 0.1**2)
    
    # parameter groupings
    #groups = fns.get_parameter_groups(pta)
    
    groups = dists.get_parameter_groups_CAW_target(pta)
    
    
    
    sampler = ptmcmc(ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, groups=groups, outDir=outdir, resume=True)
    
    
    with open(outdir+'/parameters.json', 'w') as fp:
        json.dump(pta.param_names, fp)
    
    ##from sarah
    inc_psr_term = True
    ##why this??
    if inc_psr_term:
        psr_dist = {}
        for psr in psrs:
            psr_dist[psr.name] = psr.pdist
    else:
        psr_dist = None
    
    #init = None
    
    
    
    
    
        
    jp = JumpProposal(pta, empirical_distr=empirical_distr)
    jpCW = dists.JumpProposalCW(pta, fgw=freq, psr_dist=psr_dist,empirical_distr=empirical_distr)
        
    if 'red noise' in jp.snames:
        sampler.addProposalToCycle(jp.draw_from_red_prior, 30)
    if empirical_distr:
        sampler.addProposalToCycle(jp.draw_from_empirical_distr, 30)
    
    for group in groups:
        print(np.array(pta.param_names)[np.array(group)])
    
    #sampler.addProposalToCycle(jp.draw_from_crn_prior, 20)
    cw = True
    if cw:
        sampler.addProposalToCycle(jp.draw_from_cw_prior, 20) #pick a cw param & jump
        #sampler.addProposalToCycle(jp.draw_from_cw_log_uniform_distribution, 10) #draw from uniform h
        
        sampler.addProposalToCycle(jp.draw_from_par_prior(['log10_mc']), 5) #draw from uniform Mc
        #sampler.addProposalToCycle(jp.draw_from_par_prior(['cos_gwtheta', 'gwphi']), 10)#draw skypos prior
        
        
        
        #next few have been added to CW's e_e
        # sampler.addProposalToCycle(jpCW.draw_strain_skewstep, 10)
        # sampler.addProposalToCycle(jpCW.draw_gwphi_comb, 5)
        # sampler.addProposalToCycle(jpCW.draw_gwtheta_comb, 5)
    ##these aren't needed for constant sky position
    if 'phase0' in pta.param_names and 'psi' in pta.param_names:
        sampler.addProposalToCycle(jpCW.phase_psi_reverse_jump, 1)
    # if 'psi' in pta.param_names:
    #     sampler.addProposalToCycle(jpCW.draw_strain_psi, 2)
    # if 'cos_inc' in pta.param_names:
    #     sampler.addProposalToCycle(jpCW.draw_strain_inc, 2)
        
    if inc_psr_term:
        pdist_pars = [p for p in pta.param_names if 'p_dist' in p]
        pphase_pars = [p for p in pta.param_names if 'p_phase' in p]
        sampler.addProposalToCycle(jpCW.draw_from_many_par_prior(pdist_pars, 'p_dist'), 30) #draw from p_dists
        sampler.addProposalToCycle(jpCW.draw_from_many_par_prior(pphase_pars, 'p_phase'), 30) #draw from p_phases
    eph = True
    if eph:
            
        ephempars = []
        
        for sc in pta._signalcollections:
            for signal in sc._signals:
                if signal.signal_name == 'phys_ephem':
                    ephempars.extend(signal.param_names)
        ephempars = np.unique(ephempars)
        
        jup_pars = np.array([p for p in ephempars if 'jup_orb' in p])
        print(ephempars, jup_pars)
        sampler.addProposalToCycle(jp.draw_from_ephem_prior, 20)
        sampler.addProposalToCycle(jp.draw_from_par_prior('jup_orb_elements'), 20)
    
    sampler.addAuxilaryJump(jpCW.fix_cyclic_pars)
    
    sampler.addProposalToCycle(jp.draw_from_prior, 10) #pick a param & jump
    
    
    
    
    x0 = np.hstack(p.sample() for p in pta.params)
    sampler.sample(x0, N, SCAMweight=25, AMweight=40, DEweight=0)
    
    
if os.path.isfile(outdir+'/chain_1.txt'):
    chain = np.loadtxt(outdir+'/chain_1.txt')
    l = np.shape(chain)[0]
    if l< 189900:
        N = int(N_iter-10*l)
        run_ul(N)
    else:
        print('already enough samples')
else:
    N = int(N_iter)
    run_ul(N)

