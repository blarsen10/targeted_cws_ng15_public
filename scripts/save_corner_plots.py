#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 14:41:16 2024

@author: bjornlarsen
"""

import numpy as np
import matplotlib.pyplot as plt
import json, os, pickle, glob#, natsort
import corner
import astropy.units as u
import scipy.stats as sps
import la_forge.diagnostics as dg
import la_forge.core as co

def get_prior_distr(core, param):
    idx = core.params.index(param)
    if param == 'log10_h0':
        pline = core.priors[idx-4]
    else:
        pline = core.priors[idx]
    prior_type = pline[pline.index(':')+1:pline.index('(')]
    # setup x-axis for the plot
    if prior_type == 'Uniform':
        pmin = float(pline[pline.index('pmin')+5:pline.index(', pmax')])
        pmax = float(pline[pline.index('pmax')+5:-1])
        x = np.linspace(pmin, pmax, 300)
        prior_dist = sps.uniform(loc=pmin, scale=pmax-pmin)
    if prior_type == 'Normal':
        mu = float(pline[pline.index('mu')+3:pline.index(', sigma')])
        sigma = float(pline[pline.index('sigma')+6:-1])
        prior_dist = sps.norm(loc=mu, scale=sigma)
        pmin = np.min([mu-3*sigma, np.min(core.chain[core.burn:,idx])])
        pmax = np.max([mu+3*sigma, np.max(core.chain[core.burn:,idx])])
        x = np.linspace(pmin, pmax, 300)
    y = prior_dist.pdf(x)
    return x, y

def get_mc_prior_mean(t_str):
    prior_path = f'{project_path}/priors/{t_str}_priors.json'
    with open(prior_path, 'r') as f:
        astro_priors = json.load(f)
    return astro_priors['log10_Mc']

def corner_single(t, t_str, c, vary_crn=False, save=True):
    pars_select = ['cos_inc','phase0','psi', 'log10_h0', 'log10_mc']
    #units = ['rad', 'rad', 'rad', 'calculated', r'$M_\odot$']
    titles = [r'$\cos\iota$', '$\Phi_0$', r'$\psi$', r'$\log_{10}h_0$', r'$\log_{10}\mathcal{M}_c$']
    if vary_crn:
        pars_select = ['crn_gamma', 'crn_log10_A'] + pars_select
        #units = ['', ''] + units
        titles = [r'$\gamma_{\rm{CRN}}$', r'$\log_{10}A_{\rm{CRN}}$'] + titles

    # plot
    labels = titles#[p+' \n ('+u+')' if u else p+' \n' for p,u in zip(titles, units)]
    idxs = [c[t].params.index(p) for p in pars_select]
    fig = corner.corner(c[t].chain[c[t].burn:, idxs], labels=labels,
                        #quantiles=quantiles, title_quantiles=quantiles,
                        hist_kwargs={'density':True},
                        titles=titles, show_titles=True, levels=(0.68, 0.95),
                        label_kwargs={'fontsize': 20}, title_kwargs={'fontsize': 14})

    # Extract the axes
    axes = np.array(fig.axes).reshape((len(idxs), len(idxs)))

    # Loop over the diagonal
    for i, p in enumerate(pars_select):
        ax = axes[i, i]
        x, y = get_prior_distr(c[t], p)
        ax.plot(x, y, 'C2')
        ax.set_xlim([x.min(),x.max()])
        if p == 'log10_mc':
            ax.axvline(get_mc_prior_mean(t_str), color='r')
    # Loop over the histograms
    for yi, p1 in enumerate(pars_select): # rows
        for xi, p2 in enumerate(pars_select[:yi]): # cols
            ax = axes[yi, xi]
            y, _ = get_prior_distr(c[t], p1)
            x, _ = get_prior_distr(c[t], p2)
            ax.set_xlim([x.min(),x.max()])
            ax.set_ylim([y.min(),y.max()])
    if save:
        fig.savefig(f'{save_loc}/{t_str}_log_uniform.pdf', format='pdf', dpi=600)
    return fig

project_path = '/vast/palmer/home.grace/bbl29/targeted_cws_ng15'
outdir_path = '/vast/palmer/home.grace/bbl29/project/targeted_cws_ng15/data/chains'
save_loc = f'{project_path}/reports/figures/pub_figs'
dataset = 'ng15_bipm2019_old_Tspans'
c_light = 299792458 # m/s
G = 6.67430e-11 # Nm^2/kg

# =============================================================================
# Log-uniform priors
# =============================================================================

target_paths = glob.glob(f'{outdir_path}/{dataset}/*_UL')
[target_paths.pop(target_paths.index(tp)) for tp in np.flip(target_paths) if 'altskyloc' in tp]
[target_paths.pop(target_paths.index(tp)) for tp in np.flip(target_paths) if 'altskyloc' in tp]
target_strs = [tp.replace(f'{outdir_path}/{dataset}/','').replace('_UL','') for tp in target_paths]
targets = np.copy(target_strs)
Nt = len(targets)
for i in range(Nt):
    targets[i] = targets[i].replace('_',' ')
    targets[i] = targets[i].replace('-','$-$')
    if '16nHz' in targets[i]:
        targets[i] = targets[i].replace(' 16nHz','')
targets = sorted(targets)
target_strs = sorted(target_strs)
print(Nt,'targets')
print(targets)

c = {}
for t, tp in zip(targets[len(c.keys()):], target_paths[len(c.keys()):]):
    print(t)
    c[t] = co.Core(corepath=f'{tp}/core.h5', label=t)
    
for t, tp in zip(targets, target_paths):
    prior_path = glob.glob(tp + '/*/priors.txt')[0]
    c[t].priors = np.loadtxt(prior_path, dtype=str, delimiter='\t')
    info_path = glob.glob(tp + '/*/runtime_info.txt')[0]
    c[t].runtime_info = np.loadtxt(info_path, dtype=str, delimiter='\t')
    
# append param name
for t in targets:
    if 'log10_h0' in c[t].params:
        print('skipping', t, 'params')
    else:
        print('adding', t, 'params')
        c[t].params.append('log10_h0')

# append chain and prior
for t in targets:
    if len(c[t].params) == c[t].chain.shape[1]:
        print('skipping', t, 'chain')
    else:
        print('adding', t, 'chain')
        # log10_fgw
        line = [c[t].runtime_info[i] for i in range(len(c[t].runtime_info))
                if 'log10_fgw' in c[t].runtime_info[i]][0]
        log10_fgw = float(line.replace('log10_fgw:Constant=',''))
        # log10_dL
        line = [c[t].runtime_info[i] for i in range(len(c[t].runtime_info))
                if 'log10_dL' in c[t].runtime_info[i]][0]
        log10_dL = float(line.replace('log10_dL:Constant=',''))
        log10_dL_scaled = log10_dL + np.log10(u.Mpc.to(u.m)/c_light)
        # log10_mc chain
        log10_mc = c[t]('log10_mc',to_burn=False)
        log10_mc_scaled = log10_mc + np.log10(u.Msun.to(u.kg)*G/c_light**3)
        # calculate strain
        log10_h0 = 5*log10_mc_scaled/3 + 2*log10_fgw/3 - log10_dL_scaled + np.log10(2*np.pi**(2/3))
        c[t].chain = np.vstack([c[t].chain.T,log10_h0]).T

# append prior
for t in targets:
    if len(c[t].priors) == len(c[t].params) - 4:
        print('skipping', t, 'priors')
    else:
        print('adding', t, 'priors')
        # log10_fgw
        line = [c[t].runtime_info[i] for i in range(len(c[t].runtime_info))
                if 'log10_fgw' in c[t].runtime_info[i]][0]
        log10_fgw = float(line.replace('log10_fgw:Constant=',''))
        # log10_dL
        line = [c[t].runtime_info[i] for i in range(len(c[t].runtime_info))
                if 'log10_dL' in c[t].runtime_info[i]][0]
        log10_dL = float(line.replace('log10_dL:Constant=',''))
        log10_dL_scaled = log10_dL + np.log10(u.Mpc.to(u.m)/c_light)
        # get log10_mc prior
        pline = [p for p in c[t].priors if 'log10_mc' in p][0]
        pmin = float(pline[pline.index('pmin')+5:pline.index(', pmax')])
        pmax = float(pline[pline.index('pmax')+5:-1])
        log10_mc_prior = np.array([pmin, pmax])
        log10_mc_prior_scaled = log10_mc_prior + np.log10(u.Msun.to(u.kg)*G/c_light**3)
        # calculate strain prior
        log10_h0_prior = 5*log10_mc_prior_scaled/3 + 2*log10_fgw/3 - log10_dL_scaled + np.log10(2*np.pi**(2/3))
        log10_h0_prior_str = f'log10_h0:Uniform(pmin={log10_h0_prior[0]}, pmax={log10_h0_prior[1]})'
        c[t].priors = np.concatenate([c[t].priors, [log10_h0_prior_str]])
    
for t, t_str in zip(targets, target_strs):
    print(t)
    corner_single(t, t_str, c, save=True)
    plt.close('all')
    
# =============================================================================
# Normal priors
# =============================================================================

# Use the already compiled target names for this
target_paths = [f'{outdir_path}/{dataset}/{t}' for t in target_strs]
for i in range(Nt):
    if 'NGC_3115' in target_paths[i]:
        target_paths[i] = target_paths[i] + '_16nHz'

c = {}
for t, tp in zip(targets[len(c.keys()):], target_paths[len(c.keys()):]):
    print(t)
    c[t] = co.Core(corepath=f'{tp}/core.h5', label=t)

