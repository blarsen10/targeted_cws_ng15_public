#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:03:55 2024

@author: bjornlarsen
"""

import numpy as np
import scipy.integrate as spi
import astropy.units as u

c_light = 299792458 # m/s
G = 6.67430e-11 # Nm^2/kg

def scaled_param_prior(c, param, target_param):
    '''
    Scale a params prior distribution from the one in la_forge
    
    c: la_forge core
    param: str
        Name of param whose prior you are scaling
    target_param: str
        Name of a parameter you are using the scaled prior to solve for. This will determine how to scale the prior
    '''
    # assuming log10_mc prior in units of Msun (for log10_h0 calculation)
    if param == 'log10_mc' and target_param == 'log10_h0':
        pline = [p for p in c.priors if 'log10_mc' in p][0]
        pmin = float(pline[pline.index('pmin')+5:pline.index(', pmax')])
        pmax = float(pline[pline.index('pmax')+5:-1])
        log10_mc_prior = np.array([pmin, pmax])
        log10_mc_prior_scaled = (log10_mc_prior +
                                 np.log10(u.Msun.to(u.kg)*G/c_light**3))
        log10_mc_prior_scaled
        line = [c.runtime_info[i] for i in range(len(c.runtime_info))
                if 'log10_dL' in c.runtime_info[i]][0]
        log10_dL = float(line.replace('log10_dL:Constant=',''))
        log10_dL_scaled = log10_dL + np.log10(u.Mpc.to(u.m)/c_light)
        return 5/3*log10_mc_prior_scaled - log10_dL_scaled + np.log10(2*np.pi**(2/3))
    # assuming log10_fgw prior in units of 1/s (for log10_h0 calculation)
    elif param == 'log10_fgw' and target_param == 'log10_h0':
        pline = [p for p in c.priors if 'log10_fgw' in p][0]
        pmin = float(pline[pline.index('pmin')+5:pline.index(', pmax')])
        pmax = float(pline[pline.index('pmax')+5:-1])
        log10_fgw_prior = np.array([pmin, pmax])
        return 2/3*log10_fgw_prior
    else:
        raise ValueError('Not sure what scaled prior range you want'
                         f'Param you are trying to scale: {param}'
                         f'Target param: {target_param}'
                         'either need to fix param labelling, or further code development required')

def integrand_uniform_uniform(x, z, ax, bx, ay, by):
    if (x > np.max([ax, z - by]))*(x < np.min([ay, z - bx])):
        return 1/(bx - ax)/(by - ay)
    else:
        return 0
    
def integrand_exp_uniform(x, z, ax, bx, ay, by):
    if (x > np.max([ax, z - by]))*(x < np.min([ay, z - bx])):
        return np.log(10)/(by-ay)/(10**bx - 10**ax)*10**x
    else:
        return 0
    
def integrand_norm_uniform(x, z, mu, sigma, a, b):
    if (x > z - b)*(x < z - a):
        return 1/(np.sqrt(2*np.pi**(2/3))*sigma*(b-a))*np.exp(-(x-mu)**2/(2*sigma))
    else:
        return 0

def pdf_Z(z, x_arr, x_prior, y_prior, fX='uniform'):
    #return spi.quad(integrand, x_arr[0], x_arr[-1], args=(z), limit=200)[0]
    if fX == 'uniform':
        ax, bx = x_prior
        ay, by = y_prior
        return spi.simpson([integrand_uniform_uniform(x,z,ax,bx,ay,by) for x in x_arr], x=x_arr)
    elif fX == 'normal':
        mu, sigma = x_prior
        a, b = y_prior
        return spi.simpson([integrand_norm_uniform(x,z,mu,sigma,a,b) for x in x_arr], x=x_arr)
    elif fX == 'exp':
        ax, bx = x_prior
        ay, by = y_prior
        return spi.simpson([integrand_exp_uniform(x,z,ax,bx,ay,by) for x in x_arr], x=x_arr)

def convolved_prior(c, param, detection=True, Npoints=1000):
    '''
    This point of this function is to numerically calculate the prior for a derived parameter
    which comes from the sum of two parameters sampled using a uniform distribution.
    Unlike the above, this version is agnostic to the convolved parameters and instead
    uses information from the prior directly, as long as the prior string as format:
        'param:Convolve(Uniform(pmin=#, pmax=#), Uniform(pmin=#, pmax=#))'
    Currently only convolution of either a uniform, normal, or exponential with a uniform
    distribution is supported.
    
    core: la_forge core
        This should ideally contain chains/priors for all 3 params
    param: str
        Name of param whose prior you are solving for
        
    '''
    idx = c.params.index(param)
    pline = c.priors[idx]
    p1_line = pline[pline.index('Convolve(')+len('Convolve('):pline.index(', Uniform')]
    p2_line = pline[pline.index(', Uniform')+2:-1]
    if 'Normal' in p1_line:
        fX = 'normal'
        mu = float(p1_line[p1_line.index('mu')+3:p1_line.index(', sigma')])
        sigma = float(p1_line[p1_line.index('sigma')+6:-1])
        x_prior = [mu, sigma]
        xmin = mu - 3*sigma
        xmax = mu + 3*sigma
    else:
        xmin = float(p1_line[p1_line.index('pmin')+5:p1_line.index(', pmax')])
        xmax = float(p1_line[p1_line.index('pmax')+5:-1])
        x_prior = [xmin, xmax]
        if 'Uniform' in p1_line:
            fX = 'uniform'
        elif 'LinearExp' in p1_line:
            fX = 'exp'
    ymin = float(p2_line[p2_line.index('pmin')+5:p2_line.index(', pmax')])
    ymax = float(p2_line[p2_line.index('pmax')+5:-1])
    y_prior = [ymin, ymax]
    x_arr = np.linspace(xmin, xmax, Npoints)
    z_arr = np.linspace(xmin + ymin, xmax + ymax, Npoints)
    p_z = np.array([pdf_Z(z, x_arr, x_prior, y_prior, fX=fX) for z in z_arr])
    return z_arr, p_z