#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);

from string import *
import numpy as np
from string_lib import *
from DELPHI import Qs_from_RF_param
from particle_param import *



def FCChh_param(E=50e12,Qxfrac=0.31,Qyfrac=0.32,V=16e6):

    ''' generate typical FCChh parameters, from the beam energy E in eV, the fractional parts of the tunes and the RF voltage in V.
    Outputs:
    - machine: string with machine name(here 'TLEP'+option),
    - E: same as input (beam energy in eV),
    - gamma: relativistic mass factor,
    - sigmaz: RMS bunch length in m,
    - taub: total bunch length in s (4*RMS),
    - R: machine pysical radius (circumference/(2 pi)),
    - Qx: total horizontal tune (integer + fractional parts),
    - Qxfrac: fractional horizontal tune,
    - Qy: total vertical tune (integer + fractional parts),
    - Qyfrac: fractional vertical tune,
    - Qs: synchrotron tune,
    - eta: slippage factor (alpha_p-1/gamma^2),
    - f0: revolution frequency,
    - omega0: revolution angular frequency=2pi*f0,
    - omegas: synchrotron angular frequency=Qs*omega0,
    - dphase: phase of damper w.r.t. "normal" purely resistive damper (0),
    - Estr: string with energy (e.g. '50TeV').
    '''

    e,m0,c,E0=proton_param();
    # E is the energy in eV
    Estr=float_to_str(E/1e12)+'TeV'
    machine='FCChh';
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    circ=100200; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    sigmaz=8e-2;
    Qxint=120;Qyint=120; # based on beta~R/Q with beta proportional to E^(1/3) and taking as a reference the LHC case
    Qx=Qxint+Qxfrac;
    Qy=Qyint+Qyfrac;
    alphap=1.0/110**2
    
    h=133650 # approximate harmonic number
    taub=4.*sigmaz/(beta*c); # full length in s
    eta=alphap-1./(gamma*gamma); # slip factor
    Qs=Qs_from_RF_param(V,h,gamma,eta,phis=0.,particle='proton');
    omegas=Qs*omega0;
    
    dphase=0.; # additional damper phase
    
    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr;




def FCChh_param_2016_09(E=50e12,Qxfrac=0.31,Qyfrac=0.32,V=16e6):

    ''' generate typical FCChh parameters, from the beam energy E in eV, the fractional parts of the tunes and the RF voltage in V.
    Outputs:
    - machine: string with machine name(here 'TLEP'+option),
    - E: same as input (beam energy in eV),
    - gamma: relativistic mass factor,
    - sigmaz: RMS bunch length in m,
    - taub: total bunch length in s (4*RMS),
    - R: machine pysical radius (circumference/(2 pi)),
    - Qx: total horizontal tune (integer + fractional parts),
    - Qxfrac: fractional horizontal tune,
    - Qy: total vertical tune (integer + fractional parts),
    - Qyfrac: fractional vertical tune,
    - Qs: synchrotron tune,
    - eta: slippage factor (alpha_p-1/gamma^2),
    - f0: revolution frequency,
    - omega0: revolution angular frequency=2pi*f0,
    - omegas: synchrotron angular frequency=Qs*omega0,
    - dphase: phase of damper w.r.t. "normal" purely resistive damper (0),
    - Estr: string with energy (e.g. '50TeV').
    '''

    e,m0,c,E0=proton_param();
    # E is the energy in eV
    Estr=float_to_str(E/1e12)+'TeV'
    machine='FCChh';
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    circ=101898.2192; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    sigmaz=8e-2;
    Qxint=120;Qyint=120; # based on beta~R/Q with beta proportional to E^(1/3) and taking as a reference the LHC case
    Qx=Qxint+Qxfrac;
    Qy=Qyint+Qyfrac;
    alphap=8.9E-5
    
    h=133650 # approximate harmonic number
    taub=4.*sigmaz/(beta*c); # full length in s
    eta=alphap-1./(gamma*gamma); # slip factor
    Qs=Qs_from_RF_param(V,h,gamma,eta,phis=0.,particle='proton');
    omegas=Qs*omega0;
    
    dphase=0.; # additional damper phase
    
    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr;

