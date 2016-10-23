#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from string import *
import numpy as np
import pickle as pick
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from string_lib import *
from particle_param import *
from Impedance import *

def read_coll_files(param_filename,settings_filename,beta_filename,namesref=None):
    ''' read collimator files and output parameters.'''
    
    # file with materials, angles and lengths
    if (namesref==None): namesref=read_ncol_file_identify_header(param_filename,'[nN]ame');
    names=read_ncol_file_identify_header(param_filename,'[nN]ame');
    # reorder such that the coll. names match with namesref
    ind=find_ind_names(namesref,names);
    angle=read_ncol_file_identify_header(param_filename,'[aA]ngle');
    length=read_ncol_file_identify_header(param_filename,'[lL]ength');
    material=read_ncol_file_identify_header(param_filename,'[mM]aterial');
    material=[material[i] for i in ind];angle=angle[ind];length=length[ind];
    
    # file with settings
    names=read_ncol_file_identify_header(settings_filename,'[nN]ame');
    halfgap=read_ncol_file_identify_header(settings_filename,'[hH]alfgap');
    # reorder such that the coll. names match with namesref
    ind=find_ind_names(namesref,names);
    halfgap=halfgap[ind];

    # file with beta functions
    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
    # reorder such that the coll. names match with namesref
    ind=find_ind_names(namesref,names);
    betax=betax[ind];betay=betay[ind];
    
    return namesref,material,angle,length,halfgap,betax,betay;

def read_coll_file(filename,namesref=None):
    ''' read collimator file and output parameters.'''
    
    namesref=read_ncol_file_identify_header(filename,'[nN]ame');
    names=read_ncol_file_identify_header(filename,'[nN]ame');
    angle=read_ncol_file_identify_header(filename,'[aA]ngle');
    length=read_ncol_file_identify_header(filename,'[lL]ength');
    material=read_ncol_file_identify_header(filename,'[mM]aterial');
    halfgap=read_ncol_file_identify_header(filename,'[hH]alfgap');
    betax=read_ncol_file_identify_header(filename,'[bB]etax');
    betay=read_ncol_file_identify_header(filename,'[bB]etay');
    
    return namesref,material,angle,length,halfgap,betax,betay;

def LHC_singlecoll_iw_model(name,material,halfgap,angle,gamma,length,
	wake_calc=False,coatingmat=None,coatingthickness=0,fpar=freq_param(),zpar=z_param(),
	lxplusbatch=None,comment='',dire='',queue=None):
    
    ''' construct impedance/wake model (flat chamber) for an LHC collimator 
    with a skew angle (as defined in N. Mounet PhD, p. 56)
    name is the name of the collimator, material its main material,
    angle in rad, halfgap in m
    wake_calc: flag for wake calculation
    if coatingmat is not None, first layer is a coating defined
    by this and coatingthickness
    last layer is always stainless steel 304L
    
    fpar and zpar select the frequencies and distances (for the wake) scans 
    (objects of the classes freq_param and z_param).
    
    lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus on queue 'queue'
    		   if 'retrieve' -> retrieve outputs
    comment is added to the name for IW2D output files
    dire contains the directory name where to put the outputs (default='./'=directory of IW2D)
    if queue is not None, it is the lxbatch queue where to launch the calculation.'''
    
    if (material.startswith('CU'))or(material.startswith('Cu')): material='Cu300K';
    elif (material=='C'): material='graphite';
    elif (material.startswith('HBN')): material='hBN';

    if (queue==None):
	# find on which queue to launch calculation (case when lxplusbatch=='launch')
        queues=['8nm','1nh','1nh','8nh','1nd','2nd','1nw'];
	if wake_calc:
	    if material.startswith('hBN'): iq=5;
	    elif (material=="graphite")or(material=="CFC"): iq=4;
	    else: iq=3;
	    if (gamma>=3000): iq+=1; # flat top settings -> convergence is slower
	else:      
	    nfreq=(np.mod(fpar.ftypescan,2)==0)*((np.log10(fpar.fmax)-np.log10(fpar.fmin))*fpar.nflog+1) + (fpar.ftypescan==1)*round((fpar.fmax-fpar.fmin)/fpar.fsamplin)+(fpar.ftypescan==2)*fpar.nrefine;
	    iq=max(int(np.ceil(np.log10(nfreq/500.))),0); # =0 up to 500 frequencies; 1 up to 5000; 2 up to 50000; etc.
	    iq+=1;
	    if material.startswith('hBN'): iq+=2;
	queue=queues[iq];
    print material, queue, coatingmat, comment;
    
    
    
    # some hard coded parameters
    thickness=25e-3;freqlin=1.e11;
    b=[halfgap];
    
    # construct input file
    layers=[];
    eval('layers.append('+material+'_layer(thickness=thickness))');
    layers.append(ss304L_layer());

    iw_input=impedance_wake_input(gamma=gamma,length=length,b=b,
    	layers=layers,fpar=fpar,zpar=zpar,freqlinbisect=freqlin,geometry='flat',comment=name+comment);
    print(name+comment)
    imp_mod,wake_mod=imp_model_from_IW2D(iw_input,wake_calc=wake_calc,flagrm=True,lxplusbatch=lxplusbatch,queue=queue,dire=dire);
    
    imp_mod_new=rotate_imp_wake(imp_mod,np.pi/2.-angle);
    wake_mod_new=rotate_imp_wake(wake_mod,np.pi/2.-angle);
    
    return imp_mod_new,wake_mod_new;
    

def FCChh_manycoll_iw_model(E,avbetax,avbetay,fileName,wake_calc=False,ftypescan=2,nflog=100,namesref=None,
	coatingmat=None,coatingthickness=0,lxplusbatch=None,comment='',dire='',gapScaling=1.0):

    ''' creates an impedance or wake model for all collimators.
    E is the energy in eV, avbetax and avbetay the average beta functions
    used for the weighting, param_filename is the file with all parameters
    except half-gaps and betas, settings_filename is the file with half-gaps (in m),
    beta_filename the file with beta functions (in m).
    wake_calc selects the wake computation if True, nflog are the number of frequencies
    per decade, and ftypescan the type of frequency scan (0,1 or 2: logarithmic only, 
    linear only or logarithmic with refinement around high-frequency resonance(s) ).
    namesref are the coll. names (from param_filename) to select (if None take all),
    coatingmat and coatingthickness is the info about an added coating.
    lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus on queue 'queue'
    		   if 'retrieve' -> retrieve outputs
    dire contains the directory name where to put the outputs (default='./'=directory of IW2D)
    The gaps are scaled with gapScaling'''

    e,m0,c,E0=proton_param();
    gamma=E*e/E0;
    
    # read files
    namesref,material,angle,length,halfgap,betax,betay=read_coll_file(fileName,namesref=namesref);
    for i,name in enumerate(namesref):
        halfgap[i] *= gapScaling 
        print name,material[i],angle[i],length[i],halfgap[i],betax[i],betay[i]
        

    # main loop to construct model
    imp_mod=[];wake_mod=[];
    for iname,name in enumerate(namesref):
        fminrefine=1.e11;fmaxrefine=5.e12;nrefine=5000;
        print("WakeCalc",wake_calc);
        imp,wake=LHC_singlecoll_iw_model(name,material[iname],
	        halfgap[iname],angle[iname],gamma,length[iname],wake_calc=wake_calc,
        coatingmat=coatingmat,coatingthickness=coatingthickness,fpar=freq_param(ftypescan=ftypescan,
        nflog=nflog,fminrefine=fminrefine,fmaxrefine=fmaxrefine,nrefine=nrefine),lxplusbatch=lxplusbatch,
        comment=comment+'_'+material[iname]+'_'+float_to_str(round(halfgap[iname]*1e5)/1e2)+'mm',dire=dire);

        add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
        add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);
	    
    return imp_mod,wake_mod;




def FCChh_manycoll_iw_model_refine(E,avbetax,avbetay,fileName,wake_calc=False,ftypescan=2,nflog=100,namesref=None,coatingmat=None,coatingthickness=0,lxplusbatch=None,comment='',dire='',gapScaling=1.0):

    ''' creates an impedance or wake model for all collimators.
    E is the energy in eV, avbetax and avbetay the average beta functions
    used for the weighting, param_filename is the file with all parameters
    except half-gaps and betas, settings_filename is the file with half-gaps (in m),
    beta_filename the file with beta functions (in m).
    wake_calc selects the wake computation if True, nflog are the number of frequencies
    per decade, and ftypescan the type of frequency scan (0,1 or 2: logarithmic only, 
    linear only or logarithmic with refinement around high-frequency resonance(s) ).
    namesref are the coll. names (from param_filename) to select (if None take all),
    coatingmat and coatingthickness is the info about an added coating.
    lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus on queue 'queue'
    		   if 'retrieve' -> retrieve outputs
    dire contains the directory name where to put the outputs (default='./'=directory of IW2D)
    The gaps are scaled with gapScaling'''

    e,m0,c,E0=proton_param();
    gamma=E*e/E0;
    
    # read files
    namesref,material,angle,length,halfgap,betax,betay=read_coll_file(fileName,namesref=namesref);
    for i,name in enumerate(namesref):
        halfgap[i] *= gapScaling 
        print name,material[i],angle[i],length[i],halfgap[i],betax[i],betay[i]
        

    # main loop to construct model
    imp_mod=[];wake_mod=[];
    for iname,name in enumerate(namesref):
        fminrefine=1.e11;fmaxrefine=5.e12;nrefine=5000;
        print("WakeCalc",wake_calc);
        imp,wake=LHC_singlecoll_iw_model(name,material[iname],
	        halfgap[iname],angle[iname],gamma,length[iname],wake_calc=wake_calc,
        coatingmat=coatingmat,coatingthickness=coatingthickness,fpar=freq_param(ftypescan=ftypescan,
        nflog=nflog,fminrefine=fminrefine,fmaxrefine=fmaxrefine,nrefine=nrefine),lxplusbatch=lxplusbatch,
        comment=comment+'_'+material[iname]+'_'+float_to_str(round(halfgap[iname]*1e5)/1e2)+'mm_coarse',dire=dire);

        add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
        add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);


    for iname,name in enumerate(namesref):
        fminrefine=2.9e9;fmaxrefine=3.1e9;nrefine=20000;
        print("WakeCalc",wake_calc);
        imp,wake=LHC_singlecoll_iw_model(name,material[iname],
	        halfgap[iname],angle[iname],gamma,length[iname],wake_calc=wake_calc,
        coatingmat=coatingmat,coatingthickness=coatingthickness,fpar=freq_param(ftypescan=ftypescan,
        nflog=nflog,fminrefine=fminrefine,fmaxrefine=fmaxrefine,nrefine=nrefine),lxplusbatch=lxplusbatch,
        comment=comment+'_'+material[iname]+'_'+float_to_str(round(halfgap[iname]*1e5)/1e2)+'mm_fine',dire=dire);

        add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
        add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);
    
    return imp_mod,wake_mod;
