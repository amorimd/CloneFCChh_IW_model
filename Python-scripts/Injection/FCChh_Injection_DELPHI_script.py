#!/usr/bin/python
# use "./FCChh_[...].py restore launch Results_folder_name" to restore impedance model results and launch DELPHI scans
# use "./FCChh_[...].py restore retrieve Results_folder_name" to retrieve the results of DELPHI scans
# the results are stored in Results_folder_name in the scenario folder
# if no folder is given, use a default folder

import sys
import commands

if len(sys.argv)>3: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);ResultsFolderDELPHI=str(sys.argv[3]);
elif len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);ResultsFolderDELPHI='DELPHIrun';
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;ResultsFolderDELPHI=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;ResultsFolderDELPHI=None;
print lxplusbatchImp,lxplusbatchDEL,ResultsFolderDELPHI;

from string import *
import time
import numpy as np
from copy import deepcopy
#import pylab
import os,re
path_here=os.getcwd()+"/";
from io_lib import *
from tables_lib import select_in_table
from particle_param import *
from Impedance import *
from DELPHI import *
from FCChh_param import *
#from BSH_coll import FCChh_param
import pickle as pkl
import inspect
import datetime

#def FCChh_param(E=50e12,Qxfrac=0.31,Qyfrac=0.32,V=16e6):
#
#    ''' generate typical FCChh parameters, from the beam energy E in eV, the fractional parts of the tunes and the RF voltage in V.
#    Outputs:
#    - machine: string with machine name(here 'TLEP'+option),
#    - E: same as input (beam energy in eV),
#    - gamma: relativistic mass factor,
#    - sigmaz: RMS bunch length in m,
#    - taub: total bunch length in s (4*RMS),
#    - R: machine pysical radius (circumference/(2 pi)),
#    - Qx: total horizontal tune (integer + fractional parts),
#    - Qxfrac: fractional horizontal tune,
#    - Qy: total vertical tune (integer + fractional parts),
#    - Qyfrac: fractional vertical tune,
#    - Qs: synchrotron tune,
#    - eta: slippage factor (alpha_p-1/gamma^2),
#    - f0: revolution frequency,
#    - omega0: revolution angular frequency=2pi*f0,
#    - omegas: synchrotron angular frequency=Qs*omega0,
#    - dphase: phase of damper w.r.t. "normal" purely resistive damper (0),
#    - Estr: string with energy (e.g. '50TeV').
#    '''
#
#    e,m0,c,E0=proton_param();
#    # E is the energy in eV
#    Estr=float_to_str(E/1e12)+'TeV';print(Estr)
#    machine='FCChh';
#    gamma=E*e/E0
#    beta=np.sqrt(1.-1./(gamma**2))
#    circ=101898.2192; # total circumference in m
#    R=circ/(2.*np.pi) # machine radius
#    f0=c*beta/circ # rev. frequency
#    omega0=2.*np.pi*f0;
#    sigmaz=8e-2;
#    Qxint=120;Qyint=120; # based on beta~R/Q with beta proportional to E^(1/3) and taking as a reference the LHC case
#    Qx=Qxint+Qxfrac;
#    Qy=Qyint+Qyfrac;
#    alphap=8.9E-5
#    
#    h=133650 # approximate harmonic number
#    taub=4.*sigmaz/(beta*c); # full length in s
#    eta=alphap-1./(gamma*gamma); # slip factor
#    Qs=Qs_from_RF_param(V,h,gamma,eta,phis=0.,particle='proton');
#    omegas=Qs*omega0;
#    
#    dphase=0.; # additional damper phase
#    
#    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr;



if __name__ == "__main__":

    # beam parameters
    e,m0,c,E0=proton_param();
    
    # machine parameters
    machine2save='FCChh'; 

    # subdirectory (inside DELPHI_results/[machine2save]) where to put the results
    RunDir='Injection_3p3TeV/';
    ResultDir='/afs/cern.ch/work/d/damorim/work/DELPHI_results/'+machine2save+'/'+RunDir;
    os.system("mkdir -p "+ResultDir);

    # flags for plotting and DELPHI
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    strnorm=[''];
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    flagplot=True; # to write impedance files by elements
    nevery=1; # downsampling of the impedance (take less points than in the full model)
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)


    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=200; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y','b--','r--','g--','m--','k--','c--','y--']; # colors
    linetype=['-','--',':'];

    # scan definition
    scenarioscan = np.array(['BS_coll_3p3TeV','BS_only_3p3TeV']) 
    print scenarioscan
    

    model='FCChh'
    Escan=np.array([3.3e12 for ii in scenarioscan])
    E=3.3e12
    V=12e6
    subscan=np.arange(0,len(Escan))
    print subscan


    #DELPHI convergence criterion
    crit=5e-2
    abseps=1e-4

    # setting the scans
    planes=['x'];
    Qpscan=np.arange(0,2,1);
    dampscan=np.array([0.0]); # damper gain scan
    Nbscan=np.array([1.0e11])
    Mscan=np.array([1]); # scan on number of bunches
    imp_fact=1. #impedance factor

    # Longitudinal distribution
    typelong='Gaussian'

    queue='2nd'

    simulation_parameters={\
            'Simulation_time':datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),\
            'DELPHI_version':inspect.getmodule(DELPHI_wrapper).__file__,\
            'flagdamperimp':flagdamperimp,\
            'strnorm':strnorm,\
            'flagnorm':flagnorm,\
            'flagplot':flagplot,\
            'nevery':nevery,\
            'wake_calc':wake_calc,\
            'kmax':kmax,\
            'kmaxplot':kmaxplot,\
            'crit':crit,\
            'abseps':abseps,\
            'scenarioscan':scenarioscan,\
            'model':model,\
            'Escan':Escan,\
            'subscan':subscan,\
            'planes':planes,\
            'Qpscan':Qpscan,\
            'dampscan':dampscan,\
            'Nbscan':Nbscan,\
            'Mscan':Mscan,\
            'imp_fact':imp_fact,\
            'queue':queue,\
            }

    # initialize impedance model and tune shifts
    tuneshiftQp=np.zeros((len(subscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    
    tuneshiftm0Qp=np.zeros((len(subscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios






    for iscenario,scenario in enumerate(scenarioscan[subscan]):

        root_result=ResultDir+scenarioscan[subscan[iscenario]]+'/'
        os.system("mkdir -p "+root_result);

	if machine2save=='FCChh':
		machine_str,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=FCChh_param(E=E,V=V)
	
	g,a,b=longdistribution_decomp(taub,typelong=typelong);
        avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
        simulation_parameters.update({'g':g,'a':a,'b':b,'typelong':typelong})
        print "scenario: ",scenario

        # compute imp. model
        if ((lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp.startswith('launch')))and(lxplusbatchDEL==None):

            print 'Error'
            imp_mod=[]; wake_mod=[];

            if machine2save=='FCChh':
                imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,dire=path_here+"../LHC_elements/",commentcoll=comment_coll_machine,direcoll=dircollscan[subscan[iscenario]]+'/',lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario);


        elif (lxplusbatchImp.startswith('restore')):
            print 'Loading from impedance database...'+scenario
            imp_mod=[]; wake_mod=[];
	    suffix='_'+machine2save+'_Allthemachine_'+Estr+'_'+scenario+'.dat';
            freq_mod,Z_mod=readZ(root_result+"Zxdip"+suffix);
	    Z_mod*=imp_fact

	    imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=freq_mod,func=Z_mod));

            freq_mod,Z_mod=readZ(root_result+"Zydip"+suffix);
	    Z_mod*=imp_fact
		

            imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=freq_mod,func=Z_mod));

        imp_mod_list.append(imp_mod);
        wake_mod_list.append(wake_mod);


        if (lxplusbatchImp.startswith('retrieve'))and(lxplusbatchDEL==None):

            # write Ascii files with each component
            write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+'_'+scenario,
                listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
                dire=root_result+'/')

            if (wake_calc):
            # write Ascii files with each component
                write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+'_'+scenario,
                    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
                dire=root_result+'/')
                # wake for HEADTAIL
                write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+'_'+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)
                # dip only
                write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+'_'+scenario+'_dip.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=2)

	
    if (lxplusbatchDEL!=None):
	    if (lxplusbatchDEL.startswith('launch'))or(lxplusbatchDEL.startswith('retrieve')):

		# DELPHI scans now
		for iscenario,scenario in enumerate(scenarioscan[subscan]):

		    root_result=ResultDir+scenario;
		    print 'DELPHI computation for '+scenario
		    if model=='FCChh':

		    	machine_str,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=FCChh_param(E=E,V=V);
			
                        
                        if (lxplusbatchDEL.startswith('launch')):
                            #ResultsFolderDELPHI_idx=last_folder_index+1
                            #pkl.dump(ResultsFolderDELPHI_idx, open(root_result+'/last_folder_index.pkl','wb'))
                                                            
                            with open(root_result+'/parameters_'+ResultsFolderDELPHI+'.pkl', 'w') as f:
			        #pkl.dump(machine, f)
                                pkl.dump(simulation_parameters,f)
			    f.close()


                        os.system("mkdir -p "+root_result+'/'+ResultsFolderDELPHI);
                        


		    # DELPHI run
		    tuneshiftQp[iscenario,:,:,:,:,:,:,:,:],tuneshiftm0Qp[iscenario,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod_list[iscenario],Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],omega0,Qx,Qy,gamma,eta,a,b,taub,g,planes,nevery=nevery,particle='proton',flagnorm=0,flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=True,lxplusbatch=lxplusbatchDEL,comment=machine_str+scenario+'_'+float_to_str(round(E/1e9))+'GeV_Z'+str(imp_fact),queue=queue,dire=root_result+'/'+ResultsFolderDELPHI+'/',flagQpscan_outside=True);


	    # now the most unstable modes
	    if (lxplusbatchDEL.startswith('retrieve')):
                for iscenario,scenario in enumerate(scenarioscan[subscan]):

                    #last_folder_index=pkl.load(open(root_result+'/last_folder_index.pkl','rb'))             
                    #ResultsFolderDELPHI='ResultsFolderDELPHI_'+str(ResultsFolderDELPHI_idx)
                                                                                                       
                    # output files name for data vs Qp
                    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
                    root_result=ResultDir+scenarioscan[subscan[iscenario]]+'/'+ResultsFolderDELPHI+'/';

                    np.save(root_result+'/TuneshiftQp.npy',tuneshiftQp[iscenario,:,:,:,:,:,:,:,:])
                    np.save(root_result+'/Tuneshiftm0Qp.npy',tuneshiftm0Qp[iscenario,:,:,:,:,:,:,:])


		for iplane,plane in enumerate(planes):
		    for iM,M in enumerate(Mscan):
			for idamp,damp in enumerate(dampscan):
			    for Nb in Nbscan:
				strpart=['Re','Im'];
				for ir,r in enumerate(['real','imag']):
				    for iscenario,scenario in enumerate(scenarioscan[subscan]):
                                        
                                        #last_folder_index=pkl.load(open(root_result+'/last_folder_index.pkl','rb'))
                                        #ResultsFolderDELPHI='ResultsFolderDELPHI_'+str(ResultsFolderDELPHI_idx)

					# output files name for data vs Qp
					Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
					root_result=ResultDir+scenarioscan[subscan[iscenario]]+'/'+ResultsFolderDELPHI+'/';
					fileoutdataQp=root_result+'/data_vs_Qp_'+machine_str+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane+'_Z'+str(imp_fact);
					fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine_str+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane+'_Z'+str(imp_fact);
					fileoutdata_all=root_result+'/data_vs_Qp_all_'+machine_str+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane+'_Z'+str(imp_fact);

					ts=getattr(tuneshiftQp[iscenario,iplane,iM,:,idamp,np.where(Nbscan==Nb),0,0,0],r);
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")

					tsm0=getattr(tuneshiftm0Qp[iscenario,iplane,iM,:,idamp,np.where(Nbscan==Nb),0,0],r);
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshiftm0")
					
					all_unstable_modes=getattr(tuneshiftQp[iscenario,iplane,iM,:,idamp,np.where(Nbscan==Nb),0,0,:],r);
					data=np.hstack((Qpscan.reshape((-1,1)),all_unstable_modes.reshape((-1,kmaxplot))));
					write_ncol_file(fileoutdata_all+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")


