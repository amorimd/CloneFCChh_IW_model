#!/usr/bin/python

import sys
import pickle as pick

if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;
print(lxplusbatchImp,lxplusbatchDEL);   

from string import *
import numpy as np
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure,cmap
from particle_param import *
from Impedance import *
from DELPHI import *
from FCChh_coll_imp import FCChh_manycoll_iw_model
from FCChh_param import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    
    #mu0=4e-7*np.pi;Z0=mu0*c;
    E=3.3e12;V=12e6;
    #E=50e12;V=16e6;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=FCChh_param(E=E,V=V);
    circ=2*np.pi*R;
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used in DELPHI

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    pat=['-','--'];

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=1; # 0 to avoid computing (simply plot from existing data)
    flagSach=(lxplusbatchDEL=='launch'); # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    #flagSach=0; # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    wake_calc=False;
    
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=50; # number of plotted eigenvalues
    
    
    scenario_name='BS_only_'+Estr
    root_result='/afs/cern.ch/work/d/damorim/work/DELPHI_results/'+machine+'/'+scenario_name;
    os.system("mkdir -p "+root_result);
    suffix=''; # suffix for output files 
    lmaxSach=1;
    
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # scan definitions
    materialscan=['Cu50K']; # material is cold Cu (5K should be ~same as 20K)
    scenarioscan=deepcopy(materialscan);
    hgapscan=1e-3*np.array([13]);
    Nbscan=np.arange(0e10,25.2e11,2.0e10);
    Nbscanplot=np.array([1.e9,1e10,5e10,1.e11,5.e11,1.e12]);
    dampscan=np.array([0,0.02,0.01]);
    Mscan=np.array([1,10600]); # number of equidistant bunches
    #Qpscan=np.arange(-10,21);
    #Qpplothg=np.array([0,2,5,10]); # scan in Qp for plot vs half-gap (and for TMCI plot)
    Qpscan=np.arange(-20,21,1);
    Qpplothg=np.array([0]); # scan in Qp for plot vs half-gap (and for TMCI plot)

    imp_mod_list={}; # complete list of impedance scenarios
    wake_mod_list={};# complete list of wake scenarios
    
    for imat,material in enumerate(materialscan):
    
        name='vacpipe_'+Estr+'_'+material;
        imp_mod_list[material]=[];leg=[];

        for ihalfgap,halfgap in enumerate(hgapscan):

            imp_mod=[];wake_mod=[];
    
            # compute model for vacuum pipe
            strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';
            
            # B field found from energy and radius of curvature (arc filling factor=0.78, 1400*12 m of straigth sections)
            #radius=(circ-1440*12)*0.78/2/np.pi;
            layers_iw=[eval(material+'_layer(thickness=np.inf,RRR=70,B=16.0*E/50e12)')];
            print(material,", resistivity=",layers_iw[0].rhoDC,"Ohm.m,tau=",layers_iw[0].tau,'s');

            betax=avbetax;betay=avbetay;

            fpar=freq_param(fmin=1,fmax=1e14,ftypescan=2,nflog=40,fminrefine=1e11,fmaxrefine=5e12,nrefine=2000)
            zpar=z_param(ztypescan=2,zmin=0.1,zmax=1e7,nzlog=20,zminrefine=2e-6,zmaxrefine=0.02,zsamplin=2e-6)
            waketol=1.e12;
            
            iw_input=impedance_wake_input(gamma=gamma,length=1,b=[halfgap],layers=layers_iw,
                    fpar=fpar,zpar=zpar,waketol=waketol,freqlinbisect=1e11,geometry='round',comment='_'+name+strhgap);

            imp_mod_vac,wake_mod_vac=imp_model_elliptic(iw_input,halfgap,orientation='V',
                    wake_calc=wake_calc,flagrm=True,lxplusbatch=lxplusbatchImp,queue='1nh',dire=machine+'_'+scenario_name+'_'+Estr+'/');

            multiply_impedance_wake(imp_mod_vac,circ);
            multiply_impedance_wake(wake_mod_vac,circ);

            # add to the model
            add_impedance_wake(imp_mod,imp_mod_vac,betax/avbetax,betay/avbetay);
            add_impedance_wake(wake_mod,wake_mod_vac,betax/avbetax,betay/avbetay);



            if ((lxplusbatchImp==None) or lxplusbatchImp.startswith('retrieve')):

                imp_mod_list[material].append(imp_mod);
                leg.append('r='+str(halfgap*1e3)+'mm');
                
                # dump into a file
                filemodel=open(root_result+'/impedances'+name+strhgap+'.txt','w');
                pick.dump(imp_mod,filemodel);
                filemodel.close();
                
                # write Ascii files with each component
                write_imp_wake_mod(imp_mod,"_"+machine+"_Allthemachine_"+Estr+'_'+scenario_name,
                        listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
                    dire=root_result+'/')
                write_wake_HEADTAIL(wake_mod,root_result+"/"+machine+"_Allthemachine_"+Estr+'_'+scenario_name+'.wake',ncomp=2);
                # impedance plot
                plot_compare_imp_model([imp_mod],['Vacuum pipe RW impedance'],listcomp=['Zydip'],
                        saveimp=root_result+'/plot_imp_'+name+strhgap,
                        xlim=[1e3,1e11],ylim=[1e4,1e10],legpos=3);

        if ((lxplusbatchImp==None) or lxplusbatchImp.startswith('retrieve')):

            # impedance list plot
            plot_compare_imp_model(imp_mod_list[material][::3],leg[::3],listcomp=['Zydip'],
                    saveimp=root_result+'/plot_imp_'+name,
                    xlim=[1e3,1e11],ylim=[1e4,1e10],legpos=3);

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):
        
        # DELPHI loops now
        tuneshiftQp=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
        tuneshiftm0Qp=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
        tuneshiftQpSach=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex);
        tuneshiftQpSachm0=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(Nbscan),1),dtype=complex);

        for ihalfgap,halfgap in enumerate(hgapscan):
        
            strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

            for iscenario,scenario in enumerate(scenarioscan):

                imp_mod=imp_mod_list[scenario][ihalfgap];

                for iplane,plane in enumerate(['x']):
                    # select Zxdip or Zydip
                    for iw in imp_mod:
                        if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);
                    
                    for iM,M in enumerate(Mscan):

                        # normalization factor for damper
                        dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
                                flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
                        dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

                        flag_trapz=0; # by default no trapz method

                        if (M==1): nxscan=np.array([0]);flag_trapz=1;
                        else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,M/20),np.arange(M/2-10,M/2+11),
                                np.arange(M-10,M),np.arange(0,10))));print("number of coupled-bunch modes=",len(nxscan));

                        tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
                        tuneshiftnxSach=np.zeros((len(Qpscan),len(nxscan),len(Nbscan),1,2*lmaxSach+1),dtype=complex);
                        ZeffSach=np.zeros((len(Qpscan),len(nxscan),1,2*lmaxSach+1),dtype=complex);

                        # DELPHI
                        print("half-gap=",halfgap,", M=",M)
                        tuneshiftQp[ihalfgap,iscenario,iplane,iM,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[ihalfgap,iscenario,iplane,iM,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
                                nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
                                a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
                                flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
                                kmax=kmax,kmaxplot=kmaxplot,crit=1.e-2,abseps=1e-5,flagm0=True,
                                lxplusbatch=lxplusbatchDEL,comment=machine+strhgap+'_holes_'+scenario+'_'+Estr+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
                                queue='1nd',dire=root_result+'/');

                        if flagSach:
                            # Sacherer (no damper)
                            tuneshiftQpSach[ihalfgap,iscenario,iplane,iM,:,:,:,:],tuneshiftnxSach,tuneshiftQpSachm0[ihalfgap,iscenario,iplane,iM,:,:,:],ZeffSach=sacherer(imp_mod,
                                    Qpscan,nxscan,Nbscan,[omegas],M,omega0,eval('Q'+plane),gamma,eta,taub,lmaxSach,
                                    particle='proton',modetype='sinusoidal',compname='Z'+plane+'dip');

        if flagSach:
            # save Sacherer tuneshifts
            fileSach=open(root_result+'/Sacherer_'+Estr+'.txt','w');
            pick.dump(tuneshiftQpSach,fileSach);
            pick.dump(tuneshiftQpSachm0,fileSach);
            fileSach.close();
        else:
            # load Sacherer tuneshifts
            fileSach=open(root_result+'/Sacherer_'+Estr+'.txt','r');
            tuneshiftQpSach=pick.load(fileSach);
            tuneshiftQpSachm0=pick.load(fileSach);
            fileSach.close();
        

        # now the plots (outside loop on scenarios)
        if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

            for iplane,plane in enumerate(['x']):

                for iM,M in enumerate(Mscan):

                    for idamp,damp in enumerate(dampscan):

                        for Nb in Nbscanplot:

                            if False:
                                # plots vs Q'
                                for ihalfgap,halfgap in enumerate(hgapscan):

                                    strhgapleg='halfgap '+str(1e3*halfgap)+'mm';
                                    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

                                    for iscenario,scenario in enumerate(scenarioscan):

                                        # initialize plots vs Qp
                                        figQpm0,axQpm0=init_figure(axes=[0.15,0.1,0.8,0.85]);
                                        figQp=[];axQp=[];
                                        for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

                                        # output file name for plots vs Qp
                                        fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
                                        fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

                                        strpart=['Re','Im'];
                                        for ir,r in enumerate(['real','imag']):

                                            # output file name for data vs Qp
                                            fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
                                            fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

                                            ts=getattr(tuneshiftQp[ihalfgap,iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
                                            Sachstr='';
                                            if damp==0:
                                                # compare with Sacherer most unstable mode
                                                tsSach=getattr(tuneshiftQpSach[ihalfgap,iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0,0],r);
                                                Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift"
                                                data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
                                            else:
                                                data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));

                                            write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\tDELPHI_"+strpart[ir]+"_tuneshift"+Sachstr)

                                            sgn=1;sgnstr='';
                                            if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
                                            plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+strhgapleg+', '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
                                            if damp==0:
                                                plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, '+strhgapleg+', '+scenario,'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");


                                            # real tune shift of mode 0
                                            if (ir==0):
                                                ts=getattr(tuneshiftm0Qp[ihalfgap,iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
                                                Sachstr='';

                                                if damp==0:
                                                    # compare with Sacherer most unstable mode
                                                    tsSach=getattr(tuneshiftQpSachm0[ihalfgap,iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0],r);
                                                    Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift_mode0"
                                                    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));

                                                else:
                                                    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));

                                                write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\tDELPHI_"+strpart[ir]+"_tuneshift_mode0"+Sachstr)

                                                sgn=1;sgnstr='';
                                                if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
                                                plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+strhgapleg+', '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");
                                                if damp==0:
                                                        plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, '+strhgapleg+', '+scenario,'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");

                                        # finish plots vs Qp
                                        if (ir==0):
                                            end_figure(figQpm0,axQpm0,save=flagsave*(fileoutplotQpm0+'_'+r))

                                    for ir,r in enumerate(['real','imag']):
                                        end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))

                            
                            # plot imag tune shifts (in terms of growth rates) of most unstable mode vs half-gap
                            for Qp in Qpplothg:

                                iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];
                                
                                strpart=['Re','Im'];
                                r='real';ir=0;fact=1;ylab=" $ |"+strpart[ir]+"(Q-Q_0)| $ ";
                                r='imag';ir=1;fact=omega0;ylab="Growth rate [1/s]";

                                for iscenario,scenario in enumerate(scenarioscan):

                                    # initialize plots vs half-gap
                                    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);

                                    # output file name for plots vs hg
                                    fileoutplothg=root_result+'/plot_vs_hg_m0_'+machine+'_'+Estr+'_Qp'+float_to_str(Qp)+'_'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

                                    # output file name for data vs halfgap
                                    fileoutdatahg=root_result+'/data_vs_hg_m0_'+machine+'_'+Estr+'_Qp'+float_to_str(Qp)+'_'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

                                    # tune shift of most unstable mode
                                    ts=getattr(tuneshiftQp[:,iscenario,iplane,iM,iQp,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
                                    Sachstr='';
                                    if damp==0:
                                        # compare with Sacherer mode 0
                                        tsSach=getattr(tuneshiftQpSach[:,iscenario,iplane,iM,iQp,pylab.mlab.find(Nbscan==Nb),0,0],r);
                                        Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift"
                                        data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
                                    
                                    else:
                                            data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1))));
                                    
                                    write_ncol_file(fileoutdatahg+'_'+r+'.dat',data,header="halfgap[m]\tDELPHI_"+strpart[ir]+"_tuneshift_mode0"+Sachstr)

                                    sgn=1;sgnstr='-';
                                    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
                                    
                                    plot(hgapscan*1e3,np.abs(np.squeeze(sgn*ts))*fact,'DELPHI',col[iscenario],ylab,ax,2,xlab=" half-gap [mm] ");
                                    if damp==0:
                                            plot(hgapscan*1e3,np.abs(np.squeeze(sgn*tsSach))*fact,'Sacherer','--'+col[iscenario],ylab,ax,2,xlab=" half-gap [mm] ");

                                    # add on top of the plot some lines for the feedback gain (100, 50 and 10 turns)
                                    for idamping,damping in enumerate([10.,20.,50.]):
                                        dampstr=str(int(damping))+' turns damping';
                                        plot(hgapscan*1e3,f0/damping*np.ones(len(hgapscan)),dampstr,'-'+col[iscenario+idamping+1],ylab,ax,2,xlab=" half-gap [mm] ");

                                    # finish plots vs half-gap
                                    end_figure(fig,ax,save=flagsave*(fileoutplothg+'_'+r))

                        # TMCI plots
                        Nbthres=np.zeros((len(hgapscan),len(scenarioscan),len(Qpplothg)))
                        for Qp in Qpplothg:

                            iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];

                            for ihalfgap,halfgap in enumerate(hgapscan):

                                strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

                                for iscenario,scenario in enumerate(scenarioscan):

                                    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+Estr+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(round(10.*Qp)/10.)+'_converged'+strnorm[flagnorm]+'_'+plane;
                                    patcol=['.b','b'];
                                    ylim=([-5,3],[-0.001,0.02]);

                                    for ir,r in enumerate(['real','imag']):

                                        fig,ax=init_figure();

                                        ts=tuneshiftQp[ihalfgap,iscenario,iplane,iM,iQp,idamp,:,0,0,:];

                                        plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg='DELPHI',patcol=patcol[ir],xlab='Nb [p+/b]',
                                            title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.)+", half-gap="+str(halfgap*1e3)+"mm, "+scenario,ms=1,ylim=ylim[ir]);

                                        end_figure(fig,ax,save=flagsave*(fileoutplotTMCI+'_'+r),fontsize=25);

                                    # initialize plot of threshold vs half-gap
                                    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);
                                    # output file name for threshold plot vs hg
                                    fileoutplothg=root_result+'/plot_vs_hg_thres_'+machine+'_'+Estr+scenario+'_Qp'+float_to_str(Qp)+'_'+str(M)+'b_d'+float_to_str(damp)+'_converged'+strnorm[flagnorm]+'_'+plane;

                                    # find intensity threshold
                                    Nbthres[ihalfgap,iscenario,iQp]=find_intensity_threshold(Nbscan,tuneshiftQp[ihalfgap,iscenario,iplane,iM,iQp,idamp,:,0,0,0]*omega0,thresgrowth=1.);
                                    plot(hgapscan*1e3,np.squeeze(Nbthres[:,iscenario,iQp]),'DELPHI, '+scenario,col[iscenario],"TMCI threshold (nb p/b)",ax,2,xlab=" half-gap [mm] ");
                                    ax.set_xlim([5,20]);ax.set_ylim([1e9,1e12]);
                                
                                    # finish plots vs half-gap
                                    end_figure(fig,ax,save=flagsave*(fileoutplothg+'_'+r))

    if not(flagsave): pylab.show();            #Lh = 8e-3; Wh = 1.5e-3; T = 1.075e-3; b =18.375e-3; d = 23.175e-3; eta = 0.044; rhob = 6e-7; rhod = 6e-7;nb_holes_per_cs = 8;fcutoff=5E9; # LHCparameters
