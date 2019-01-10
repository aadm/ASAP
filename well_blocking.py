# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 10:35:05 2018

@author: ANDD
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
import matplotlib 
#%matplotlib qt
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)

Wnam = '15_5_5';top_blocks = [2050,2154, 2200,2245]
#Wnam = '155_6';top_blocks = [2056,2160, 2222,2267]
#Wnam = '155_3';top_blocks = [2076,2180, 2240,2285]
infile = Wnam + '.mat'
Wdict = sio.loadmat(infile)

keys = ['Md','Vp_b', 'Vs_b', 'Rhob_b', 'porosity', 'Vsh', 'Vp_o', 'Vs_o', 'Rhob_o']
cleaned_wdict = {key: Wdict[key].ravel() for key in keys}
df_in = pd.DataFrame.from_dict(cleaned_wdict)
df_in = df_in.dropna()

Md = df_in.Md

Wdict=df_in
def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]
  
def logplot(Vsh,Md,Vp,Vs,Rho,MD,Vp_out,Vs_out,Rho_out):   
    plt.figure(figsize = (20,10))
    plt.subplot(1,4,1)
    plt.plot(Vsh,Md)
    plt.ylim(np.max(MD)+2,np.min(MD)-2)
    plt.xlabel('Vsh (fraction)')
    plt.ylabel('Measured depth (m)')
    
    plt.subplot(1,4,2)
    plt.plot(Vp,Md)
    plt.plot(Vp_out,MD, color = 'red',linewidth=4.0 )
    plt.ylim(np.max(MD)+2,np.min(MD)-2)
    plt.xlabel('Vp (m/s)')
    plt.gca().set_yticks([])
    
    plt.subplot(1,4,3)
    plt.plot(Vs,Md)
    plt.plot(Vs_out,MD, color = 'red' ,linewidth=4.0)
    plt.ylim(np.max(MD)+2,np.min(MD)-2)
    plt.xlabel('Vs (m/s)')
    plt.gca().set_yticks([])
    
    plt.subplot(1,4,4)
    plt.plot(Rho,Md)
    plt.plot(Rho_out,MD, color = 'red',linewidth=4.0 )
    plt.ylim(np.max(MD)+2,np.min(MD)-2)
    plt.xlabel('Density (g/ccm)')
    plt.gca().set_yticks([])
    
def generate_wells(Wdict,top_blocks,Wnam):
    
    Vp = Wdict['Vp_b'].ravel();
    Vs = Wdict['Vs_b'].ravel()
    Rho = Wdict['Rhob_b'].ravel()
    Por = Wdict['porosity'].ravel()
    Vsh = Wdict['Vsh'].ravel()
    Md = Wdict['Md'].ravel()
    Vp_o = Wdict['Vp_o'].ravel()
    Vs_o = Wdict['Vs_o'].ravel()
    Rho_o = Wdict['Rhob_o'].ravel()

    Vp_bl = np.zeros(np.size(Vp))
    Vs_bl=np.zeros(np.size(Vp))
    Rho_bl = np.zeros(np.size(Vp))
    
    
    Vp_const = np.zeros((np.shape(top_blocks)))
    Vs_const = np.zeros(np.shape(top_blocks))
    Rho_const = np.zeros(np.shape(top_blocks))
    Por_const = np.zeros(np.shape(top_blocks))
    
    Vp_const_o = np.zeros((np.shape(top_blocks)))
    Vs_const_o = np.zeros(np.shape(top_blocks))
    Rho_const_o = np.zeros(np.shape(top_blocks))
    
    for ii, top in enumerate(top_blocks):
        if ii== len(top_blocks)-1:
            continue
        ind = indices(Md, lambda x: x>= top_blocks[ii] and x< top_blocks[ii+1])
        Vp_bl[ind] = np.mean(Vp[ind])
        Vs_bl[ind] = np.mean(Vs[ind])
        Rho_bl[ind] = np.mean(Rho[ind])
       
        Vp_const[ii]=np.mean(Vp[ind]) 
        Vs_const[ii]=np.mean(Vs[ind]) 
        Rho_const[ii]=np.mean(Rho[ind])
        Por_const[ii]=np.mean(Por[ind])
        
        Vp_const_o[ii]=np.mean(Vp_o[ind]) 
        Vs_const_o[ii]=np.mean(Vs_o[ind]) 
        Rho_const_o[ii]=np.mean(Rho_o[ind]) 
        
        
    
    

        
    d = {'Depth':Md, 'Vp': Vp_bl, 'Vs': Vs_bl, 'Rho': Rho_bl}
    df = pd.DataFrame(data=d)
    df2 = df[(df.Depth>=top_blocks[0]) & (df.Depth<=top_blocks[3])]
#    df2.to_csv('W15_5_5')
    
    Z= np.linspace(0,60,4)
    Sw= np.linspace(0,1,5)
    dPor = [-.1, -.05, 0, .05]
    MD = df2.Depth.ravel()
     
    ind_top = indices(MD, lambda x: x>= top_blocks[0] and x<= top_blocks[1])
    for jj in Z:
        Vp_z = np.ones(np.shape(MD))*Vp_const[2]
        Vs_z = np.ones(np.shape(MD))*Vs_const[2]
        Rho_z = np.ones(np.shape(MD))*Rho_const[2]
        
        Vp_z_o = np.ones(np.shape(MD))*Vp_const_o[2]
        Vs_z_o = np.ones(np.shape(MD))*Vs_const_o[2]
        Rho_z_o = np.ones(np.shape(MD))*Rho_const_o[2]
        
        ind_base = indices(MD, lambda x: x> top_blocks[1] and x< top_blocks[1]+jj)
        Vp_z[ind_top]=Vp_const[0];Vp_z[ind_base]=Vp_const[1]
        Vs_z[ind_top]=Vs_const[0];Vs_z[ind_base]=Vs_const[1]
        Rho_z[ind_top]=Rho_const[0];Rho_z[ind_base]=Rho_const[1]
        
        Vp_z_o[ind_top]=Vp_const_o[0];Vp_z_o[ind_base]=Vp_const_o[1]
        Vs_z_o[ind_top]=Vs_const_o[0];Vs_z_o[ind_base]=Vs_const_o[1]
        Rho_z_o[ind_top]=Rho_const_o[0];Rho_z_o[ind_base]=Rho_const_o[1]
        
        
        
        for kk in Sw:
            Vp_sw = kk*Vp_z + (1-kk)*Vp_z_o
            Vs_sw = kk*Vs_z + (1-kk)*Vs_z_o
            Rho_sw = kk*Rho_z + (1-kk)*Rho_z_o
            
            
            for pp in dPor:
                dVp_br = -pp*4500
                dVp_o = -pp*4500
                dVs_br = -pp*3100
                dVs_o = -pp*3050
                dVs = kk*dVs_br + (1-kk)*dVs_o
                dVp = kk*dVp_br + (1-kk)*dVp_o
                
                dRho_br = (2.65 - 1)*-pp
                dRho_o  = (2.65 -0.8)*-pp
                dRho = kk*dRho_br + (1-kk)*dRho_o
                
                #Vp_out=Vp_sw; Vs_out=Vs_sw; Rho_out=Rho_sw; 
                Vp_out=np.copy(Vp_sw); Vs_out=np.copy(Vs_sw); Rho_out=np.copy(Rho_sw);
                Vp_out[ind_base] = Vp_sw[ind_base] + dVp
                Vs_out[ind_base] = Vs_sw[ind_base] + dVs
                Rho_out[ind_base] = Rho_sw[ind_base] + dRho
               
#                if jj==60 and kk==1 and pp>0:
#                    print(Rho_out.shape, MD.shape, Rho_sw.shape)
#                    print('dRho is: ' +str(dRho),' and dPor is:  ',str(pp))
#                    logplot(Vsh,Md,Vp,Vs,Rho,MD,Vp_out,Vs_out,Rho_out)
            
                d_out = {'Depth':MD, 'Vp': Vp_out, 'Vs': Vs_out, 'Rho': Rho_out}
                df_out = pd.DataFrame(data=d_out)
#                name = Wnam + '_Z' + str(int(np.round(jj))) + '_Sw' + str(int(100*kk)) + '_Por' + str(int((0.3+pp)*100))
                name = 'W{:s}_Z{:02d}_Sw{:03d}_Por{:02d}'.format(Wnam,int(np.round(jj)),int(100*kk),int((0.3+pp)*100))
#                name.replace('Z0','Z00')
#                name.replace('Sw0','Sw00')
                df_out.to_csv(name)
                if pp==-.10 and jj == 40 and kk == 1:
                    logplot(Vsh,Md,Vp,Vs,Rho,MD,Vp_out,Vs_out,Rho_out)
                    plt.savefig('Pseudo_well_logs_porchanve.png')

generate_wells(df_in,top_blocks,Wnam)       
    
 
    
