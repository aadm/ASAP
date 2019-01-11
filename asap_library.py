import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import bruges as b
from scipy.interpolate import interp1d
import os

topbox = dict(boxstyle='round', ec='none', fc='w', alpha=0.6)
format_tops={'fontsize':10, 'color':'blue', 'ha':'right', 'bbox':topbox}
format_title={'fontsize':14, 'weight':'bold'}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
def contactcement(K0, G0, phi, phi_c=0.4, Cn=8.6, Kc=37, Gc=45, scheme=2):
    PR0=(3*K0-2*G0)/(6*K0+2*G0)
    PRc = (3*Kc-2*Gc)/(6*Kc+2*Gc)
    if scheme == 1: # scheme 1: cement deposited at grain contacts
        alpha = ((phi_c-phi)/(3*Cn*(1-phi_c))) ** (1/4)
    else: # scheme 2: cement evenly deposited on grain surface
        alpha = ((2*(phi_c-phi))/(3*(1-phi_c)))**(1/2)
    LambdaN = (2*Gc*(1-PR0)*(1-PRc)) / (np.pi*G0*(1-2*PRc))
    N1 = -0.024153*LambdaN**-1.3646
    N2 = 0.20405*LambdaN**-0.89008
    N3 = 0.00024649*LambdaN**-1.9864
    Sn = N1*alpha**2 + N2*alpha + N3
    LambdaT = Gc/(np.pi*G0)
    T1 = -10**-2*(2.26*PR0**2+2.07*PR0+2.3)*LambdaT**(0.079*PR0**2+0.1754*PR0-1.342)
    T2 = (0.0573*PR0**2+0.0937*PR0+0.202)*LambdaT**(0.0274*PR0**2+0.0529*PR0-0.8765)
    T3 = 10**-4*(9.654*PR0**2+4.945*PR0+3.1)*LambdaT**(0.01867*PR0**2+0.4011*PR0-1.8186)
    St = T1*alpha**2 + T2*alpha + T3
    K_DRY = 1/6*Cn*(1-phi_c)*(Kc+(4/3)*Gc)*Sn
    G_DRY = 3/5*K_DRY+3/20*Cn*(1-phi_c)*Gc*St
    return K_DRY, G_DRY

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
def constantcement(K0, G0, phi, phi_cem=0.38, phi_c=0.4, Cn=8.6, Kc=37, Gc=45, scheme=2):
    # contact cement model
    K_HI, G_HI = contactcement(K0, G0, phi, phi_c=phi_c, Cn=Cn, Kc=Kc, Gc=Gc, scheme=scheme)
    # lower bound Hashin-Shtrikman starting from phi_cem
    Kcc, Gcc = contactcement(K0, G0, phi_cem, phi_c=phi_c, Cn=Cn, Kc=Kc, Gc=Gc, scheme=scheme)
    K_LO = -4/3*Gcc + (((phi/phi_cem)/(Kcc+4/3*Gcc)) + ((1-phi/phi_cem)/(K0+4/3*Gcc)))**-1
    tmp = Gcc/6*((9*Kcc+8*Gcc) / (Kcc+2*Gcc))
    G_LO= -tmp + ((phi/phi_cem)/(Gcc+tmp) + ((1-phi/phi_cem)/(G0+tmp)))**-1
    # initialize empty vectors for K and G dry
    K_DRY, G_DRY=(np.full(phi.size, np.nan) for _ in range(2))
    # for porosities>phi_cem use [K,G]_HI = contact cement model
    # for porosities<=phi_cem use [K,G]_LO = constant cement model
    K_DRY[phi>phi_cem]=K_HI[phi>phi_cem]
    K_DRY[phi<=phi_cem]=K_LO[phi<=phi_cem]
    G_DRY[phi>phi_cem]=G_HI[phi>phi_cem]
    G_DRY[phi<=phi_cem]=G_LO[phi<=phi_cem]
    return K_DRY, G_DRY

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
def inccement(K0, G0, phi, phi_cem=0.38, phi_c=0.4, Cn=8.6, Kc=37, Gc=45, scheme=2):
    Kcc, Gcc = contactcement(K0, G0, phi_cem, phi_c=phi_c, Cn=Cn, Kc=Kc, Gc=Gc, scheme=scheme)
    K_DRY = -4/3*G0 + (((phi/phi_cem)/(Kcc+4/3*G0)) + ((1-phi/phi_cem)/(K0+4/3*G0)))**-1
    tmp = G0/6*((9*K0+8*G0) / (K0+2*G0))
    G_DRY = -tmp + ((phi/phi_cem)/(Gcc+tmp) + ((1-phi/phi_cem)/(G0+tmp)))**-1
    return K_DRY, G_DRY

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
def vels(k_dry,mu_dry,k_min,rho_min,k_fl,rho_fl,phi):
    # converts all inputs to SI (density in kg/m3 and moduli in Pa)
    KD = k_dry*1e9
    GD = mu_dry*1e9
    K0 = k_min*1e9
    D0 = rho_min*1e3
    Kf = k_fl*1e9
    Df = rho_fl*1e3
    rho = D0*(1-phi)+Df*phi
    K = KD + (1-KD/K0)**2 / ( (phi/Kf) + ((1-phi)/K0) - (KD/K0**2) )
    vp = np.sqrt((K+4/3*GD)/rho)
    vs = np.sqrt(GD/rho)
    return vp, vs, rho/1e3, K/1e9

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
def blok_plot(REF,BLK,ztop=None,zbot=None):
    if ztop is None: ztop = BLK.index.min()
    if zbot is None: zbot = BLK.index.max()
    REF=REF[(REF.index>=ztop) & (REF.index<=zbot)]
    REF.columns=REF.columns.str.upper()
    BLK.columns=BLK.columns.str.upper()

    f, ax = plt.subplots(nrows=1,ncols=5,sharey=True,figsize=(8,5))
    ax[0].plot(REF.VSH, REF.index, color='.5')
    ax[0].set_xlabel('Vsh', color='.5')
    ax[1].plot(REF.PHIE, REF.index, color='.5')
    ax[1].set_xlabel('phi', color='.5')
    ax[2].plot(REF.VP_FRMB, REF.index, color='black')
    ax[2].plot(BLK.VP, BLK.index, color='red', lw=4)
    ax[2].set_xlabel('Vp [m/s]')
    ax[3].plot(REF.VS_FRMB, REF.index, color='black')
    ax[3].plot(BLK.VS, BLK.index, color='red', lw=4)
    ax[3].set_xlabel('Vs [m/s]')
    ax[4].plot(REF.RHO_FRMB, REF.index, color='black')
    ax[4].plot(BLK.RHO, BLK.index, color='red', lw=4)
    ax[4].set_xlabel('Density [g/cc]')
    ax[0].set_ylim(zbot,ztop)   
    plt.tight_layout()
 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def generate_wells(WW,top_blocks,name,output_dir):
    # top_blocks is an array with 4 values defining 3 blocks:
    # block 1 = shallow block, caprock
    # block 2 = central block, reservoir
    # block 3 = deep block
    z = WW.index.values
    nblks = np.size(top_blocks)-1

    # compute and store average vp,vs,rho for fluids 1 (brine) and 2 (oil)
    vp1_k,vs1_k,rho1_k,vp2_k,vs2_k,rho2_k,phi_k = (np.zeros(nblks) for _ in range(7))
    for nn in range(nblks):
        rr = (WW.index>=top_blocks[nn]) & (WW.index<top_blocks[nn+1])
        vp1_k[nn] = WW.VP_FRMB[rr].mean()
        vs1_k[nn] = WW.VS_FRMB[rr].mean()
        rho1_k[nn] = WW.RHO_FRMB[rr].mean()
        vp2_k[nn] = WW.VP_FRMO[rr].mean()
        vs2_k[nn] = WW.VS_FRMO[rr].mean()
        rho2_k[nn] = WW.RHO_FRMO[rr].mean()
    
    # range of variations for thickness (vTH), water saturation (vSW) and porosity (vPH)
    vTH= np.linspace(0,60,4)
    vSW= np.linspace(0,1,5)
    vPH = [-.1, -.05, 0, .05]

    # output depths will be within top block 1 and base block 3
    z_out = z[(z>=top_blocks[0]) & (z<=top_blocks[3])]
    # mm0 = to select depths above reservoir
    mm0 = (z_out>=top_blocks[0]) & (z_out<top_blocks[1])

    # loop over thickness variations
    for tt in vTH:
        # initialize vp,vs,rho for both fluids with values taken from block 3
        vp1_TH  = np.ones(z_out.size)*vp1_k[2]
        vs1_TH  = np.ones(z_out.size)*vs1_k[2]
        rho1_TH = np.ones(z_out.size)*rho1_k[2]
        vp2_TH  = np.ones(z_out.size)*vp2_k[2]
        vs2_TH  = np.ones(z_out.size)*vs2_k[2]
        rho2_TH = np.ones(z_out.size)*rho2_k[2]

        # sets vp,vs,rho for block 1
        vp1_TH[mm0] = vp1_k[0]
        vs1_TH[mm0] = vs1_k[0]
        rho1_TH[mm0] = rho1_k[0]
        vp2_TH[mm0] = vp2_k[0]
        vs2_TH[mm0] = vs2_k[0]
        rho2_TH[mm0] = rho2_k[0]

        # mm1 = to select depths in block 2=reservoir (base varies with thickness variations)
        mm1 = (z_out>=top_blocks[1]) & (z_out<top_blocks[1]+tt)

        # sets vp,vs,rho for block 2
        vp1_TH[mm1] = vp1_k[1]
        vs1_TH[mm1] = vs1_k[1]
        rho1_TH[mm1] = rho1_k[1]
        vp2_TH[mm1] = vp2_k[1]
        vs2_TH[mm1] = vs2_k[1]
        rho2_TH[mm1] = rho2_k[1]
        # loop over water saturation variations
        for ss in vSW:
            # modify velocities and density according to saturation changes
            # assume a linear change from brine to oil filled rock
            vp_SW = ss*vp1_TH + (1-ss)*vp2_TH
            vs_SW = ss*vs1_TH + (1-ss)*vs2_TH
            rho_SW = ss*rho1_TH + (1-ss)*rho2_TH

            # loop over porosity variations
            for pp in vPH:
                # modify velocities and density according to porosity changes
                # assume linear porosity-velocity ratios in the limited range of study
                # relations are given by rock physics study
                dvp_b = -pp*4500
                dvp_o = -pp*4500
                dvs_b = -pp*3100
                dvs_o = -pp*3050
                dvs   = ss*dvs_b + (1-ss)*dvs_o
                dvp   = ss*dvp_b + (1-ss)*dvp_o
                drho_b = (2.65-1.0)*-pp
                drho_o = (2.65-0.8)*-pp
                drho   = ss*drho_b + (1-ss)*drho_o
                # calculate new parameters with varying saturations 
                vp_out=np.copy(vp_SW); vs_out=np.copy(vs_SW); rho_out=np.copy(rho_SW);
                vp_out[mm1] = vp_SW[mm1] + dvp
                vs_out[mm1] = vs_SW[mm1] + dvs
                rho_out[mm1] = rho_SW[mm1] + drho
                d_out = {'DEPTH':z_out, 'VP': vp_out, 'VS': vs_out, 'RHO': rho_out}
                df_out = pd.DataFrame(data=d_out)
                filename = '{:s}_Z{:02d}_Sw{:03d}_Por{:02d}'.format(name,int(np.round(tt)),int(100*ss),int((0.3+pp)*100))
                df_out.to_csv(output_dir+'/'+filename, index=False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_well_files(wells_dir, name):
   well_files = []
   for dirpath, _, filenames in os.walk(wells_dir):
       well_files += [os.path.join(dirpath, f) for f in filenames if f.startswith(name)]
   return well_files

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_nears_fars(df):
   sample_size = 64
   number_of_splits = abs(df['NEAR'].size/64)
   nears = np.array_split(df['NEAR'].values, number_of_splits)
   fars = np.array_split(df['FAR'].values, number_of_splits)
   nears = np.asarray(nears).transpose()
   fars = np.asarray(fars).transpose()
   return nears, fars

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_twt(tdr,z):
    tt=tdr[:,1]
    zz=tdr[:,0] 
    d2t = interp1d(zz, tt, kind='linear', bounds_error=False, fill_value='extrapolate')
    return d2t(z)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_synt(WW,ang,wavelet,method='shuey'):
    '''
    WW: Pandas dataframe with VP, VS, RHO (optionally EPSILON and DELTA)
    ang: angle range, define with ang=np.arange(0,50,1)
    wavelet
    method: 'shuey' (Shuey 2-terms), 'shuey3' (Shuey 3-terms),
            'aki' (Aki-Richards)
    '''
    WW.columns=WW.columns.str.upper()
    uvp, lvp   = WW.VP.values[:-1], WW.VP.values[1:]
    uvs, lvs   = WW.VS.values[:-1], WW.VS.values[1:]
    urho, lrho = WW.RHO.values[:-1], WW.RHO.values[1:]
    z=WW.index.values  # z is two-way-time
    synt  = np.zeros((z.size,ang.size))
    #--> calculate reflectivities with AVO equation,
    #--> convolve with input wavelet and fill in traces of synthetic seismogram
    for i,alpha in enumerate(ang):
        if method is 'shuey':
            RC = shuey(uvp,uvs,urho,lvp,lvs,lrho,alpha)
        elif method is 'shuey3':
            RC = shuey(uvp,uvs,urho,lvp,lvs,lrho,alpha,approx=False)
        else:
            RC = akirichards(uvp,uvs,urho,lvp,lvs,lrho,alpha)
        RC = np.append(np.nan, RC)
        RC = np.nan_to_num(RC)
        synt[:,i] = np.convolve(RC, wavelet, mode='same')
    return RC, synt

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_synt(WW,synt,ztop,zbot,gain=10):
    '''
    WW: Pandas dataframe with VP, VS, RHO (in time)
    synt: synthetic seismogram computed with make_synt (in time)
    ztop,zbot: display window
    gain: multiplier to be applied to wiggles (default=5)
    method: 'shuey' (Shuey 2-terms), 'shuey3' (Shuey 3-terms),
            'aki' (Aki-Richards)
    '''
    WW.columns=WW.columns.str.upper()
    it1=np.abs(WW.index-ztop).argmin()
    it2=np.abs(WW.index-zbot).argmin()
    ss = synt[it1:it2,:]
    clip=np.abs(synt.max())
    f,ax=plt.subplots(nrows=1,ncols=5)
    opz1={'color':'k','linewidth':.5}
    opz2={'linewidth':0, 'alpha':0.6}
    ax[0].plot(WW.VP*WW.RHO,WW.index,'-k')
    ax[1].plot(WW.VP/WW.VS,WW.index,'-k')
    ax[2].plot(synt[:,0],WW.index,'-k')
    ax[3].plot(synt[:,1],WW.index,'-k')
    im=ax[4].imshow(ss,interpolation='none', cmap='Greys',aspect='auto')
    cbar=plt.colorbar(im, ax=ax[4])
    ax[0].set_xlabel('AI [m/s*g/cc]')
    ax[0].set_ylabel('TWT [s]')
    ax[1].set_xlabel('Vp/Vs')
    ax[2].set_title('Near')
    ax[3].set_title('Far')
    ax[4].set_title('Near|Far')
    for aa in ax[:4]:
        aa.set_ylim(zbot,ztop)
        aa.grid()
    for aa in ax[1:]:
        aa.set_yticklabels([])
    for aa in ax[2:4]:
        aa.set_xlim(-clip,+clip)
        aa.set_xticklabels([])
    for aa in ax[:2]:
        aa.xaxis.tick_top()
        plt.setp(aa.xaxis.get_majorticklabels(), rotation=90)
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def td(WW,sr=0.1524,KB=0,WD=0,repl_vel=1600):
    '''
    td (C) aadm 2016
    Calculates time-depth relation by sonic integration.

    INPUT
    WW: Pandas dataframe
    sr: depth sampling rate (m, default 0.1524 m = half-foot)
    KB: kelly bushing elevation (m, default 0)
    WD: water depth (m, default 0)
    repl_vel: replacement velocity for overburden i.e. interval between seafloor and beginning of the logs
              (m/s, default 1600)

    OUTPUT
    numpy array with 2 columns: 0=depth, 1=twt (secs)
    '''
    WW.columns=WW.columns.str.upper()
    if 'TVD' in WW.columns:
        depth = WW.TVD.values
    else:
        depth = WW.index.values
    WW.VP.interpolate(inplace=True)
    sonic=1/WW.VP.values # VP in m/s
    start = depth.min()
    water_vel = 1480
    wb_twt = 2.0*WD/water_vel
    sonic_start=depth[np.isfinite(sonic)].min()
    sonic_start_twt=2*(sonic_start-KB-WD)/repl_vel + wb_twt
    scaled_sonic = sr*sonic[depth>=sonic_start]
    twt = 2*np.cumsum(scaled_sonic) + sonic_start_twt
    print('[TD] water bottom two-way-time: {:.3f} [s]'.format(wb_twt))
    print('[TD] sonic log start: {:.3f} [m] = {:.3f} [s]'.format(sonic_start, sonic_start_twt))
    print('[TD] computed twt scale range: {:.3f}-{:.3f} [s]'.format(twt.min(),twt.max()))
    return np.column_stack((depth[depth>=sonic_start],twt))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def welltime(WW,tdr,dt=0.001,ztop=None,zbot=None,name=None,tops=None,qcplot=True):
    '''
    welltime (C) aadm 2016
    Converts logs sampled in depth to time using a reference time-depth function.
    Load existing t-d with tdr=np.loadtxt('TD.dat', delimiter=',')
    or use squit.well.td to create one.

    INPUT
    WW: Pandas dataframe
    tdr: time-depth table: numpy array shape Nx2, column 0 = MD [m], column 1 = TWT [s]
    dt: sample rate in seconds
    ztop,zbot: display window (defaults to min,max depth)
    name: well name (or anything else) to print
    tops: dictionary containing stratigraphic tops, e.g.: tops={'Trias': 1200,'Permian': 2310}
          or Pandas Series, e.g: tops=pd.Series({'Trias': 1200,'Permian': 2310})

    OUTPUT
    Pandas dataframe with logs sampled in time
    ztop, zbot converted to TWT
    '''
    WW.columns=WW.columns.str.upper()
    flagfrm=True if 'VP_FRMB' in WW.columns else False
    z = WW.index
    if ztop is None: ztop = z.min()
    if zbot is None: zbot = z.max()
    flagtops=False if tops is None else True
    if flagtops:
        if not isinstance(tops, pd.Series):
            tops=pd.Series(tops)
        tops=tops.dropna().sort_values()

    #--> load depth-time relationship (depth is MD)
    tt=tdr[:,1]
    zz=tdr[:,0]
    #-->  twt reference log sampled like depth reference log
    # twt_log = np.interp(z, zz, tt, left=np.nan, right=np.nan)
    ff=(z>=zz.min()) & (z<=zz.max())
    twt_log = np.interp(z[ff], zz, tt, left=np.nan, right=np.nan)
    #-->  interpolant to convert depths to times on the fly (e.g., tops)
    d2t = interp1d(zz, tt, kind='linear', bounds_error=False, fill_value='extrapolate')
    if qcplot:
        print('[WELLTIME] plot window top, bottom [m]:{:.0f}-{:.0f}, [s]:{:.4f}-{:.4f}'.format(ztop,zbot,float(d2t(ztop)),float(d2t(zbot))))
    #-->  regularly-sampled twt scale and its depth (MD) equivalent on the basis of depth-time rel.
    twt = np.arange(0, tt.max(), dt)
    zt = np.interp(x=twt, xp=tt, fp=zz, left=np.nan, right=np.nan)

    #-->  resample logs to twt
    WWt=pd.DataFrame(data=zt, columns=['DEPTH'], index=twt)
    WWt.index.rename('TWT',inplace=True)
    loglist=WW.columns
    for i in loglist:
        tmp = np.interp(x=twt, xp=twt_log, fp=WW[i][ff].values, left=np.NaN, right=np.NaN)
        WWt=pd.concat([WWt,pd.Series(tmp, index=twt, name=i)],axis=1)
        WWt.interpolate(inplace=True)
        WWt.fillna(method = 'bfill',inplace=True)

    #--> QC plot with IP in depth and time
    if qcplot:
        tmp_IP = WW['VP']*WW['RHO']
        tmp_IPt = WWt['VP']*WWt['RHO']
        plotmax = tmp_IP[(z>=ztop) & (z<=zbot)].max()
        plotmin = tmp_IP[(z>=ztop) & (z<=zbot)].min()
        plotmax += plotmax*.1
        plotmin -= plotmin*.1

        f, ax = plt.subplots(nrows=1,ncols=2,figsize=(5,5), facecolor='w')
        ax[0].plot(tmp_IP, z, '-k')
        ax[1].plot(tmp_IPt, WWt.index, '-k')
        ax[0].set_xlabel('AI [m/s*g/cc]'), ax[0].set_ylabel('MD [m]')
        ax[1].set_xlabel('AI [m/s*g/cc]'), ax[1].set_ylabel('TWT [s]')
        ax[1].yaxis.set_label_position('right')
        ax[1].yaxis.set_ticks_position('right')
        ax[0].set_ylim(ztop,zbot)
        ax[1].set_ylim(d2t(ztop),d2t(zbot))
        for aa in ax.flatten():
            aa.invert_yaxis()
            aa.grid()
            aa.xaxis.tick_top()
            plt.setp(aa.xaxis.get_majorticklabels(), rotation=90, fontsize=8)
            # aa.set_xlim(plotmin,plotmax)
        if flagtops: # plot top markers on all columns
            for topn,topz in tops.iteritems():
                if (topz>=ztop) & (topz<=zbot):
                    ax[0].axhline(y=topz,color=color_top,alpha=alpha_top)
                    ax[0].text(x=plotmax,y=topz,s=topn,**format_tops)
                    ax[1].axhline(y=d2t(topz),color=color_top,alpha=alpha_top)
                    ax[1].text(x=plotmax,y=d2t(topz),s=topn,**format_tops)
        if name is not None:
            plt.suptitle(name, **format_title)
        plt.tight_layout()
    return WWt,float(d2t(ztop)),float(d2t(zbot))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta, approx=True, terms=False):
    '''
    shuey (C) aadm 2016
    Calculates P-wave reflectivity with Shuey's equation

    reference:
    Avseth et al. (2005), Quantitative Seismic Interpretation, Cambridge University Press (p.182)

    INPUT
    vp1, vs1, rho1: P-, S-wave velocity (m/s) and density (g/cm3) of upper medium
    vp2, vs2, rho2: P-, S-wave velocity (m/s) and density (g/cm3) of lower medium
    theta: angle of incidence (degree)
    approx: returns approximate (2-term) form (default: True)
    terms: returns reflectivity, intercept and gradient (default: False)

    OUTPUT
    reflectivity (and optionally intercept, gradient; see terms option) at angle theta
    '''
    a = np.radians(theta)
    dvp = vp2-vp1
    dvs = vs2-vs1
    drho = rho2-rho1
    vp  = np.mean([vp1,vp2])
    vs  = np.mean([vs1,vs2])
    rho = np.mean([rho1,rho2])
    R0 = 0.5*(dvp/vp + drho/rho)
    G  = 0.5*(dvp/vp) - 2*(vs**2/vp**2)*(drho/rho+2*(dvs/vs))
    F =  0.5*(dvp/vp)
    if approx:
        R = R0 + G*np.sin(a)**2
    else:
        R = R0 + G*np.sin(a)**2 + F*(np.tan(a)**2-np.sin(a)**2)
    if terms:
        return R,R0,G
    else:
        return R

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def akirichards(vp1, vs1, rho1, vp2, vs2, rho2, theta):
    '''
    Aki-Richards (C) aadm 2017
    Calculates P-wave reflectivity with Aki-Richards approximate equation
    only valid for small layer contrasts.

    reference:
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.182)

    INPUT
    vp1, vs1, rho1: P-, S-wave velocity (m/s) and density (g/cm3) of upper medium
    vp2, vs2, rho2: P-, S-wave velocity (m/s) and density (g/cm3) of lower medium
    theta: angle of incidence (degree)

    OUTPUT
    reflectivity at angle theta
    '''
    a = np.radians(theta)
    p = np.sin(a)/vp1
    dvp = vp2-vp1
    dvs = vs2-vs1
    drho = rho2-rho1
    vp  = np.mean([vp1,vp2])
    vs  = np.mean([vs1,vs2])
    rho = np.mean([rho1,rho2])
    A = 0.5*(1-4*p**2*vs**2)*drho/rho
    B = 1/(2*np.cos(a)**2) * dvp/vp
    C = 4*p**2*vs**2*dvs/vs
    R = A + B - C
    return R

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def quicklook(WW,ztop=None,zbot=None,name=None,tops=None):
    '''
    quicklook (C) aadm 2015-2018
    Summary well plot with raw and processed logs.

    INPUT
    WW: Pandas dataframe with VP, [VS], RHO, IP, [VPVS or PR], [SWE], PHIE, VSH
    ztop,zbot: depth range to plot (defaults to min,max depth)
    name: well name (or anything else) to print
    tops: dictionary containing stratigraphic tops, e.g.: tops={'Trias': 1200,'Permian': 2310}
          or Pandas Series, e.g: tops=pd.Series({'Trias': 1200,'Permian': 2310})
    '''
    if ztop is None: ztop = WW.index.min()
    if zbot is None: zbot = WW.index.max()
    WW=WW[(WW.index>=ztop) & (WW.index<=zbot)]
    WW.columns=WW.columns.str.upper()
 
    flagtops=False if tops is None else True
    if flagtops:
        if not isinstance(tops, pd.Series):
            tops=pd.Series(tops)
        tops=tops.dropna().sort_values()

    f, ax = plt.subplots(nrows=1,ncols=4,sharey=True,figsize=(8,5))
    ax[0].plot(WW.VSH, WW.index, color='.5')
    ax[0].set_xlabel('Vsh', color='.5')
    ax[1].plot(WW.PHIE, WW.index, color='.5')
    ax[1].set_xlabel('phi', color='.5')
    ax[2].plot(WW.VP*WW.RHO, WW.index, 'black')
    ax[2].set_xlabel('AI [m/s*g/cc]')
    ax[3].plot(WW.VP/WW.VS, WW.index, 'black')
    ax[3].set_xlabel('Vp/Vs')
    rlims = ax[3].get_xlim()
    for i,aa in enumerate(ax):
        if flagtops: # plot top markers on all columns
            for topn,topz in tops.iteritems():
                if (topz>=ztop) & (topz<=zbot):
                    aa.axhline(y=topz,color='blue',alpha=.8)
                    if i is 3: # plot top name also on last column
                        aa.text(x=rlims[1]+rlims[1]*.01,y=topz,s=topn,**format_tops)
    if name is not None:
        plt.suptitle(name, **format_title)
    ax[0].set_ylim(zbot,ztop)   
