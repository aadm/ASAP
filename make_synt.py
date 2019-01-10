import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import bruges as b
import hackathon2018_library as hl
%load_ext autoreload
%autoreload 2

%matplotlib qt

# import scipy.io as sio
# w2_new 
# temp = sio.loadmat('155_5.mat')
# pd datframe from dict (... orient = index)


w1=pd.read_csv('qsiwell1.csv', index_col=0)
w2=pd.read_csv('qsiwell2.csv', index_col=0)
w3=pd.read_csv('qsiwell3.csv', index_col=0)
w4=pd.read_csv('qsiwell4.csv', index_col=0)
w5=pd.read_csv('qsiwell5.csv', index_col=0)


wells=    [ w2,       w3,       w5]
names=     ['15/5-5','15/5-3', '15/5-6']
tops_w2={'Lista': 2127, 'Heimdal': 2154,'OWC': 2187}
tops_w3={'Lista': 2153, 'Heimdal': 2180}
tops_w5={'Lista': 2145, 'Heimdal': 2172,'OWC': 2185}
tops=pd.DataFrame()
for i,val in enumerate([tops_w2,tops_w3,tops_w5]):
    tempdf=pd.DataFrame.from_dict(val, orient='index')
    tempdf.columns=[names[i]]
    tops=pd.concat([tops, tempdf],axis=1)

for i,p in enumerate(wells):
    z1 = tops[names[i]]['Lista']-100
    z2 = tops[names[i]]['Heimdal']+100
    hl.quicklook(p,ztop=z1,zbot=z2,name=names[i],tops=tops[names[i]])
    plt.savefig('FIGS/{}_original_wells.png'.format(names[i].replace('/','_')))

# forward modeling parameters
# sample rate: 4 ms
# wavelet: Ricker 30 Hz
dt=.004
wavl=b.filters.ricker(.25, dt, 30)

# angle range
ang=np.linspace(5,30,2)

# calculate time-depth tables
tdr_w2=hl.td(w2,KB=26,WD=109)
topres_w2=hl.get_twt(tdr_w2,tops['15/5-5']['Heimdal'])

tdr_w3=hl.td(w3,KB=26,WD=109)
topres_w3=hl.get_twt(tdr_w2,tops['15/5-3']['Heimdal'])

tdr_w5=hl.td(w5,KB=26,WD=109)
topres_w5=hl.get_twt(tdr_w5,tops['15/5-6']['Heimdal'])


kkk, tdr, topres = 'W155_5', tdr_w2, topres_w2
kkk, tdr, topres = 'W155_X', tdr_w2, topres_w2
kkk, tdr, topres = 'W155_3', tdr_w3, topres_w3
kkk, tdr, topres = 'W155_6', tdr_w5, topres_w5


well_files = get_well_files(wells_dir='SYNTWELLS', name=kkk)
t1=topres-0.048
t2=topres+0.208
print(t1,t2,t2-t1)
# then, output everything to a pandas dataframe
megasynt = pd.DataFrame()
for ff in well_files:
    nn=pd.read_csv(ff, usecols=[1,2,3,4], index_col=0)
    nn.columns=['VP','VS','RHO']
    nnt,_,_ = hl.welltime(nn,tdr,dt=dt,qcplot=False)
    rc,synt=hl.make_synt(nnt,ang,wavl,method='aki')
    outdf = nnt[(nnt.index>=t1) & (nnt.index<=t2)].copy()
    it1=np.abs(nnt.index-t1).argmin()
    it2=np.abs(nnt.index-t2).argmin()
    synt_cut = synt[it1:it2,:]
    outdf['NEAR'] = synt_cut[:,0]
    outdf['FAR'] = synt_cut[:,1]
    outdf['WELL'] = ff[10:-16] # will pick up the well name
    outdf['ID'] = ff[-15:] 
    megasynt = megasynt.append(outdf)
    # hl.plot_synt(nnt,synt,t1,t2,20)

megasynt['WELL']=megasynt['WELL'].astype('category')
megasynt['ID']=megasynt['ID'].astype('category')
megasynt.to_csv('syntdataset_{}.csv'.format(kkk))




#~~~~~ figures for Lucy

dt=.001  # just to get nicer plots with smaller sample rate
wavl=b.filters.ricker(.25, dt, 30)

# angle range
ang=np.linspace(5,30,2)

files=[
'W155_5_Z60_Sw000_Por35',
'W155_5_Z60_Sw025_Por35',
'W155_5_Z60_Sw075_Por35',
'W155_5_Z60_Sw100_Por35']

for ff in files:
    wellname = '15/5-5 Glytne'
    nn=pd.read_csv(ff, usecols=[1,2,3,4], index_col=0)
    nn.columns=['VP','VS','RHO']
    nnt,_,_ = hl.welltime(nn,tdr_w2,dt=dt,tops=tops_w2,qcplot=False)
    rc,synt=hl.make_synt(nnt,ang,wavl,method='aki')
    hl.plot_synt(nnt,synt,t1,t2,20)
    plt.savefig('FIGS/{}_synth_variable_Sw.png'.format(names[i].replace('/','_')))




syntw2=pd.read_csv('syntdataset_W155_5.csv')
syntw2soft=pd.read_csv('syntdataset_W155_X.csv')

uuu='15_5-5'
nears, fars = get_nears_fars(syntw2)

uuu='15_5-5_soft'
nears, fars = get_nears_fars(syntw2soft)


clip=abs(np.percentile([nears,fars], 0.80))

f,ax = plt.subplots(nrows=2,ncols=1,figsize=(12,10))
im0=ax[0].imshow(nears,cmap='Greys',vmax=clip,vmin=-clip,aspect='auto', interpolation='bicubic')
im1=ax[1].imshow(fars,cmap='Greys',vmax=clip,vmin=-clip,aspect='auto', interpolation='bicubic')
cax = f.add_axes([0.2, 0.05, 0.6, 0.025])
f.colorbar(im0, cax=cax, orientation='horizontal')
ax[0].set_title('15/5-5 Near')
ax[1].set_title('15/5-5 Far')
ax[1].set_xlabel('Variation of petrophysical properties')
plt.savefig('FIGS/{}_synth_variable_Sw.png'.format(uuu))


test_samples = np.load("test_out_49.npy")

opt={'linewidth':4, 'color':'black'}

fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(8, 6))
# real
ax[0].plot(test_samples[0,5,0,:], range(64), **opt) # near
ax[1].plot(test_samples[0,5,1,:], range(64), **opt) # far
# reconstructed
ax[2].plot(test_samples[1,5,0,:], range(64), **opt) # near
ax[3].plot(test_samples[1,5,1,:], range(64), **opt) # far
for aa in [ax[0],ax[2]]:
    aa.set_xlabel('Near')
for aa in [ax[1],ax[3]]:
    aa.set_xlabel('Far')
for aa in ax[:3]:
    aa.set_title('ORIGINAL')
for aa in ax[2:]:
    aa.set_title('RECONSTRUCTED')
for aa in ax:
    aa.set_yticklabels([])
    aa.set_xticklabels([])
    aa.invert_yaxis()
    aa.set_xlim(-3,3)
plt.savefig('FIGS/original-reconstructed.png')


