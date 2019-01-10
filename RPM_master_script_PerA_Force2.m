%%% RPM model: 1) Contact cement model; 2) Constant cement model; 3)Increasing cement model
%%%----------------------------------------------------------------------------------------
%% 1) Contact cement modeling ****

%% input parameters
Ft=input('Reduced shear effect: 0:no friction, 1:no slip (Typical value at 2km burial depth: 0.5)')

Koil=1.1e9;
Kbrine=2.75e9;
Rhooil=0.8;
Rhobrine=1.09;
Sw=input('Brine saturation?');
So=1-Sw;
Kfl=1./(So/Koil+Sw/Kbrine);
Rhofl=So*Rhooil+Sw*Rhobrine;

Clay=input('What is the fraction of clay?');
Qz=1-Clay;
Kmin=((1./(Qz/36.8e9+Clay/17.5e9))+(Qz*36.8e9+Clay*17.5e9))./2;
Gmin=((1./(Qz/44e9+Clay/7.5e9))+(Qz*44e9+Clay*7.5e9))./2;
rhomin=2.65;

Kcem=36.8e9;
Gcem=44e9;
rhoc=2.65;

por_max=input('Max porosity in fraction? (Critical porosity)')
Vcem=input('Volume of cement? (fraction)')
por_hs=por_max-Vcem;

%% Estimation of alpha and stiffnesses in Dvorkin-Nur contact cement model
por=[0.001:0.001:por_max];

for i=1:length(por)


alpha(i)=(2*(por_max-por(i))/(3*(1-por_max)) )^0.5;

K=Kmin./1e9;
G=Gmin./1e9;
Kc=Kcem./1e9;
Gc=Gcem./1e9;
Poiss=(3*K-2*G)./(2*(3*K+G));
PoissC=(3*Kc-2*Gc)./(2*(3*Kc+Gc));

AAn=(2*Gc/(3.14*G))* ((1-Poiss)*(1-PoissC)/(1-2*PoissC));
AAt=Gc/(3.14*G);
Ct=10^(-4)*(9.654*Poiss^2+4.945*Poiss+3.1)*AAt^(0.01867*(Poiss^2)+0.4011*Poiss-1.8186);
Bt=(0.0573*Poiss^2 + 0.0937*Poiss +0.202)*AAt^(0.0274*Poiss^2 + 0.0529*Poiss - 0.8765);
At=-10^(-2)*(2.26*Poiss^2+2.07*Poiss+2.3) * AAt^(0.079*Poiss^2 + 0.1754*Poiss - 1.342); 
St(i)=Ft.*(At*alpha(i)^2 + Bt*alpha(i) + Ct);
Cn=0.00024649*AAn^(-1.9846);
Bn=0.20405*AAn^(-0.89008);
An=-0.024153*AAn^(-1.3646);
Sn(i)=An*alpha(i)^2+Bn*alpha(i)+Cn;


%% Calculate dry rock properties (Contact theory)

n=20-34*por_max+14*por_max.^2;

Vpc=sqrt((Kc+4/3*Gc)./rhoc);
Vsc=sqrt(Gc./rhoc);
Mc=rhoc*Vpc^2;
KeffC(i)=(1/6)*n*(1-por_max)*Mc*Sn(i);
GeffC(i)=(3/5)*KeffC(i)+(3/20)*n*(1-por_max)*Gc*St(i);

rhodC(i)=rhomin*(1-por(i));
Vpd(i)=sqrt( (KeffC(i)+ (4/3)*GeffC(i) )/rhodC(i));
Vsd(i)=sqrt(  GeffC(i)/rhodC(i) );

%% Calculate saturated rock properties (Gassmann theory)

FactorB(i)=((KeffC(i)).*1e9)./(Kmin-(KeffC(i).*1e9)) + Kfl./(por(i)*(Kmin-Kfl));

KsatC(i)=Kmin*FactorB(i)/(1+FactorB(i));
GsatC(i)=GeffC(i).*1e9;
rhosC(i)=rhomin*(1-por(i))+Rhofl*por(i);

VpsatC(i)=0.001*sqrt((KsatC(i) + (4/3)*GsatC(i) )./(rhosC(i).*1000) );
VssatC(i)=0.001*sqrt(GsatC(i)./(rhosC(i).*1000));

end
b_p=por_max*1000;
a_p=por_hs*1000;
K_output=KeffC(round(a_p))*1e9;
G_output=GeffC(round(a_p))*1e9;
figure
plot(por(round(a_p):1:round(b_p)),VpsatC(round(a_p):1:round(b_p)),'b--')
hold on
xlabel('Porosity'),ylabel('Vp (km/s)')

%%%--------------------------------------------------------------------
%% 2) Constant cement model (Modified Hashin-Shtrikman lower bound)****

Keff=[];
Geff=[];
Khs=K_output;
Ghs=G_output;

porvec=[0:0.001:por_hs];

for i=1:length(porvec)

%% Calculate dry rock properties

    Ktemp1_CC(i)=   (porvec(i)/por_hs)/(Khs+(4/3)*Ghs);
    Ktemp2_CC(i)= ( 1-(porvec(i)/por_hs) ) / (Kmin+(4/3)*Ghs);
    Keff_CC(i)= ( Ktemp1_CC(i) + Ktemp2_CC (i) )^(-1)  - (4/3)*Ghs;
    
    Gtemp1_CC(i) = ( porvec(i)/por_hs) / ( Ghs + (Ghs/6)*( (9*Khs+8*Ghs)/(Khs+2*Ghs)) );
    Gtemp2_CC(i) = (1-(porvec(i)/por_hs) ) / ( Gmin +(Ghs/6)* (9*Khs+8*Ghs)/(Khs+2*Ghs) );
    Gtemp3_CC(i) = (Ghs/6)*( (9*Khs+8*Ghs)/(Khs+2*Ghs));
    Geff_CC(i)= ( Gtemp1_CC(i) + Gtemp2_CC(i) )^(-1) - Gtemp3_CC(i);
    
    rhod(i)=2650*(1-porvec(i));
    Vpd_CC(i)=sqrt( (Keff_CC(i)+ (4/3)*Geff_CC(i) )/rhod(i));
    Vsd_CC(i)=sqrt(  Geff_CC(i)/rhod(i) );
    
%% Calculate saturated rock properties
    
    FactorB_CC(i)=Keff_CC(i)/(Kmin-Keff_CC(i)) + Kfl/(porvec(i)*(Kmin-Kfl));
    Ksat_CC(i)=Kmin*FactorB_CC(i)/(1+FactorB_CC(i));
    Gsat_CC(i)=Geff_CC(i);
    rhos(i)=2.65e3*(1-porvec(i))+Rhofl*1000*porvec(i);
    Vpsat_CC(i)=0.001*sqrt((Ksat_CC(i) + (4/3)*Gsat_CC(i) )/rhos(i) );
    Vssat_CC(i)=0.001*sqrt((Gsat_CC(i))./rhos(i));
    Vpsat_CC(1)=0.001*sqrt((Kmin+4/3*Gmin)./2650);
end

plot (porvec,Vpsat_CC,'b')

%% Regression of high-porosity "leg" of constant cement model
P_RC=polyfit(porvec(250:(round(por_hs*1000)))',Vpsat_CC(250:(round(por_hs*1000)))',1);
Vp_reg=P_RC(2)+P_RC(1)*porvec;

S_RC=polyfit(porvec(250:(round(por_hs*1000)))',Vssat_CC(250:(round(por_hs*1000)))',1);
Vs_reg=S_RC(2)+S_RC(1)*porvec;

plot (porvec,Vp_reg,'r+')

%%%--------------------------------------------------------------------
%%  3) Increasing cement model (Modified Hashin-Shtrikman upper bound) ****


for i=1:length(porvec)

Ktemp1_IC(i)=( 1-(porvec(i)/por_hs) ) / (Kmin+(4/3)*Gmin);

Ktemp2_IC(i)=(porvec(i)/por_hs)/(Khs+(4/3)*Gmin);

Keff_IC(i)= ( Ktemp1_IC(i) + Ktemp2_IC (i) )^(-1)  - (4/3)*Gmin;

Gtemp1_IC(i) = ( 1-(porvec(i)/por_hs)) / ( Gmin + (Gmin/6)*( (9*Kmin+8*Gmin)/(Kmin+2*Gmin)) );

Gtemp2_IC(i) = (porvec(i)/por_hs) / ( Ghs +(Gmin/6)* (9*Kmin+8*Gmin)/(Kmin+2*Gmin) );

Gtemp3_IC(i) = (Gmin/6)*( (9*Kmin+8*Gmin)/(Kmin+2*Gmin));

Geff_IC(i)= ( Gtemp1_IC(i) + Gtemp2_IC(i) )^(-1) - Gtemp3_IC(i);

rhod3(i)=2650*(1-porvec(i));
Vpd3_IC(i)=sqrt( (Keff_IC(i)+ (4/3)*Geff_IC(i) )/rhod(i));
Vsd3_IC(i)=sqrt(  Geff_IC(i)/rhod(i) );

FactorB_IC(i)=Keff_IC(i)/(Kmin-Keff_IC(i)) + Kfl/(porvec(i)*(Kmin-Kfl));
Ksat_IC(i)=Kmin*FactorB_IC(i)/(1+FactorB_IC(i));
Gsat_IC(i)=Geff_IC(i);

rhos(i)=2.65e3*(1-porvec(i))+Rhofl*1000*porvec(i);

Vpsat_IC(i)=0.001*sqrt((Ksat_IC(i) + (4/3)*Gsat_IC(i) )/rhos(i) );
Vssat_IC(i)=0.001*sqrt((Gsat_IC(i))./rhos(i));
Vpsat_IC(1)=0.001*sqrt((Kmin+4/3*Gmin)./2650);
end
plot (porvec,Vpsat_IC,'k')

%% load data---------------------

Vp_data=input('Vp data vector (km/s):');
Por_data=input('Porosity data vector (fraction):');

plot(Por_data, Vp_data, 'co')


%%%--------------------------------------------------------
%% ------------create pdf of data-----------

%trdat1=input('vector 1?');
%trdat2=input('vector 2?');
trdat1=Vp_data;
trdat2=Por_data;
%nbin1=input('number of bins for vector 1?');
%nbin2=input('number of bins for vector 2?');
nbin1=20;
nbin2=20;
%smfilt1=input('smoothing filter size? ');
%smfilt2=input('smoothing filter std? ');
smfilt1=5;
smfilt2=3;
smfilt=[smfilt1 smfilt2];
trdat=[trdat1 trdat2];

%% remove NaNs
trdat(isnan(trdat(:,1)),:)=[]; trdat(isnan(trdat(:,2)),:)=[]; 

%% select good default bin size and filter size
optbin=floor(sqrt(length(trdat)/0.03));
if nargin==1, nbin1=optbin; nbin2=optbin;end;
if prod(size(nbin1))==1, nbin1=linspace(min(trdat(:,1)),max(trdat(:,1)),nbin1); end;  
if prod(size(nbin2))==1, nbin2=linspace(min(trdat(:,2)),max(trdat(:,2)),nbin2); end;  
n1=length(nbin1); n2=length(nbin2); ndat=length(trdat(:,1));

%% raw histogram
pdftbl=zeros(n1,n2);  
for k=1:ndat
 i1=sum(trdat(k,1)>=nbin1); i2=sum(trdat(k,2)>=nbin2); 
pdftbl(i1,i2)=pdftbl(i1,i2)+1;
end;

%% smoothed pdfs
pdftbl=filter2(smfilt,pdftbl);
cpdftbl=pdftbl;
cpdftbl=cpdftbl./sum(sum(cpdftbl));

%% plotting pdfs/contours
%imagesc(nbin2,nbin1,cpdftbl); axis xy; colormap hot;colorbar;
 contour(nbin2,nbin1,cpdftbl,15,'-k'); 
%plot(trdat(:,2),trdat(:,1),'.c','markersize',1); hold off;
