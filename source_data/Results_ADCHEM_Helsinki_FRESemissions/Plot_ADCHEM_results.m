clear all
close all
clc
sim = 10;
year=17;
month=12;
day=5; % 9,14
date1=year*1000000+month*10000+day*100;


R = 8.3145;
Na = 6.022*10^23;


input_names=sprintf('species_names.dat');
fid = fopen(input_names);
spc_name_load = textscan(fid,'%s'); % MCM names
fclose(fid); 
spc_name = spc_name_load{1,1};
spc_name=strtrim(spc_name);
NSPEC=length(spc_name);

fid=fopen('MCM_RO2.txt');
RO2_name=textscan(fid,'%s');
RO2_name = RO2_name{1,1};
RO2_name=strtrim(RO2_name);
fclose(fid); 


input=sprintf('Diameter_dry.dat');
dg=load(input); % Dry geometric particle diameter for each size bin 1.07 nm to 2.5 um

input=sprintf('Height_levels.dat');
height=load(input); 
nz=length(height)

input=sprintf('dlogDp.dat');
dlogDp=load(input);


comp_prop=load('MCM33_allvoc_ELVOC_compound_prop_20180912.dat'); % column 1: molar mass (g/mol), column 2 and 3: Nannoolal coeff. (a and b resp)
Morg = comp_prop(:,1); 
A_Nannoolal = comp_prop(:,2); 
B_Nannoolal = comp_prop(:,3);
psat_org=10.^(A_Nannoolal-B_Nannoolal/298.0)*101325D0; % Vapor pressures condensable organic comp at 298 K (Pa).
Cprim=psat_org/(R*298).*Morg*1D6; % C* in ug/m^3 

fid=fopen('MCM33_allvoc_ELVOC_names_20180912.dat')
VOC_names=textscan(fid,'%s');
VOC_names = VOC_names{1,1};
VOC_names=strtrim(VOC_names);
fclose(fid); 
NCOND=length(VOC_names); % Number of condensable vapours

VBS_298K=zeros(25,7,nz); % VBS gas-phase distribution in all height levels
VBS_298K_surf=zeros(25,nz); % VBS gas-phase distribution at the surface
for jj=1:25 % 25 simulations for each date
if jj<=2    
date=date1-100+22+(jj-1);
else
date=date1+(jj-3);
end

% Load size resolved particle concentrations at all 20 height levels:
input=sprintf('pm_bins_NVOC_%d%s',date,'.dat');
pm_bins_NVOC=load(input);

input=sprintf('pm_bins_SVOC_%d%s',date,'.dat');
pm_bins_SVOC=load(input);

input=sprintf('pm_bins_SO4_%d%s',date,'.dat');
pm_bins_SO4=load(input);

input=sprintf('pm_bins_NO3_%d%s',date,'.dat');
pm_bins_NO3=load(input);

input=sprintf('pm_bins_Cl_%d%s',date,'.dat');
pm_bins_Cl=load(input);

input=sprintf('pm_bins_Na_%d%s',date,'.dat');
pm_bins_Na=load(input);

input=sprintf('pm_bins_NH4_%d%s',date,'.dat');
pm_bins_NH4=load(input);

input=sprintf('pm_bins_BC_%d%s',date,'.dat');
pm_bins_BC=load(input);

input=sprintf('pm_bins_metals_%d%s',date,'.dat');
pm_bins_metals=load(input);

input=sprintf('pn_bins_%d%s',date,'.dat');
pn_bins=load(input);

% Load gas-phase concentrations at all 20 height levels 
input=sprintf('SO2_%d%s',date,'.dat');
SO2(jj,:)=load(input);

input=sprintf('NO2_%d%s',date,'.dat');
NO2(jj,:)=load(input);

input=sprintf('NO_%d%s',date,'.dat');
NO(jj,:)=load(input);

input=sprintf('OH_%d%s',date,'.dat');
OH(jj,:)=load(input);

input=sprintf('HO2_%d%s',date,'.dat');
HO2(jj,:)=load(input);

input=sprintf('O3_%d%s',date,'.dat');
O3(jj,:)=load(input);

input=sprintf('CO_%d%s',date,'.dat');
CO(jj,:)=load(input);

input=sprintf('HNO3_%d%s',date,'.dat');
HNO3(jj,:)=load(input);

input=sprintf('H2SO4_%d%s',date,'.dat');
H2SO4(jj,:)=load(input);

input=sprintf('NH3_%d%s',date,'.dat');
NH3(jj,:)=load(input);

input=sprintf('Alkanes_%d%s',date,'.dat');
Alkanes(jj,:)=load(input);

input=sprintf('Carbonyls_%d%s',date,'.dat');
Carbonyls(jj,:)=load(input);

input=sprintf('RO2pool_%d%s',date,'.dat');
RO2pool(jj,:)=load(input);

input=sprintf('Gas_conc_%d%s',date,'.dat');
conc=load(input); % Concentration of all gases listed in spc_name 

% Create a VBS distribution of the condensable gases:
VOCs=zeros(NCOND,nz);
for j=1:NCOND
VOC_ind = find(strcmp(VOC_names(j),spc_name) == 1);
VOCs(j,:) = conc(:,VOC_ind)/Na*Morg(j)*1D12; % VOC gas-phase concentrations in ug/m^3
end
Cprim_VBS=[0.001,0.01,0.1,1,10,100,1000];
VBS_298K(jj,1,:)=sum(VOCs(Cprim<0.005,:)); % All compounds with C* <0.005 ug/m^3
VBS_298K(jj,2,:)=sum(VOCs(Cprim>=0.005 & Cprim<0.05,:)); % All compounds with 0.005 ug/m^3 < C* <0.05 ug/m^3
VBS_298K(jj,3,:)=sum(VOCs(Cprim>=0.05 & Cprim<0.5,:)); % All compounds with 0.05 ug/m^3 < C* <0.5 ug/m^3
VBS_298K(jj,4,:)=sum(VOCs(Cprim>=0.5 & Cprim<5,:)); % All compounds with 0.5 ug/m^3 < C* <5 ug/m^3
VBS_298K(jj,5,:)=sum(VOCs(Cprim>=5 & Cprim<50,:)); % All compounds with 5 ug/m^3 < C* <50 ug/m^3
VBS_298K(jj,6,:)=sum(VOCs(Cprim>=50 & Cprim<500,:)); % All compounds with 50 ug/m^3 < C* <500 ug/m^3
VBS_298K(jj,7,:)=sum(VOCs(Cprim>=500,:)); % All compounds with 500 ug/m^3 < C*

VBS_298K_surf(jj,1:7)=VBS_298K(jj,:,1); % VBS distribution at the surface


% Sum particle mass concentrations in each size bin 1.07 nm to 2.5 um:
PM2_5_NVOC(jj,:)=sum(pm_bins_NVOC,2);
PM2_5_SVOC(jj,:)=sum(pm_bins_SVOC,2);
PM2_5_SO4(jj,:)=sum(pm_bins_SO4,2);
PM2_5_NO3(jj,:)=sum(pm_bins_NO3,2);
PM2_5_NH4(jj,:)=sum(pm_bins_NH4,2);
PM2_5_BC(jj,:)=sum(pm_bins_BC,2);
PM2_5_Na(jj,:)=sum(pm_bins_Na,2);
PM2_5_Cl(jj,:)=sum(pm_bins_Cl,2);
PM2_5_metals(jj,:)=sum(pm_bins_metals,2);

PNtot(jj,:)=sum(pn_bins,2);
PNtot_3nm(jj,:)=sum(pn_bins(:,dg>3E-9),2);
PNtot_10nm(jj,:)=sum(pn_bins(:,dg>10E-9),2);

dNdlogDp_surf(jj,:)=pn_bins(1,:)/dlogDp; % Particle number size distribution at the surface (5 m)
dNdlogDp_45m(jj,:)=pn_bins(3,:)/dlogDp; % Particle number size distribution at the surface (45 m)
dNdlogDp_80m(jj,:)=pn_bins(4,:)/dlogDp; % Particle number size distribution at the surface (80 m)

pm_surf_NVOC(jj,:)=pm_bins_NVOC(1,:);
pm_surf_SVOC(jj,:)=pm_bins_SVOC(1,:);
pm_surf_SO4(jj,:)=pm_bins_SO4(1,:);
pm_surf_NO3(jj,:)=pm_bins_NO3(1,:);
pm_surf_NH4(jj,:)=pm_bins_NH4(1,:);
pm_surf_BC(jj,:)=pm_bins_BC(1,:);
pm_surf_Na(jj,:)=pm_bins_Na(1,:);
pm_surf_Cl(jj,:)=pm_bins_Cl(1,:);
pm_surf_metals(jj,:)=pm_bins_metals(1,:);

input=sprintf('HYSPLIT/HYSPLIT_Helsinki_%d%s',date,'.txt');
HYSPLIT=load(input);
Lat=HYSPLIT(1:169,6);
Long=HYSPLIT(1:169,7);
Altitude=HYSPLIT(1:169,8);

if jj==1
figure(1)    
set(gcf,'paperorientation','portrait')
set(gcf,'paperunits','centimeters');
set(gcf,'papertype','a4letter');
set(gcf,'paperposition',[0 0 20 10]);
gray2=[1:-.1:0; 1:-0.1:0; 1:-0.1:0;]';

  latlim = [40 80];
  lonlim = [-60 50];
  worldmap(latlim, lonlim)
 geoshow('landareas.shp')%, 'FaceColor', [1.0 1.0 1.0]);
 hold on
end
figure(1)
geoshow(Lat,Long,'DisplayType', 'Line', 'Color', [1-jj/25, 0, jj/25])

end


T=273.15;
p=1.01325D5;
M_air = 28.98*10^-3; %kg/mol
cM = 10^-6*p*Na./(R*T); % Air molecules cm^-3 at STP

% Gas-concentrations in ppbv at standard temperature and pressure (273.15 K, p=101325 Pa)
SO2_ppb=SO2/cM*1D9;
NO2_ppb=NO2/cM*1D9;
NO_ppb=NO/cM*1D9;
O3_ppb=O3/cM*1D9;
CO_ppb=CO/cM*1D9;
HNO3_ppb=HNO3/cM*1D9;
H2SO4_ppb=H2SO4/cM*1D9;
NH3_ppb=NH3/cM*1D9;
Alkanes_ppb=Alkanes/cM*1D9;
Carbonyls_ppb=Carbonyls/cM*1D9;

% Size resolved chemical composition at the surface:

pm_surf_tot=pm_surf_NVOC+pm_surf_SVOC+pm_surf_SO4+pm_surf_NO3+pm_surf_NH4+...
pm_surf_BC+pm_surf_Na+pm_surf_Cl+pm_surf_metals;

dmdlogDp_surf=pm_surf_tot/dlogDp;
dmdlogDp_surf_NVOC=pm_surf_NVOC/dlogDp;
dmdlogDp_surf_SVOC=pm_surf_SVOC/dlogDp;
dmdlogDp_surf_SO4=pm_surf_SO4/dlogDp;
dmdlogDp_surf_NO3=pm_surf_NO3/dlogDp;
dmdlogDp_surf_NH4=pm_surf_NH4/dlogDp;
dmdlogDp_surf_BC=pm_surf_BC/dlogDp;
dmdlogDp_surf_Na=pm_surf_Na/dlogDp;
dmdlogDp_surf_Cl=pm_surf_Cl/dlogDp;
dmdlogDp_surf_metals=pm_surf_metals/dlogDp;

% Mass fraction in each size bin:
x_surf_NVOC=pm_surf_NVOC/pm_surf_tot;
x_surf_SVOC=pm_surf_SVOC/pm_surf_tot;
x_surf_SO4=pm_surf_SO4/pm_surf_tot;
x_surf_NO3=pm_surf_NO3/pm_surf_tot;
x_surf_NH4=pm_surf_NH4/pm_surf_tot;
x_surf_BC=pm_surf_BC/pm_surf_tot;
x_surf_Na=pm_surf_Na/pm_surf_tot;
x_surf_Cl=pm_surf_Cl/pm_surf_tot;
x_surf_metals=pm_surf_metals/pm_surf_tot;



% Load SMEAR III and SMEAR II data:

GasesSMEARIII=load('smearIII_gas_data_2017.txt');
year_meas=GasesSMEARIII(:,1);
month_meas=GasesSMEARIII(:,2);
day_meas=GasesSMEARIII(:,3);
hour_meas=GasesSMEARIII(:,4);
min_meas=GasesSMEARIII(:,5);
index_meas=find(month_meas==month & day_meas==day);
time_meas=hour_meas(index_meas)/24+min_meas(index_meas)/(24*60);

NO2_meas=GasesSMEARIII(index_meas,7);
NO_meas=GasesSMEARIII(index_meas,8);
O3_meas=GasesSMEARIII(index_meas,9);
SO2_meas=GasesSMEARIII(index_meas,10);
CO_meas=GasesSMEARIII(index_meas,11);


SMPS_SMEARIII=load('smearIII_SMPS_2017.txt');
year_meas3=SMPS_SMEARIII(:,1);
month_meas3=SMPS_SMEARIII(:,2);
day_meas3=SMPS_SMEARIII(:,3);
hour_meas3=SMPS_SMEARIII(:,4);
min_meas3=SMPS_SMEARIII(:,5);
index_meas=find(month_meas3==month & day_meas3==day);
time_meas3=hour_meas3(index_meas)/24+min_meas3(index_meas)/(24*60);
dNdlogDp_SMPS=SMPS_SMEARIII(index_meas,16:end-2);
Dp_SMPS=[2.82,3.16,3.55,3.98,4.47,5.01,5.62,6.31,7.08,7.94,8.91,10,11.2,12.6,...
    14.1,15.8,17.8,20,22.4,25.1,28.2,31.6,35.5,39.8,44.7,50.1,56.2,63.1,70.8,79.4,...
    89.1,100,112,126,141,158.8,178,200,224,251,282,316,355,398,447,501,562,631,708,...
    794];
Vp_SMPS=(Dp_SMPS*1D-3).^3*pi/6; % um^3

SMPS_SMEARII=load('smearII_SMPS_2017.txt');
year_meas2=SMPS_SMEARII(:,1);
month_meas2=SMPS_SMEARII(:,2);
day_meas2=SMPS_SMEARII(:,3);
hour_meas2=SMPS_SMEARII(:,4);
min_meas2=SMPS_SMEARII(:,5);
index_meas=find(month_meas2==month & day_meas2==day);
time_meas2=hour_meas2(index_meas)/24+min_meas2(index_meas)/(24*60);
dNdlogDp_SMPS_SMEARII=SMPS_SMEARIII(index_meas,16:end-2);


dlogDp=mean(log10(Dp_SMPS(2:end))-log10(Dp_SMPS(1:end-1)));

    Nbins_SMPS=dNdlogDp_SMPS.*dlogDp;
    Nbins_SMPS_SMEARII=dNdlogDp_SMPS_SMEARII.*dlogDp;
    N=length(dNdlogDp_SMPS(:,1));
    for j=1:N
    Vbins_SMPS(j,:)=Nbins_SMPS(j,:).*Vp_SMPS;
    end
    N2=length(dNdlogDp_SMPS_SMEARII(:,1));
    for j=1:N2
    Vbins_SMPS_SMEARII(j,:)=Nbins_SMPS_SMEARII(j,:).*Vp_SMPS;
    end
    
    PN_SMPS=sum(Nbins_SMPS,2);
    PN_SMPS_SMEARII=sum(Nbins_SMPS_SMEARII,2);
    
    PV_SMPS=sum(Vbins_SMPS,2);
    PV_SMPS_SMEARII=sum(Vbins_SMPS_SMEARII,2);

% Plot some of the results:


time=[0:jj-1]/24;

figure(2)
subplot(4,3,1)
surf(time,height,SO2_ppb')
view(2)
shading flat
title('SO_2 (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,2)
surf(time,height,NO2_ppb')
view(2)
shading flat
title('NO_2 (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,3)
surf(time,height,NO_ppb')
view(2)
shading flat
title('NO (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,4)
surf(time,height,NH3_ppb')
view(2)
shading flat
title('NH_3 (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,5)
surf(time,height,O3_ppb')
view(2)
shading flat
title('O_3 (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,6)
surf(time,height,OH')
view(2)
shading flat
title('OH (molec cm^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,7)
surf(time,height,HO2')
view(2)
shading flat
title('HO_2 (molec cm^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,8)
surf(time,height,RO2pool')
view(2)
shading flat
title('\Sigma(RO_2) (molec cm^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,9)
surf(time,height,Alkanes_ppb')
view(2)
shading flat
title('\Sigma(RH) (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,10)
surf(time,height,Carbonyls_ppb')
view(2)
shading flat
title('\Sigma(RCHO) (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,11)
surf(time,height,HNO3_ppb')
view(2)
shading flat
title('HNO_3 (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(4,3,12)
surf(time,height,H2SO4')
view(2)
shading flat
title('H2SO4 (molec cm^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')



figure(3)
subplot(2,2,1)
plot(time,O3_ppb(:,[1,5]),time_meas,O3_meas,'.')
xlabel('Time (days)')
ylabel('[O_3] (ppb_v)')
legend('Model 5 m',' Model 125 m','Meas.')
subplot(2,2,2)
plot(time,NO_ppb(:,[1,5])+NO2_ppb(:,[1,5]),time_meas,NO2_meas+NO_meas,'.')
xlabel('Time (days)')
ylabel('[NO_x] (ppb_v)')
legend('Model 5 m',' Model 125 m','Meas.')
subplot(2,2,3)
plot(time,SO2_ppb(:,[1,5]),time_meas,SO2_meas,'.')
xlabel('Time (days)')
ylabel('[SO_2] (ppb_v)')
legend('Model 5 m',' Model 125 m','Meas.')
subplot(2,2,4)
plot(time,CO_ppb(:,[1,5]),time_meas,CO_meas,'.')
xlabel('Time (days)')
ylabel('[CO] (ppb_v)')
legend('Model 5 m',' Model 125 m','Meas.')

% Sum particle mass concentrations in each size bin 1.07 nm to 2.5 um:

figure(4)
subplot(3,3,1)
surf(time,height,PM2_5_NVOC')
view(2)
shading flat
title('PM_2_._5 NVOC (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,2)
surf(time,height,PM2_5_SVOC')
view(2)
shading flat
title('PM_2_._5 SVOC (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,3)
surf(time,height,PM2_5_BC')
view(2)
shading flat
title('PM_2_._5 BC (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,4)
surf(time,height,PM2_5_SO4')
view(2)
shading flat
title('PM_2_._5 SO_4 (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,5)
surf(time,height,PM2_5_NO3')
view(2)
shading flat
title('PM_2_._5 NO_3 (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,6)
surf(time,height,PM2_5_NH4')
view(2)
shading flat
title('PM_2_._5 NH_4 (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,7)
surf(time,height,PM2_5_Na')
view(2)
shading flat
title('PM_2_._5 Na (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,8)
surf(time,height,PM2_5_Cl')
view(2)
shading flat
title('PM_2_._5 Cl (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(3,3,9)
surf(time,height,PM2_5_metals')
view(2)
shading flat
title('PM_2_._5 Metals & other insoluble mat. (\mug m^-^3)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')

figure(5)
semilogx(dg,median(dmdlogDp_surf),'k--')
xlabel('Diameter (m)')
ylabel('dM/lodgDp')

figure(6)
loglog(dg,median(dmdlogDp_surf),'k.')
hold on
loglog(dg,median(dmdlogDp_surf_NVOC),'g',...
    dg,median(dmdlogDp_surf_SVOC),'g--',...
    dg,median(dmdlogDp_surf_SO4),'r',...
    dg,median(dmdlogDp_surf_NO3),'b',...
    dg,median(dmdlogDp_surf_NH4),'y',...
    dg,median(dmdlogDp_surf_BC),'k',...
    dg,median(dmdlogDp_surf_Na),'m--',...
    dg,median(dmdlogDp_surf_Cl),'c--',...
    dg,median(dmdlogDp_surf_metals),'--')
xlabel('Diameter (m)')
ylabel('dM/lodgDp (\mug m^-^3)')
legend('Tot.','NVOC','SVOC','SO_4','NO_3','NH_4','BC','Na','Cl','Metals & other mat.')
axis([1E-8 2.5E-6 1E-3 12])

figure(7)
plot(median(PNtot)*1D-6,height)
hold on
plot(median(PNtot_3nm)*1D-6,height)
hold on
plot(median(PNtot_10nm)*1D-6,height)
xlabel('PN (cm^-^3)')
ylabel('Altitude (m)')
legend('Tot.','PN(D_p>3 nm)','PN(D_p>10 nm)')


figure(8)
subplot(2,4,1)
semilogx(dg*1E9,median(dNdlogDp_surf(1:4,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(1:4,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(1:6*3,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(1:6*3,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('00:00-03:00')
xlabel('Diameter (nm)')
ylabel('dN/dlogDp (cm^-^3)')
subplot(2,4,2)
semilogx(dg*1E9,median(dNdlogDp_surf(4:7,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(4:7,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*3+1:6*6,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*3+1:6*6,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('03:00-06:00')
subplot(2,4,3)
semilogx(dg*1E9,median(dNdlogDp_surf(7:10,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(7:10,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*6+1:6*9,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*6+1:6*9,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('06:00-09:00')
subplot(2,4,4)
semilogx(dg*1E9,median(dNdlogDp_surf(10:13,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(10:13,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*9+1:6*12,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*9+1:6*12,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('09:00-12:00')
subplot(2,4,5)
semilogx(dg*1E9,median(dNdlogDp_surf(13:16,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(13:16,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*12+1:6*15,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*12+1:6*15,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('12:00-15:00')
subplot(2,4,6)
semilogx(dg*1E9,median(dNdlogDp_surf(16:19,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(16:19,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*15+1:6*18,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*15+1:6*18,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('15:00-18:00')
subplot(2,4,7)
semilogx(dg*1E9,median(dNdlogDp_surf(19:22,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(19:22,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*18+1:6*21,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*18+1:6*21,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('18:00-21:00')
subplot(2,4,8)
semilogx(dg*1E9,median(dNdlogDp_surf(22:25,:))*1D-6,'k',dg*1E9,median(dNdlogDp_45m(22:25,:))*1D-6,'k--')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*21+1:end,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*21+1:end,:)),'-o')
legend('Model 5 m','Model 45 m','Meas. SMEAR III','Meas. SMEAR II')
title('21:00-24:00')

figure(9)
subplot(4,1,1)
h=gca
surf(time,dg,dNdlogDp_surf'*1D-6)
view(2)
shading interp
colorbar
caxis([0 2E4])
colormap jet
xlabel('Time (days)')
ylabel('Diameter (m)')
title('Modelled 5 m')
axis([0 1 1E-9 1E-6]) 
set(h,'yscale','log')
subplot(4,1,2)
h=gca
surf(time,dg,dNdlogDp_45m'*1D-6)
view(2)
shading interp
colorbar
caxis([0 2E4])
colormap jet
xlabel('Time (days)')
ylabel('Diameter (m)')
title('Modelled 45 m')
axis([0 1 1E-9 1E-6]) 
set(h,'yscale','log')
subplot(4,1,3)
surf(time_meas3,Dp_SMPS*1D-9,dNdlogDp_SMPS')
h=gca
view(2)
shading interp
colorbar
caxis([0 2E4])
colormap jet
xlabel('Time (days)')
ylabel('Diameter (m)')
title('Meas SMEAR III')
axis([0 1 1E-9 1E-6])
set(h,'yscale','log')
subplot(4,1,4)
surf(time_meas2,Dp_SMPS*1D-9,dNdlogDp_SMPS_SMEARII')
h=gca
view(2)
shading interp
colorbar
caxis([0 2E4])
colormap jet
xlabel('Time (days)')
ylabel('Diameter (m)')
title('Meas SMEAR II')
axis([0 1 1E-9 1E-6])
set(h,'yscale','log')

figure(10)
h=gca
bar([0:24],VBS_298K_surf,'stacked')
legend('log_1_0(C*[\mug/m^3])=0.001','log_1_0(C*[\mug/m^3])=0.01','log_1_0(C*[\mug/m^3])=0.1','log_1_0(C*[\mug/m^3])=1.0','log_1_0(C*[\mug/m^3])=10.0','log_1_0(C*[\mug/m^3])=100','log_1_0(C*[\mug/m^3])=1000')
ylabel('\mug/m^3')
xlabel('All 25 cases from 00:00 to 24:00')
set(h,'yscale','log')
