clear all
close all
clc

cd('/home/monakurp/makelankatu_simulations/Results_ADCHEM_Helsinki_FRESemissions/')

%%

year=17;
month= 6;
day=9;
date1=year*1000000+month*10000+day*100;

if month==6 && day==9
  T = 281.15;
elseif month==6 && day==14
  T = 287.0;
elseif month==12 && day==5
  T = 273.15;  
elseif month==12 && day==7
  T = 270.15;  
end
p=1.01325D5;
Na = 6.022*10^23;
R = 8.3145;
M_air = 28.98*10^-3; %kg/mol
cM = 10^-6*p*Na./(R*T); % Air molecules cm^-3 at STP

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
nz=length(height);

input=sprintf('dlogDp.dat');
dlogDp=load(input);

% VBS:

fid=fopen('MCM33_allvoc_ELVOC_names_20180912.dat');
VOC_names=textscan(fid,'%s');
VOC_names = VOC_names{1,1};
VOC_names=strtrim(VOC_names);
fclose(fid); 
NCOND=length(VOC_names); % Number of condensable vapours


comp_prop=load('MCM33_allvoc_ELVOC_compound_prop_20180912.dat'); % column 1: molar mass (g/mol), column 2 and 3: Nannoolal coeff. (a and b resp)
Morg = comp_prop(:,1); 
A_Nannoolal = comp_prop(:,2); 
B_Nannoolal = comp_prop(:,3);
psat_org=10.^(A_Nannoolal-B_Nannoolal/T )*p; % Vapor pressures condensable organic comp (Pa)
Cprim=psat_org/(R * T ).* Morg * 1D6; % C* in ug/m^3 

% Create a VBS distribution of the condensable gases:
VOCs=zeros(NCOND,nz);
for i=1:NCOND
  VOC_ind(i)=find(strcmp(VOC_names(i),spc_name)==1); % Index of condensabel organic vapours
end


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

  dNdlogDp_surf(jj,:)=pn_bins(3,:)/dlogDp; % Particle number size distribution at the surface (45 m)
  % Particle number size distributions (dN/dlogDp) at different altitudes:
  dNdlogDp_5m(jj,:)=pn_bins(1,:)'./dlogDp;
  dNdlogDp_20m(jj,:)=pn_bins(2,:)'./dlogDp;
  dNdlogDp_45m(jj,:)=pn_bins(3,:)'./dlogDp;
  dNdlogDp_80m(jj,:)=pn_bins(4,:)'./dlogDp;
  dNdlogDp_125m(jj,:)=pn_bins(5,:)'./dlogDp;
  dNdlogDp_180m(jj,:)=pn_bins(6,:)'./dlogDp;
  dNdlogDp_245m(jj,:)=pn_bins(7,:)'./dlogDp;
  dNdlogDp_320m(jj,:)=pn_bins(8,:)'./dlogDp;
  dNdlogDp_405m(jj,:)=pn_bins(9,:)'./dlogDp;
  dNdlogDp_500m(jj,:)=pn_bins(10,:)'./dlogDp;  

  pm_surf_NVOC(jj,:)=pm_bins_NVOC(1,:);
  pm_surf_SVOC(jj,:)=pm_bins_SVOC(1,:);
  pm_surf_SO4(jj,:)=pm_bins_SO4(1,:);
  pm_surf_NO3(jj,:)=pm_bins_NO3(1,:);
  pm_surf_NH4(jj,:)=pm_bins_NH4(1,:);
  pm_surf_BC(jj,:)=pm_bins_BC(1,:);
  pm_surf_Na(jj,:)=pm_bins_Na(1,:);
  pm_surf_Cl(jj,:)=pm_bins_Cl(1,:);
  pm_surf_metals(jj,:)=pm_bins_metals(1,:);
  
  for j = 1:NCOND
    VOCs(j,:) = conc(:,VOC_ind(j))/Na*Morg(j)*1D12; % VOC gas-phase concentrations in ug/m^3
  end

  % Group condensable vapour concentrations into VBS bins:
  C_VBS=[1D-3,1D-2,1D-1,1D0,1D1,1D2,1D3,1D4]; % Saturation concentration bins VBS (ug/m^3)
  VBS_1(jj,:)=sum(VOCs(Cprim<=5D-3,:));
  VBS_2(jj,:)=sum(VOCs(Cprim>5D-3 & Cprim<=5D-2,:));
  VBS_3(jj,:)=sum(VOCs(Cprim>5D-2 & Cprim<=5D-1,:));
  VBS_4(jj,:)=sum(VOCs(Cprim>5D-1 & Cprim<=5D0, :));
  VBS_5(jj,:)=sum(VOCs(Cprim>5D0  & Cprim<=5D1, :));
  VBS_6(jj,:)=sum(VOCs(Cprim>5D1  & Cprim<=5D2, :));
  VBS_7(jj,:)=sum(VOCs(Cprim>5D2  & Cprim<=5D3, :));
  VBS_8(jj,:)=sum(VOCs(Cprim>5D3  & Cprim<=5D4, :));
  
  % Separate OCNV and OCSV and change the units bac to #/cm3
  OCNV(jj,:) = sum( VOCs(Cprim<1,:) ); % C* < 1 ug/m^3
  is_OCNV = Cprim<1;
  sep_OCNV = zeros(length(is_OCNV),nz);
  i = 1;
  for j = 1:NCOND
    if is_OCNV(j)
      sep_OCNV(i,:) = VOCs(j,:) * Na / Morg(j) * 1e-12;
      i = i+1;
    end
  end
  OCNV(jj,:) = sum( sep_OCNV );
  
  OCSV(jj,:) = sum( VOCs(Cprim>=1 & Cprim<=1000,:) ); % 1 < C* < 1000 ug/m^3
  is_OCSV = Cprim>=1 & Cprim<=1000;
  sep_OCSV = zeros(length(is_OCSV),nz);
  i = 1;
  for j = 1:NCOND
    if is_OCSV(j)
      sep_OCSV(i,:) = VOCs(j,:) * Na / Morg(j) * 1e-12;
      i = i+1;
    end
  end
  OCSV(jj,:) = sum( sep_OCSV );
  
  
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
    lonlim = [-60 40];
    worldmap(latlim, lonlim)
    geoshow('landareas.shp')%, 'FaceColor', [1.0 1.0 1.0]);
    hold on
  end
  figure(1)
  geoshow(Lat,Long,'DisplayType', 'Line', 'Color', [1-jj/25, 0, jj/25])
end
cmap = [linspace( 1, 0, 25 ); zeros(1,25); linspace( 1/25, 1, 25 )];
colormap( cmap' )
c = colorbar( gca, 'XTickLabel',...
  {'0:00','3:00','6:00','9:00','12:00','15:00','18:00','21:00','24:00'}, ...
  'XTick', (0:3:24)/24 );
c.Label.String = 'Time UTC (HH:MM)';

filepath = '/home/monakurp/Manuscripts/Kurppa_et_al_makelankatu_evaluation/git_repo/fig';
filename = sprintf('%s/adchem_trajectory_%02d%02d.pdf',filepath,month,day);
set(gcf,'Units','inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,filename,'-dpdf','-r0')

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
OCNV_ppb = OCNV/cM*1D9;
OCSV_ppb = OCSV/cM*1D9;

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


VBS_surf(:,1)=VBS_1(:,1);
VBS_surf(:,2)=VBS_2(:,1);
VBS_surf(:,3)=VBS_3(:,1);
VBS_surf(:,4)=VBS_4(:,1);
VBS_surf(:,5)=VBS_5(:,1);
VBS_surf(:,6)=VBS_6(:,1);
VBS_surf(:,7)=VBS_7(:,1);
VBS_surf(:,8)=VBS_8(:,1);

%% Plot some of the results:


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
plot(time,NO_ppb(:,[1,5])+NO2(:,[1,5]),time_meas,NO2_meas+NO_meas,'.')
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
semilogx(dg*1E9,median(dNdlogDp_surf(1:4,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(1:6*3,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(1:6*3,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('00:00-03:00')
xlabel('Diameter (nm)')
ylabel('dN/dlogDp (cm^-^3)')
subplot(2,4,2)
semilogx(dg*1E9,median(dNdlogDp_surf(4:7,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*3+1:6*6,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*3+1:6*6,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('03:00-06:00')
subplot(2,4,3)
semilogx(dg*1E9,median(dNdlogDp_surf(7:10,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*6+1:6*9,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*6+1:6*9,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('06:00-09:00')
subplot(2,4,4)
semilogx(dg*1E9,median(dNdlogDp_surf(10:13,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*9+1:6*12,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*9+1:6*12,:)),'-o')
legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('09:00-12:00')
subplot(2,4,5)
semilogx(dg*1E9,median(dNdlogDp_surf(13:16,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*12+1:6*15,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*12+1:6*15,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('12:00-15:00')
subplot(2,4,6)
semilogx(dg*1E9,median(dNdlogDp_surf(16:19,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*15+1:6*18,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*15+1:6*18,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('15:00-18:00')
subplot(2,4,7)
semilogx(dg*1E9,median(dNdlogDp_surf(19:22,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*18+1:6*21,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*18+1:6*21,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('18:00-21:00')
subplot(2,4,8)
semilogx(dg*1E9,median(dNdlogDp_surf(22:25,:))*1D-6,'k')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS(6*21+1:end,:)),'-+')
hold on
semilogx(Dp_SMPS,median(dNdlogDp_SMPS_SMEARII(6*21+1:end,:)),'-o')
%legend('Model','Meas. SMEAR III','Meas. SMEAR II')
title('21:00-24:00')

filepath = '/home/monakurp/Manuscripts/Kurppa_et_al_makelankatu_evaluation/git_repo/fig';
filename = sprintf('%s/adchem_psd_%02d%02d.pdf',filepath,month,day);
set(gcf, 'Position', get(0, 'Screensize'));
set(gcf,'PaperUnits','normalized');
set(gcf,'PaperPosition', [0 0 1 1]);
set(gcf,'PaperOrientation','landscape');
print(gcf,filename,'-dpdf','-r0')

figure(9)
subplot(3,1,1)
h=gca;
surf(time,dg,dNdlogDp_surf'*1D-6)
view(2)
shading interp
colorbar
caxis([0 2E4])
colormap jet
xlabel('Time (days)')
ylabel('Diameter (m)')
title('Modelled')
axis([0 1 1E-9 1E-6]) 
set(h,'yscale','log')
subplot(3,1,2)
surf(time_meas3,Dp_SMPS*1D-9,dNdlogDp_SMPS')
h=gca;
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
subplot(3,1,3)
surf(time_meas2,Dp_SMPS*1D-9,dNdlogDp_SMPS_SMEARII')
h=gca;
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
subplot(2,1,1)
surf(time,height,OCNV_ppb')
view(2)
shading flat
title('OCNV (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')
subplot(2,1,2)
surf(time,height,OCSV_ppb')
view(2)
shading flat
title('OCSV (ppb_v)')
colorbar
xlabel('Time (days)')
ylabel('Altitude (m)')

figure(11)
semilogx( dg*1E9, dNdlogDp_20m )

figure(12)
h=gca;
bar([0:24],VBS_surf,'stacked')
legend('log_1_0(C*[\mug/m^3])=0.001','log_1_0(C*[\mug/m^3])=0.01','log_1_0(C*[\mug/m^3])=0.1','log_1_0(C*[\mug/m^3])=1.0','log_1_0(C*[\mug/m^3])=10.0','log_1_0(C*[\mug/m^3])=100','log_1_0(C*[\mug/m^3])=1000')
ylabel('\mug/m^3')
xlabel('All 25 cases from 00:00 to 24:00')
set(h,'yscale','log')

%% Input data for PALM
% mechanism: salsa + simple

% ------------------------------------------------------------------------%
% FILL: from here...

date    = sprintf('%02d%02d20%02d', day, month, year);
fileout = sprintf('/home/monakurp/makelankatu_simulations/source_data/ADCHEM_data/%s_FRES.mat',...
                  date );
datamat = matfile( fileout, 'Writable',true );
tstart = 0.5;
tend   = 23.5;
dt     = 1;
zi     = height<501; % m

% nbin = 2, 8; reglim = 2.5e-9, 10.0e-9, 2.5e-6 
% dmid  = [3.535535e-09, 7.071069e-09, 1.412119e-08, 2.815878e-08,...
%          5.615085e-08, 1.119693e-07, 2.232756e-07, 4.452294e-07,...
%          8.878229e-07, 1.770389e-06];
% dlims = [2.50000050e-09, 5.00000096e-09, 1.00000018e-08, 1.99408000e-08,...
%          3.97635432e-08, 7.92916717e-08, 1.58113908e-07, 3.15291723e-07,...
%          6.28716805e-07, 1.25371138e-06, 2.50000000e-06];

% nbin = 2, 8; reglim = 2.5e-9, 10.0e-9, 1.0e-6 
% dmid  = [3.53553459e-09, 7.07106914e-09, 1.33352167e-08, 2.37137412e-08,...
%          4.21696575e-08, 7.49894332e-08, 1.33352164e-07, 2.37137407e-07,...
%          4.21696565e-07, 7.49894315e-07];
% dlims = [2.50000050e-09, 5.00000096e-09, 1.00000018e-08, 1.77827973e-08,...
%          3.16227821e-08, 5.62341419e-08, 1.00000016e-07, 1.77827969e-07,...
%          3.16227813e-07, 5.62341406e-07, 1.00000000e-06]; 

% nbin = 3, 7; reglim = 2.5e-9, 15.0e-9, 1.0e-6 
dmid  = [3.37001604e-09, 6.12372552e-09, 1.11275477e-08, 2.02474592e-08,...
         3.68917214e-08, 6.72182667e-08, 1.22474507e-07, 2.23153698e-07,...
         4.06595415e-07, 7.40833931e-07];
dlims = [2.50000050e-09, 4.54280235e-09, 8.25481966e-09, 1.50000027e-08,...
         2.73306353e-08, 4.97975659e-08, 9.07332577e-08, 1.65319808e-07,...
         3.01219638e-07, 5.48834838e-07, 1.00000000e-06];        
% ... until here

% ------------------------------------------------------------------------%
% time
datamat.year   = year+2000;
datamat.smonth = month;
datamat.sday   = day;
datamat.shour  = 0;
datamat.emonth = month;
datamat.eday   = day;
datamat.ehour  = 23;
plusUTC        = 3;

si = 1;
ei = 1;


% ------------------------------------------------------------------------%
% gases:
gas_names        = {'NO','NO2','O3','OH','RH',     'RO2',    'RCHO',     'HO2','H2SO4','HNO3','NH3','OCNV','OCSV'};
gas_names_ADCHEM = {'NO','NO2','O3','OH','Alkanes','RO2pool','Carbonyls','HO2','H2SO4','HNO3','NH3','OCNV','OCSV'};

% aerosol mass fractions
mass_frac_names = {'SO4','OC','BC','NO3','NH4'};
pmsum = ( PM2_5_SO4 + PM2_5_SVOC + PM2_5_NVOC + PM2_5_BC + PM2_5_NO3 + ...
          PM2_5_NH4 + PM2_5_metals );

% ------------------------------------------------------------------------%

time   = tstart:dt:tend;

% ------------------------------------------------------------------------%

% Save height levels (m):
datamat.z = height(zi);

% Save gas concentrations:
for i = 1:length( gas_names )
  conc = eval( gas_names_ADCHEM{i} ) / cM * 1e6;
  datamat.( gas_names{i} ) = conc(si:end-ei,zi);
end

% Save aerosol mass fractions:
mass_fracs = zeros( length( time ), sum( zi ), length( mass_frac_names ) );
for i = 1:length( mass_frac_names )
  if strcmp( mass_frac_names{i}, 'OC' ) == 1
    conc = PM2_5_SVOC + PM2_5_NVOC + PM2_5_metals;
  else
    conc = eval( sprintf( 'PM2_5_%s', mass_frac_names{i} ) );
  end
  mass_fracs(:,:,i) =  conc(si:end-ei,zi) ./ pmsum(si:end-ei,zi);
end
datamat.mass_fracs = mass_fracs;
datamat.composition_name = mass_frac_names;

% Save aerosol number concentrations:
% First adjust to PALM-SALSA bins
psd = zeros( length( time ), sum( zi ), length( dmid ) );
for k = 1:sum(zi) 
  conc = eval( sprintf( 'dNdlogDp_%dm', height(k) ) ) .* dlogDp; % 1/m3
  for b = 1:length( dmid )
    ind = ( dg(1:end-1) > dlims(b) ) & ( dg(2:end) <= dlims(b+1) );
    psd(:,k,b) = sum( conc(si:end-ei,ind), 2  );
    %fprintf('%d, %d, %d \n', k,b,sum( any(  psd(:,k,b)== 0) ) ) 
  end
  % Normalise by the total aerosol number
  ind = (dg < dlims(end) );
  %psd(:,k,:) = psd(:,k,:) .* ( sum( conc(si:end-ei,:), 2 ) ./ sum( psd(:,k,:), 3 ) ); 
end

% Then save
datamat.psd = psd;
datamat.dlims = dlims;
datamat.dmid = dmid;

fprintf('%s saved\n', fileout)

%%
figure(13)
loglog( dg, dNdlogDp_5m(4,:), 'k--' )
hold on
for k = 1:4
  loglog( dmid, reshape( psd(4,k,:),1,[]) ./ ( log10(dlims(2:end))-log10(dlims(1:end-1) ) ) )
end
legend('ADCHEM 5m','5m','20m','45m','80m')

%% 
figure(14)
hold on
for k = 1:4
  loglog( dmid, reshape( psd(4,k,:),1,[]) ./ ( log10(dlims(2:end))-log10(dlims(1:end-1) ) ) )
  legend('5m','20m','45m','80m')
end
set(gca,'xscale','log','yscale','log')