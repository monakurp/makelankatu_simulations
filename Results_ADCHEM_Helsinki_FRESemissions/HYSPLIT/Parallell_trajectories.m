% Code to calculate parallel trajectories along HYSPLIT air mass
% trajectories. Enable 2D-emissions along trajectories.
close all
clear all
clc
for jjjj=1
date = 2007032000+(jjjj-1); % chose trajectory 2806 2810
input = sprintf('HYSPLIT_Cambridge_%d%s',date,'.txt');
HYSPLIT=load(input);

nm=km2nm([0.5:1:19.5]);
lat_traj=HYSPLIT(:,6);
long_traj=HYSPLIT(:,7);

for j=1:length(lat_traj)
for i=1:20    
if j<length(lat_traj)
AZ = azimuth(lat_traj(j),long_traj(j),lat_traj(j+1),long_traj(j+1)); % computes the great circle bearing AZ (direction of trajectory)
end

deg = nm2deg(nm); % Degree distance of 1 km
deg_perp1=AZ-90; % degree perpendicular to the trajecory
deg_perp2=AZ+90; % degree perpendicular to the trajecory

[lat1(j,i),long1(j,i)]=reckon(lat_traj(j),long_traj(j),deg(i),deg_perp1);
[lat2(j,i),long2(j,i)]=reckon(lat_traj(j),long_traj(j),deg(i),deg_perp2);
end
end
Cambridge=[52.202738,0.120202]
London=[51.509865, -0.118092];    
Lat_2D=[lat1(:,20:-1:1),lat2];
Long_2D=[long1(:,20:-1:1),long2];
N=length(Lat_2D(:,1));
traj_1h=0:N-1; % length of 1 hour resolution data
traj_min=0:1/60:max(traj_1h); % length of 1 min resolution data [hours]
Lat_2D_high_res=interp1(traj_1h,Lat_2D,traj_min);
Long_2D_high_res=interp1(traj_1h,Long_2D,traj_min);


figure(19)
if jjjj==1
set(gcf,'paperorientation','portrait')
set(gcf,'paperunits','centimeters');
set(gcf,'papertype','a4letter');
set(gcf,'paperposition',[0 0 20 10]);
gray2=[1:-.1:0; 1:-0.1:0; 1:-0.1:0;]';

  latlim = [50 80];
  lonlim = [-20 20];
  worldmap(latlim, lonlim)
 geoshow('landareas.shp')%, 'FaceColor', [1.0 1.0 1.0]);
 hold on
 geoshow(Cambridge(1), Cambridge(2), 'DisplayType', 'Point', 'Marker', 'o', 'Color', 'black')
 hold on
 geoshow(London(1), London(2), 'DisplayType', 'Point', 'Marker', 'o', 'Color', 'b')
end
geoshow(Lat_2D_high_res(1:7381,20),Long_2D_high_res(1:7381,20),'DisplayType', 'Line', 'Color', 'r')
hold on
geoshow(Lat_2D_high_res(6600:7381,1),Long_2D_high_res(6600:7381,1),'DisplayType', 'Line', 'Color', 'k')
hold on
geoshow(Lat_2D_high_res(6600:7381,40),Long_2D_high_res(6600:7381,40),'DisplayType', 'Line', 'Color', 'k')
hold on
geoshow([Lat_2D_high_res(6600,40),Lat_2D_high_res(6600,1)],[Long_2D_high_res(6600,40),Long_2D_high_res(6600,1)],'DisplayType', 'Line', 'Color', 'k')
hold on
geoshow([Lat_2D_high_res(7381,40),Lat_2D_high_res(7381,1)],[Long_2D_high_res(7381,40),Long_2D_high_res(7381,1)],'DisplayType', 'Line', 'Color', 'k')

bla
outp=sprintf('HYSPLIT_Cambridge_2D_lat_long_%d%s',date,'.txt');
fid = fopen(outp,'w');
formatspeclength=repmat('\t %3.5f',1,79);
fprintf(fid,['\n %3.5f',formatspeclength],[lat1(:,20:-1:1),lat2(:,:),long1(:,20:-1:1),long2(:,:)]');
fclose(fid);
end
