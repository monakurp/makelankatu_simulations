clc
close all
clear all

for jjjj=1
if jjjj<=2
date = 17060821+jjjj
else
date = 17060900+jjjj-3
end    
input=sprintf('HYSPLIT_Helsinki_%d%s',date,'.txt');
HYSPLIT=load(input)
Lat=HYSPLIT(1:169,6);
Long=HYSPLIT(1:169,7);
Altitude=HYSPLIT(1:169,8);

if jjjj==1
figure(1)    
set(gcf,'paperorientation','portrait')
set(gcf,'paperunits','centimeters');
set(gcf,'papertype','a4letter');
set(gcf,'paperposition',[0 0 20 10]);
gray2=[1:-.1:0; 1:-0.1:0; 1:-0.1:0;]';

  latlim = [40 90];
  lonlim = [0 180];
  worldmap(latlim, lonlim)
 geoshow('landareas.shp')%, 'FaceColor', [1.0 1.0 1.0]);
 hold on
end
figure(1)
geoshow(Lat,Long,'DisplayType', 'Line', 'Color', [1-jjjj/25, 0, jjjj/25])

figure(2)
plot([1:length(Altitude)],Altitude,'color',[1-jjjj/25, 0, jjjj/25])
hold on
end
%end
