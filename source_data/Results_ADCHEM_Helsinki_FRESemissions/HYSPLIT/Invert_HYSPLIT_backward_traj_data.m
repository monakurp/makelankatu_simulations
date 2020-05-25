for j=0:23
date=17061400+j

input = sprintf('C:/hysplit4/working/HYSPLIT_HELSINKI_2017/HYSPLIT_HELSINKI_BW_%d',date);
fid=fopen(input);
HYSPLIT=textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines', 20); %,'delimiter','tab',
Lat=HYSPLIT{10};
fclose(fid);
N=169;
if length(Lat)<N
nr_headerlines=20-N+length(Lat);
input = sprintf('C:/hysplit4/working/HYSPLIT_HELSINKI_2017/HYSPLIT_HELSINKI_BW_%d',date);
fid=fopen(input);
HYSPLIT=textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines', nr_headerlines); %,'delimiter','tab',
fclose(fid);
end

data_bw=[HYSPLIT{3},HYSPLIT{4},HYSPLIT{5},HYSPLIT{6},HYSPLIT{9},HYSPLIT{10},HYSPLIT{11},HYSPLIT{12},HYSPLIT{13},HYSPLIT{14},HYSPLIT{15},HYSPLIT{16}];

input = sprintf('C:/hysplit4/working/HYSPLIT_HELSINKI_2017/HYSPLIT_HELSINKI_FW_%d',date);
fid=fopen(input);
HYSPLIT=textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines', 20); %,'delimiter','tab',
Lat=HYSPLIT{10};
fclose(fid);
N2=25;
if length(Lat)<N2
nr_headerlines=20-N2+length(Lat);
input = sprintf('C:/hysplit4/working/HYSPLIT_HELSINKI_2017/HYSPLIT_HELSINKI_FW_%d',date);
fid=fopen(input);
HYSPLIT=textscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines', nr_headerlines); %,'delimiter','tab',
fclose(fid);
end

data_fw=[HYSPLIT{3},HYSPLIT{4},HYSPLIT{5},HYSPLIT{6},HYSPLIT{9},HYSPLIT{10},HYSPLIT{11},HYSPLIT{12},HYSPLIT{13},HYSPLIT{14},HYSPLIT{15},HYSPLIT{16}];
data=data_bw(N:-1:1,:);
data(N+1:N+N2-1,:)=data_fw(2:N2,:);
Lat=data(:,6);
Long=data(:,7);

outp=sprintf('HYSPLIT_Helsinki_%d%s',date,'.txt');
fid=fopen(outp,'w');
formatspeclength=repmat('\t %12.3f',1,11);
fprintf(fid,['\n %12.3f',formatspeclength],data'); 
fclose(fid);
end