%calculating annual maximum of consecutive x-month mean 
%vertically integrated moisture flux and convergence in the
%Yangtze River basin from ModE-Sim
clc;
clear;
%% read and concatenate data
q1420 = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','q');
q1850 = ncread('ModE-Sim_set_1850-1_ensmean_q_10000_abs_1850-2009_mon.nc','q');
u1420 = ncread('ModE-Sim_set_1420-3_ensmean_u_10000_abs_1420-1849_mon.nc','u');
u1850 = ncread('ModE-Sim_set_1850-1_ensmean_u_10000_abs_1850-2009_mon.nc','u');
v1420 = ncread('ModE-Sim_set_1420-3_ensmean_v_10000_abs_1420-1849_mon.nc','v');
v1850 = ncread('ModE-Sim_set_1850-1_ensmean_v_10000_abs_1850-2009_mon.nc','v');
lat = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','lat');
lon = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','lon');
lev = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','plev');
lev = double(lev);
lev = lev/100; %convert to hPa
[LON,LAT] = meshgrid(lon,lat);
q = cat(4,q1420,q1850);
u = cat(4,u1420,u1850);
v = cat(4,v1420,v1850);
%% transpose to lat-lon-lev-time
q = permute(q,[2 1 3 4]);
u = permute(u,[2 1 3 4]);
v = permute(v,[2 1 3 4]);
%% calculate JJA mean
time = 1:size(q,4);
L = numel(time)/12; %number of years
Q = reshape(q,numel(lat),numel(lon),numel(lev),12,L);
U = reshape(u,numel(lat),numel(lon),numel(lev),12,L);
V = reshape(v,numel(lat),numel(lon),numel(lev),12,L);
Q_mean = mean(Q,5);
U_mean = mean(U,5);
V_mean = mean(V,5);
q_prim = Q - Q_mean;
u_prim = U - U_mean;
v_prim = V - V_mean;
% calculate mean moisture flux and convergence
conv_mean = zeros(numel(lat),numel(lon),12);
for j = 1:12 
    [~,~,~,~,conv_mean(:,:,j)] = moisture_transport(lon,lat,lev,U_mean(:,:,:,j),V_mean(:,:,:,j),Q_mean(:,:,:,j));
end
%x = 1; %annual maximum monthly flux and convergence
x = 3; %consecutive 3-month mean flux and convergence
% calculate variation of moisture flux and convergence
max_totlm = zeros(numel(lat),numel(lon),L);
max_vprim = zeros(numel(lat),numel(lon),L);
max_qprim = zeros(numel(lat),numel(lon),L);
max_qvprm = zeros(numel(lat),numel(lon),L);
for i = 1:L
    c_totlm = zeros(numel(lat),numel(lon),12);
    c_vprim = zeros(numel(lat),numel(lon),12);
    c_qprim = zeros(numel(lat),numel(lon),12);
    c_qvprm = zeros(numel(lat),numel(lon),12);
    for j = 1:12
        [~,~,~,~,c_totlm(:,:,j)] = moisture_transport(lon,lat,lev,U(:,:,:,j,i),V(:,:,:,j,i),Q(:,:,:,j,i));
        [~,~,~,~,c_vprim(:,:,j)] = moisture_transport(lon,lat,lev,u_prim(:,:,:,j,i),v_prim(:,:,:,j,i),Q_mean(:,:,:,j));
        [~,~,~,~,c_qprim(:,:,j)] = moisture_transport(lon,lat,lev,U_mean(:,:,:,j),V_mean(:,:,:,j),q_prim(:,:,:,j,i));
        [~,~,~,~,c_qvprm(:,:,j)] = moisture_transport(lon,lat,lev,u_prim(:,:,:,j,i),v_prim(:,:,:,j,i),q_prim(:,:,:,j,i));
    end   
    %calculate consecutive x-month mean moisture flux and convegence
    c_totlm = c_totlm - conv_mean; % total minus mean
    c_totlm_month = movmean(c_totlm,x,3,'omitnan'); 
    c_vprim_month = movmean(c_vprim,x,3,'omitnan'); 
    c_qprim_month = movmean(c_qprim,x,3,'omitnan'); 
    c_qvprm_month = movmean(c_qvprm,x,3,'omitnan'); 
    % find annual maximum of consecutive x-month mean flux and convergence
    max_totlm(:,:,i) = max(c_totlm_month,[],3); %find max along months
    max_vprim(:,:,i) = max(c_vprim_month,[],3); %find max along months
    max_qprim(:,:,i) = max(c_qprim_month,[],3); %find max along months
    max_qvprm(:,:,i) = max(c_qvprm_month,[],3); %find max along months
 end
%% calculate area-weighted mean
mask = [110 25;120 25;120 35;110 35;100 25];
years = 1420:2009;
years = years';
conv_totlm = zeros(L,1);
conv_vprim = zeros(L,1);
conv_qprim = zeros(L,1);
conv_qvprm = zeros(L,1);
for i = 1:L
    [lon_small,lat_small,tl_small] = data_mask2D(max_totlm(:,:,i),lon,lat,mask);
    [~,~,vp_small] = data_mask2D(max_vprim(:,:,i),lon,lat,mask);
    [~,~,qp_small] = data_mask2D(max_qprim(:,:,i),lon,lat,mask);
    [~,~,qv_small] = data_mask2D(max_qvprm(:,:,i),lon,lat,mask);
    conv_totlm(i) = area_mean(lon_small,lat_small,tl_small);
    conv_vprim(i) = area_mean(lon_small,lat_small,vp_small);
    conv_qprim(i) = area_mean(lon_small,lat_small,qp_small);
    conv_qvprm(i) = area_mean(lon_small,lat_small,qv_small);
end
%%
%% plot results
figure(1)
load YRBasin.txt;
rivers = shaperead('worldrivers', 'UseGeoCoords', true);  
area = [111 28; 115 28; 115 31; 111 31; 111 28]; %study area
m_proj('miller','long',[90 150],'lat',[15 45]); 
[~,h] = m_contourf(LON,LAT,mean(conv_mean(:,:,5:6),3),256);      %mean state   
set(h,'edgecolor','none');
colormap(flipud(m_colmap('diverging',256)));
hold on
%m_quiver(LON,LAT,mean(MAXfx,3),mean(MAXfy,3));
m_gshhs_l('color','k','linewidth',0.5);
hold on
% plot Yellow and Yangtze Rivers
for i = 1:length(rivers)
      if strcmpi(rivers(i).Name,'Yellow')||strcmpi(rivers(i).Name,'Yangtze')
          m_plot(rivers(i).Lon,rivers(i).Lat,'b','linewi',0.5);
          hold on
      end
end
hold on
x0 = YRBasin(:,1);
y0 = YRBasin(:,2);
m_plot(x0,y0,'y-','linewi',1); 
hold on
m_plot(area(:,1),area(:,2),'r-','linewi',1);
m_grid('linewi',1,'tickdir','in',...
     'xtick',[90 100 110 120 130 140 150],'ytick',[15 25 35 45]);
h = colorbar;
caxis([-8e-4 8e-4]);
h.Label.String = 'JJA VIMFC (kg m-2 s-1)';

%% bootstrap and smooth the results
sim = 1; 
L = 10;
nsamples = 5000;
id = years>=1490;
years = years(id);
Btotlm = bootstrap_series(movmean(conv_totlm(id),11),sim,L,nsamples);
Bvprim = bootstrap_series(movmean(conv_vprim(id),11),sim,L,nsamples);
Bqprim = bootstrap_series(movmean(conv_qprim(id),11),sim,L,nsamples);
Bqvprm = bootstrap_series(movmean(conv_qvprm(id),11),sim,L,nsamples);
%
totlm_CI = prctile(Btotlm',[2.5 25 50 75 97.5])';
vprim_CI = prctile(Bvprim',[2.5 25 50 75 97.5])';
qprim_CI = prctile(Bqprim',[2.5 25 50 75 97.5])';
qvprm_CI = prctile(Bqvprm',[2.5 25 50 75 97.5])';
totlm_mean = mean(Btotlm,2);
vprim_mean = mean(Bvprim,2);
qprim_mean = mean(Bqprim,2);
qvprm_mean = mean(Bqvprm,2);
figure(3)
subplot(4,1,1); plot(years,totlm_CI);
xlabel('years (CE)');
ylabel('Total-mean');
xlim([1490 2010]);
grid on;
subplot(4,1,2); plot(years,vprim_CI);
xlabel('years (CE)');
ylabel('Circulation');
xlim([1490 2010]);
grid on;
subplot(4,1,3); plot(years,qprim_CI);
xlabel('years (CE)');
xlim([1490 2010]);
grid on;
ylabel('Moisture');
subplot(4,1,4); plot(years,qvprm_CI);
xlabel('years (CE)');
ylabel('Interaction');
xlim([1490 2010]);
grid on;
YRS = [years;flipud(years);years(1)];
tot95 = [totlm_CI(:,1);flipud(totlm_CI(:,5));totlm_CI(1,1)];
tot50 = [totlm_CI(:,2);flipud(totlm_CI(:,4));totlm_CI(1,2)];
v95 = [vprim_CI(:,1);flipud(vprim_CI(:,5));vprim_CI(1,1)];
v50 = [vprim_CI(:,2);flipud(vprim_CI(:,4));vprim_CI(1,2)];
q95 = [qprim_CI(:,1);flipud(qprim_CI(:,5));qprim_CI(1,1)];
q50 = [qprim_CI(:,2);flipud(qprim_CI(:,4));qprim_CI(1,2)];
qv95 = [qvprm_CI(:,1);flipud(qvprm_CI(:,5));qvprm_CI(1,1)];
qv50 = [qvprm_CI(:,2);flipud(qvprm_CI(:,4));qvprm_CI(1,2)];