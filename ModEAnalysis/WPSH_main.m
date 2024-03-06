%% Calculating annual maximum consecutive x-month mean stream function
%% from the production of ModE  
%% https://www.wdc-climate.de/ui/project?acronym=ModE
clear;
clc;
%% read ModE-RA data
lat = ncread('ModE-RA_ensmean_u_85000_anom_wrt_1901-2000_1421-2008_mon.nc','latitude');
lon = ncread('ModE-RA_ensmean_u_85000_anom_wrt_1901-2000_1421-2008_mon.nc','longitude');
time_ra = ncread('ModE-RA_ensmean_u_85000_anom_wrt_1901-2000_1421-2008_mon.nc','time');
uwnd_ra = ncread('ModE-RA_ensmean_u_85000_anom_wrt_1901-2000_1421-2008_mon.nc','u');
vwnd_ra = ncread('ModE-RA_ensmean_v_85000_anom_wrt_1901-2000_1421-2008_mon.nc','v');
ptot_ra = ncread('ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc','totprec');
u_clim = ncread('ModE-RA_ensmean_u_85000_clim_1901-2000_mon.nc','u');
v_clim = ncread('ModE-RA_ensmean_v_85000_clim_1901-2000_mon.nc','v');
p_clim = ncread('ModE-RA_ensmean_totprec_clim_1901-2000_mon.nc','totprec');
%% read ModE-RAclim data
uwnd_cl = ncread('ModE-RAclim_ensmean_u_85000_anom_1421-2008_mon.nc','u');
vwnd_cl = ncread('ModE-RAclim_ensmean_v_85000_anom_1421-2008_mon.nc','v');
%% read ModE-Sim data
time_sm1 = ncread('ModE-Sim_set_1420-3_ensmean_u_85000_abs_1420-1849_mon.nc','time');
time_sm2 = ncread('ModE-Sim_set_1850-1_ensmean_u_85000_abs_1850-2009_mon.nc','time');
uwnd_sm1 = ncread('ModE-Sim_set_1420-3_ensmean_u_85000_abs_1420-1849_mon.nc','u');
uwnd_sm2 = ncread('ModE-Sim_set_1850-1_ensmean_u_85000_abs_1850-2009_mon.nc','u');
vwnd_sm1 = ncread('ModE-Sim_set_1420-3_ensmean_v_85000_abs_1420-1849_mon.nc','v');
vwnd_sm2 = ncread('ModE-Sim_set_1850-1_ensmean_v_85000_abs_1850-2009_mon.nc','v');
%% transpose to lat-lon-time
ptot_ra = permute(ptot_ra,[2 1 3]);
uwnd_ra = permute(uwnd_ra,[2 1 3]);
vwnd_ra = permute(vwnd_ra,[2 1 3]);
uwnd_cl = permute(uwnd_cl,[2 1 3]);
vwnd_cl = permute(vwnd_cl,[2 1 3]);
uwnd_sm1 = permute(uwnd_sm1,[2 1 4 3]);
uwnd_sm2 = permute(uwnd_sm2,[2 1 4 3]);
vwnd_sm1 = permute(vwnd_sm1,[2 1 4 3]);
vwnd_sm2 = permute(vwnd_sm2,[2 1 4 3]);
p_clim = permute(p_clim,[2 1 3]);
u_clim = permute(u_clim,[2 1 3]);
v_clim = permute(v_clim,[2 1 3]);
%% convert lon from -180:180 to 0:360 
lon = [lon(97:end);180+(180+lon(1:96))]; 
[LON,LAT] = meshgrid(lon,lat);
days = [31 28 31 30 31 30 31 31 30 31 30 31];
L_ra = numel(time_ra)/12; %number of years
L_sm1 = numel(time_sm1)/12; %number of years
L_sm2 = numel(time_sm2)/12; %number of years
%% reshape data to lat-lon-month-year order
ptot_ra = reshape(ptot_ra,numel(lat),numel(lon),12,L_ra);
uwnd_ra = reshape(uwnd_ra,numel(lat),numel(lon),12,L_ra);
vwnd_ra = reshape(vwnd_ra,numel(lat),numel(lon),12,L_ra);
uwnd_cl = reshape(uwnd_cl,numel(lat),numel(lon),12,L_ra);
vwnd_cl = reshape(vwnd_cl,numel(lat),numel(lon),12,L_ra);
uwnd_sm1 = reshape(uwnd_sm1,numel(lat),numel(lon),12,L_sm1);
uwnd_sm2 = reshape(uwnd_sm2,numel(lat),numel(lon),12,L_sm2);
vwnd_sm1 = reshape(vwnd_sm1,numel(lat),numel(lon),12,L_sm1);
vwnd_sm2 = reshape(vwnd_sm2,numel(lat),numel(lon),12,L_sm2);
uwnd_sm = cat(4,uwnd_sm1(:,:,:,2:end),uwnd_sm2(:,:,:,1:end-1)); %combine ModE-Sim sub-datasets
vwnd_sm = cat(4,vwnd_sm1(:,:,:,2:end),vwnd_sm2(:,:,:,1:end-1)); %combine ModE-Sim sub-datasets
%% re-arrange data along 0-360 longitude
ptot_ra = [ptot_ra(:,97:end,:,:) ptot_ra(:,1:96,:,:)]; 
uwnd_ra = [uwnd_ra(:,97:end,:,:) uwnd_ra(:,1:96,:,:)]; 
vwnd_ra = [vwnd_ra(:,97:end,:,:) vwnd_ra(:,1:96,:,:)]; 
uwnd_cl = [uwnd_cl(:,97:end,:,:) uwnd_cl(:,1:96,:,:)]; 
vwnd_cl = [vwnd_cl(:,97:end,:,:) vwnd_cl(:,1:96,:,:)]; 
uwnd_sm = [uwnd_sm(:,97:end,:,:) uwnd_sm(:,1:96,:,:)]; 
vwnd_sm = [vwnd_sm(:,97:end,:,:) vwnd_sm(:,1:96,:,:)]; 
p_clim = [p_clim(:,97:end,:) p_clim(:,1:96,:)];
u_clim = [u_clim(:,97:end,:) u_clim(:,1:96,:)];
v_clim = [v_clim(:,97:end,:) v_clim(:,1:96,:)];
%% add climatology back to ModE-RA and ModE-RAClim 
ptot_ra = ptot_ra + p_clim; 
uwnd_ra = uwnd_ra + u_clim; 
vwnd_ra = vwnd_ra + v_clim; 
uwnd_cl = uwnd_cl + u_clim; 
vwnd_cl = vwnd_cl + v_clim; 
%% convert daily to monthly precipitation
for i = 1:L_ra
    for j = 1:12
        ptot_ra(:,:,j,i) = ptot_ra(:,:,j,i)*days(j)*24*60*60;
    end
end
%% calculate JJA mean and annual max consecutive 3-month mean stream function 
x = 3; %consecutive 3-month mean
max_stream_ra = zeros(numel(lat),numel(lon),L_ra);
max_stream_cl = zeros(numel(lat),numel(lon),L_ra);
max_stream_sm = zeros(numel(lat),numel(lon),L_ra);
jja_stream_ra = zeros(numel(lat),numel(lon),L_ra);
jja_stream_cl = zeros(numel(lat),numel(lon),L_ra);
jja_stream_sm = zeros(numel(lat),numel(lon),L_ra);
jja_ptot_ra = zeros(numel(lat),numel(lon),L_ra);
for i = 1:L_ra
    %calculate monthly stream function
    stream_ra_mon = zeros(numel(lat),numel(lon),12);
    stream_cl_mon = zeros(numel(lat),numel(lon),12);
    stream_sm_mon = zeros(numel(lat),numel(lon),12);
    for j = 1:12
        [~,~,stream_ra_mon(:,:,j)] = stream_fun(lon,lat,uwnd_ra(:,:,j,i),vwnd_ra(:,:,j,i)); %calculate stream function
        [~,~,stream_cl_mon(:,:,j)] = stream_fun(lon,lat,uwnd_cl(:,:,j,i),vwnd_cl(:,:,j,i)); %calculate stream function
        [~,~,stream_sm_mon(:,:,j)] = stream_fun(lon,lat,uwnd_sm(:,:,j,i),vwnd_sm(:,:,j,i)); %calculate stream function
    end
    %calculate JJA mean 
    jja_stream_ra(:,:,i) = mean(stream_ra_mon(:,:,6:8),3,'omitnan');
    jja_stream_cl(:,:,i) = mean(stream_cl_mon(:,:,6:8),3,'omitnan');
    jja_stream_sm(:,:,i) = mean(stream_sm_mon(:,:,6:8),3,'omitnan');
    %calculate 3-month moving mean
    movmean_ra = movmean(stream_ra_mon,x,3,'omitnan'); 
    movmean_cl = movmean(stream_cl_mon,x,3,'omitnan'); 
    movmean_sm = movmean(stream_sm_mon,x,3,'omitnan'); 
    %find annual max
    max_stream_ra(:,:,i) = max(movmean_ra(:,:,6:7),[],3);
    max_stream_cl(:,:,i) = max(movmean_cl(:,:,6:7),[],3);
    max_stream_sm(:,:,i) = max(movmean_sm(:,:,6:7),[],3);
end
%% plot mean state and variability of the WPSH
%calculate mean and normalized std
JJA_mean = mean(jja_stream_ra,3);
JJA_stdd = std(jja_stream_ra,0,3);
zonal_mean = repmat(mean(JJA_stdd,2),1,numel(lon));
zonal_stdd = repmat(std(JJA_stdd,0,2),1,numel(lon));
JJA_stdd = (JJA_stdd - zonal_mean)./zonal_stdd;
figure(1)
%plot mean
rivers = shaperead('worldrivers', 'UseGeoCoords', true);  
subplot(2,1,1)
m_proj('miller','long',[90 240],'lat',[0 60]); 
[~,h] = m_contourf(LON,LAT,JJA_mean,256);      %mean state   
set(h,'edgecolor','none');
colormap(m_colmap('diverging',256));
hold on
m_gshhs_l('color','k','linewidth',0.5);
hold on
% plot Yellow and Yangtze Rivers
for i = 1:length(rivers)
      if strcmpi(rivers(i).Name,'Yellow')||strcmpi(rivers(i).Name,'Yangtze')
          m_plot(rivers(i).Lon,rivers(i).Lat,'b','linewi',0.5);
          hold on
      end
end
m_grid('linewi',1,'tickdir','in',...
     'xtick',[90 120 150 180 210 240],'ytick',[0 20 40 60]);
h = colorbar;
h.Label.String = 'Stream function (m^2/s)';
caxis([-1e+7,2e+7]);
title('JJA 850 hPa stream function (1421-2008 CE), mean');
% plot 1 s.d.
subplot(2,1,2)
m_proj('miller','long',[90 240],'lat',[0 60]);
[~,h] = m_contourf(LON,LAT,JJA_stdd,256);     % variability
set(h,'edgecolor','none');
colormap(m_colmap('diverging',256));
hold on
m_gshhs_l('color','k','linewidth',0.5);
hold on
% plot Yellow and Yangtze Rivers
for i = 1:length(rivers)
      if strcmpi(rivers(i).Name,'Yellow')||strcmpi(rivers(i).Name,'Yangtze')
          m_plot(rivers(i).Lon,rivers(i).Lat,'b','linewi',0.5);
          hold on
      end
end
% plot area for EOF analysis
hold on
area = [110 10;150 10;150 40; 110 40; 110 10];
m_plot(area(:,1),area(:,2),'linewi',0.5,'color','k');
m_grid('linewi',1,'tickdir','in',...
     'xtick',[90 120 150 180 210 240],'ytick',[0 20 40 60]);
h = colorbar;
h.Label.String = 'Stream function (m^2/s)';
%caxis([0.5e+6,1.5e+6]);
caxis([-3,3]);
title('JJA 850 hPa stream function (1421-2008 CE), normalized standard deviation');%
%
%% conduct EOF analysis
WEIT_LAT = sqrt(cosd(LAT));
stream_anom_ra = zeros(numel(lat),numel(lon),L_ra);
stream_anom_cl = zeros(numel(lat),numel(lon),L_ra);
stream_anom_sm = zeros(numel(lat),numel(lon),L_ra);
for i = 1:L_ra
    stream_anom_ra(:,:,i) = (max_stream_ra(:,:,i) - mean(max_stream_ra,3)).*WEIT_LAT;
    stream_anom_cl(:,:,i) = (max_stream_cl(:,:,i) - mean(max_stream_cl,3)).*WEIT_LAT;
    stream_anom_sm(:,:,i) = (max_stream_sm(:,:,i) - mean(max_stream_sm,3)).*WEIT_LAT;
end
%mask = geomask(LAT,LON,[10 40],[110 180]); %smaller area for EOF
mask = geomask(LAT,LON,[10 30],[120 140]); %smaller area for EOF
M = 2; %keep first two modes
[~,pc_ra,expvar_ra] = eof(stream_anom_ra,M,'mask',mask); %eof analysis
[~,pc_cl,expvar_cl] = eof(stream_anom_cl,M,'mask',mask); %eof analysis
[~,pc_sm,expvar_sm] = eof(stream_anom_sm,M,'mask',mask); %eof analysis
%% normalize PCs
for i = 1:M
    pc_ra(i,:) = pc_ra(i,:)/sqrt(pc_ra(i,:)*pc_ra(i,:)');%normalize PCs to 1
    pc_cl(i,:) = pc_cl(i,:)/sqrt(pc_cl(i,:)*pc_cl(i,:)');%normalize PCs to 1
    pc_sm(i,:) = pc_sm(i,:)/sqrt(pc_sm(i,:)*pc_sm(i,:)');%normalize PCs to 1
end
years = 1421:1:2008;
years = years';
id1 = years<1850;
id2 = years>=1850;
pc11_ra = (pc_ra(1,id1) - mean(pc_ra(1,id1)))/std(pc_ra(1,id1));
pc12_ra = (pc_ra(1,id2) - mean(pc_ra(1,id2)))/std(pc_ra(1,id2));
pc1_ra = [pc11_ra pc12_ra];
% t-test of trend 
[H_sm,P_sm] = ttest2(pc_sm(1,id1),pc_sm(1,id2));
[H_cl,P_cl] = ttest2(pc_cl(1,id1),pc_cl(1,id2));
[H_ra,P_ra] = ttest2(pc_ra(1,id1),pc_ra(1,id2));
%% regress stream function anomalies onto the PCs of EOFs
reg_eof = zeros(numel(lat),numel(lon),M);
for i = 1:numel(lat)
    for j = 1:numel(lon)
        data = permute(stream_anom_ra(i,j,:),[3,1,2]);
        for k = 1:M
            rho = corr(data,pc_ra(k,:)');
            reg_eof(i,j,k) = rho*std(stream_anom_ra(i,j,:))/std(pc_ra(k,:));
        end
    end
end
%
%% plot EOF modes
area = [111 28; 115 28; 115 31; 111 31; 111 28]; % Jianghan Plain
load YRBasin.txt;
figure(2)
subplot(2,1,1)
m_proj('miller','long',[90 210],'lat',[0 60]);
hold on
[~,h] = m_contourf(LON,LAT,reg_eof(:,:,1),256);
set(h, 'edgecolor','none');
colormap(m_colmap('diverging',256));
hold on
m_gshhs_l('color','k','linewidth',0.5);
hold on
%plot Yangtze River basin shape
x0 = YRBasin(:,1);
y0 = YRBasin(:,2);
m_plot(x0,y0,'y-','linewi',2); 
hold on
% plot Yellow and Yangtze Rivers
for j = 1:length(rivers)
      if strcmpi(rivers(j).Name,'Yellow')||strcmpi(rivers(j).Name,'Yangtze')
          m_plot(rivers(j).Lon,rivers(j).Lat,'b','linewi',0.5);
          hold on
      end
end
% plot study area    
hold on
m_plot(area(:,1),area(:,2),'r-','linewi',1);    
m_grid('linewi',1,'tickdir','in',...
 'xtick',[90 120 150 180 210 240],'ytick',[0 20 40 60]);
h = colorbar;
h.Label.String = 'm^2/s';
caxis([0.5e+7,2.5e+7]);
title(['JJA 850 hPa stream function anomalies regressed onto PC ',...
    num2str(1)]);
%
subplot(2,1,2)
m_proj('miller','long',[90 210],'lat',[0 60]);
hold on
[~,h] = m_contourf(LON,LAT,reg_eof(:,:,2),256);
set(h, 'edgecolor','none');
colormap(m_colmap('diverging',256));
hold on
m_gshhs_l('color','k','linewidth',0.5);
hold on
%plot Yangtze River basin shape
x0 = YRBasin(:,1);
y0 = YRBasin(:,2);
m_plot(x0,y0,'g-','linewi',1); 
hold on
% plot Yellow and Yangtze Rivers
for j = 1:length(rivers)
      if strcmpi(rivers(j).Name,'Yellow')||strcmpi(rivers(j).Name,'Yangtze')
          m_plot(rivers(j).Lon,rivers(j).Lat,'b','linewi',0.5);
          hold on
      end
end
% plot study area    
hold on
m_plot(area(:,1),area(:,2),'r-','linewi',1);    
m_grid('linewi',1,'tickdir','in',...
 'xtick',[90 120 150 180 210 240],'ytick',[0 20 40 60]);
h = colorbar;
h.Label.String = 'm^2/s';
caxis([-2e+7,2e+7]);
title(['JJA 850 hPa stream function anomalies regressed onto PC ',...
    num2str(2)]);
%% plot PCs
figure(3)
subplot(2,2,1) %plot ModE-RA
plot(years,pc_ra(1,:));
hold on
plot(years(id1),ones(numel(years(id1)),1).*mean(pc_ra(1,id1)));
hold on
plot(years(id2),ones(numel(years(id2)),1).*mean(pc_ra(1,id2)));
xlim([1490 2010]);
ylim([-0.2 0.2]);
xlabel('Year (CE)');
ylabel('PC1');
grid on
subplot(2,2,2) %plot ModE-Sim
plot(years,pc_sm(1,:));
hold on
plot(years(id1),ones(numel(years(id1)),1).*mean(pc_sm(1,id1)));
hold on
plot(years(id2),ones(numel(years(id2)),1).*mean(pc_sm(1,id2)));
xlim([1490 2010]);
ylim([-0.2 0.2]);
xlabel('Year (CE)');
ylabel('PC1');
grid on
subplot(2,2,3) %plot ModE-RAClim
plot(years,pc_cl(1,:));
hold on
plot(years(id1),ones(numel(years(id1)),1).*mean(pc_cl(1,id1)));
hold on
plot(years(id2),ones(numel(years(id2)),1).*mean(pc_cl(1,id2)));
xlim([1490 2010]);
ylim([-0.2 0.2])
xlabel('Year (CE)');
ylabel('PC1');
grid on
subplot(2,2,4); % plot pc1 z-score of ModE-RA
plot(years,pc1_ra);
xlim([1490 2010])
ylim([-4 4])
xlabel('Year (CE)');
ylabel('PC1 z-score');
grid on
%% calculate strong WPSH frequency from PC1
[pkt,yrs] = findpeaks(pc1_ra,years,'MinPeakHeight',1e-5,'Threshold',1e-5,'MinPeakDistance',2);
age_scale ='BCE/CE';
bsamples = 5000;
A = 1480;
B = 2020;
age = (A:B)';
W = 31;
FREQ_PC = zeros(numel(age),bsamples);
pc_boot = bootstrap(yrs,age_scale,bsamples);
for i = 1:bsamples
    f = unique(pc_boot(:,i));
    [~,freq] = flood_freq(pc_boot(:,i),A,B,W);
    FREQ_PC(:,i) = movmean(freq,41);
end
freq_CI = prctile(FREQ_PC',[2.5 25 50 75 97.5]);
freq_CI = freq_CI';
YRS = [age; flipud(age);age(1)];
CI_95 = [freq_CI(:,1); flipud(freq_CI(:,5));freq_CI(1,1)];
CI_50 = [freq_CI(:,2); flipud(freq_CI(:,4));freq_CI(1,2)];
%
figure(4)
fill(YRS,CI_95,'r');
hold on
fill(YRS,CI_50,'b');
hold on
plot(age,freq_CI(:,3),'y');
xlim([1480 2020]);
ylim([0 0.15]);
grid on
xlabel('Year (CE)');
ylabel('Frequency');
%% calculate moisture flux and convergence from ModE-Sim
% read and concatenate data
q1420 = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','q');
q1850 = ncread('ModE-Sim_set_1850-1_ensmean_q_10000_abs_1850-2009_mon.nc','q');
u1420 = ncread('ModE-Sim_set_1420-3_ensmean_u_10000_abs_1420-1849_mon.nc','u');
u1850 = ncread('ModE-Sim_set_1850-1_ensmean_u_10000_abs_1850-2009_mon.nc','u');
v1420 = ncread('ModE-Sim_set_1420-3_ensmean_v_10000_abs_1420-1849_mon.nc','v');
v1850 = ncread('ModE-Sim_set_1850-1_ensmean_v_10000_abs_1850-2009_mon.nc','v');
lat_sm = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','lat');
lon_sm = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','lon');
lev_sm = ncread('ModE-Sim_set_1420-3_ensmean_q_10000_abs_1420-1849_mon.nc','plev');
lev_sm = double(lev_sm);
lev_sm = lev_sm/100; %convert to hPa
[LON_sm,LAT_sm] = meshgrid(lon_sm,lat_sm);
q = cat(4,q1420,q1850);
u = cat(4,u1420,u1850);
v = cat(4,v1420,v1850);
%% transpose to lat-lon-lev-time
q = permute(q,[2 1 3 4]);
u = permute(u,[2 1 3 4]);
v = permute(v,[2 1 3 4]);
time_sm = 1:size(q,4);
L_sm = numel(time_sm)/12; %number of years
Q = reshape(q,numel(lat_sm),numel(lon_sm),numel(lev_sm),12,L_sm);
U = reshape(u,numel(lat_sm),numel(lon_sm),numel(lev_sm),12,L_sm);
V = reshape(v,numel(lat_sm),numel(lon_sm),numel(lev_sm),12,L_sm);
%% calculate jja mean moisture flux and convergence
jja_fluxx = zeros(numel(lat_sm),numel(lon_sm),L_sm);
jja_fluxy = zeros(numel(lat_sm),numel(lon_sm),L_sm);
jja_convg = zeros(numel(lat_sm),numel(lon_sm),L_sm);
for i = 1:L_sm
    fluxx_mon = zeros(numel(lat_sm),numel(lon_sm),12);
    fluxy_mon = zeros(numel(lat_sm),numel(lon_sm),12);
    convg_mon = zeros(numel(lat_sm),numel(lon_sm),12);
    for j = 1:12
        [~,~,fluxx_mon(:,:,j),fluxy_mon(:,:,j),convg_mon(:,:,j)] = moisture_transport(lon_sm,lat_sm,lev_sm,U(:,:,:,j,i),V(:,:,:,j,i),Q(:,:,:,j,i));
    end   
    %calculate jja mean
    jja_fluxx(:,:,i) = mean(fluxx_mon(:,:,6:8),3,'omitnan'); 
    jja_fluxy(:,:,i) = mean(fluxy_mon(:,:,6:8),3,'omitnan'); 
    jja_convg(:,:,i) = mean(convg_mon(:,:,6:8),3,'omitnan'); 
end
%% calculate composite jja moisture flux and convergence anonalies
id = ismember(years,yrs); 
jja_fluxx_anom = jja_fluxx - mean(jja_fluxx,3);
jja_fluxy_anom = jja_fluxy - mean(jja_fluxy,3);
jja_convg_anom = jja_convg - mean(jja_convg,3);
composite_fluxx_anom = mean(jja_fluxx_anom(:,:,id),3);
composite_fluxy_anom = mean(jja_fluxy_anom(:,:,id),3);
composite_convg_anom = mean(jja_convg_anom(:,:,id),3);
% add quiver scale
fx_rf = zeros(size(LON_sm));
fy_rf = zeros(size(LAT_sm));
fx_rf(10,30) = 5; % 2 kg/m/s
%% calculate JJA total precipitation anomalies from ModE-RA
jja_ptot = zeros(numel(lat),numel(lon),L_ra);
for i = 1:L_ra
    jja_ptot(:,:,i) = sum(ptot_ra(:,:,6:8,i),3);
end
jja_ptot_anom = jja_ptot - mean(jja_ptot,3);
composite_ptot_anom = mean(jja_ptot_anom(:,:,id),3);
%
figure(5)
subplot(1,2,1)
m_proj('miller','long',[90 150],'lat',[15 45]);
hold on
[~,h] = m_contourf(LON_sm,LAT_sm,composite_convg_anom,256);
set(h, 'edgecolor','none');
colormap(flipud(m_colmap('diverging',256)));
hold on
m_quiver(LON_sm,LAT_sm,composite_fluxx_anom,composite_fluxy_anom,'k','linewidth',0.5);
hold on
m_quiver(LON_sm,LAT_sm,fx_rf,fy_rf,'r','linewidth',1); %add quiver scale
hold on
m_gshhs_l('color','k','linewidth',0.5);
hold on
%plot Yangtze River basin shape
x0 = YRBasin(:,1);
y0 = YRBasin(:,2);
m_plot(x0,y0,'g-','linewi',1); 
hold on
% plot Yellow and Yangtze Rivers
for j = 1:length(rivers)
      if strcmpi(rivers(j).Name,'Yellow')||strcmpi(rivers(j).Name,'Yangtze')
          m_plot(rivers(j).Lon,rivers(j).Lat,'b','linewi',0.5);
          hold on
      end
end
% plot study area    
hold on
m_plot(area(:,1),area(:,2),'r-','linewi',1);    
m_grid('linewi',1,'tickdir','in',...
 'xtick',[90 100 110 120 130 140 150],'ytick',[15 25 35 45]);
h = colorbar('southoutside');
h.Label.String = 'Moisture convergence anomalies (10^-5 kg m^-2 s^-1)';
caxis([-2e-5,2e-5]);
title('Composite JJA mean moisture flux and convergence anomalies');
%
subplot(1,2,2)
m_proj('miller','long',[90 150],'lat',[15 45]);
hold on
[~,h] = m_contourf(LON,LAT,composite_ptot_anom,256);
set(h, 'edgecolor','none');
colormap(flipud(m_colmap('diverging',256)));
hold on
m_gshhs_l('color','k','linewidth',0.5);
hold on
%plot Yangtze River basin shape
x0 = YRBasin(:,1);
y0 = YRBasin(:,2);
m_plot(x0,y0,'g-','linewi',1); 
hold on
% plot Yellow and Yangtze Rivers
for j = 1:length(rivers)
      if strcmpi(rivers(j).Name,'Yellow')||strcmpi(rivers(j).Name,'Yangtze')
          m_plot(rivers(j).Lon,rivers(j).Lat,'b','linewi',0.5);
          hold on
      end
end
% plot study area    
hold on
m_plot(area(:,1),area(:,2),'r-','linewi',1);    
m_grid('linewi',1,'tickdir','in',...
 'xtick',[90 100 110 120 130 140 150],'ytick',[15 25 35 45]);
h = colorbar('southoutside');
h.Label.String = 'Precipitation anomalies (mm)';
caxis([-30,30]);
title('Composite JJA mean precipitation anomalies');