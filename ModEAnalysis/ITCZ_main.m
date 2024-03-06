%% Calculating annual maximum consecutive x-month mean latitude of 
%% the ITCZ from the production of ModE  
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
%% convert to monthly total precipitation
for i = 1:L_ra
    for j = 1:12
        ptot_ra(:,:,j,i) = ptot_ra(:,:,j,i)*days(j)*24*60*60;
    end
end
%% calculate annual max consecutive 3-month mean divergence and total precipitation
x = 3; %consecutive 3-month mean
max_conv_ra = zeros(numel(lat),numel(lon),L_ra);
max_conv_cl = zeros(numel(lat),numel(lon),L_ra);
max_conv_sm = zeros(numel(lat),numel(lon),L_ra);
for i = 1:L_ra
    %calculate monthly convergence
    conv_ra_mon = zeros(numel(lat),numel(lon),12);
    conv_cl_mon = zeros(numel(lat),numel(lon),12);
    conv_sm_mon = zeros(numel(lat),numel(lon),12);
    for j = 1:12
        [~,~,conv_ra_mon(:,:,j)] = div_spherical(lon,lat,uwnd_ra(:,:,j,i),vwnd_ra(:,:,j,i)); %calculate stream function
        [~,~,conv_cl_mon(:,:,j)] = div_spherical(lon,lat,uwnd_cl(:,:,j,i),vwnd_cl(:,:,j,i)); %calculate stream function
        [~,~,conv_sm_mon(:,:,j)] = div_spherical(lon,lat,uwnd_sm(:,:,j,i),vwnd_sm(:,:,j,i)); %calculate stream function
    end
    conv_ra_mon = -conv_ra_mon;
    conv_cl_mon = -conv_cl_mon;
    conv_sm_mon = -conv_sm_mon;
    %calcualte 3-month moving mean 
    movmean_ra = movmean(conv_ra_mon,x,3,'omitnan'); 
    movmean_cl = movmean(conv_cl_mon,x,3,'omitnan'); 
    movmean_sm = movmean(conv_sm_mon,x,3,'omitnan'); 
    %find annual max
    max_ra = max(movmean_ra(:,:,6:7),[],3);
    max_cl = max(movmean_cl(:,:,6:7),[],3);
    max_sm = max(movmean_sm(:,:,6:7),[],3);
    %spatial smoothing
    smooth_ra = conv2(max_ra,ones(3,1)/3,'same'); 
    smooth_cl = conv2(max_cl,ones(3,1)/3,'same'); 
    smooth_sm = conv2(max_sm,ones(3,1)/3,'same'); 
    %remove noise
    max_conv_ra(:,:,i) = wiener2(smooth_ra,[3 3]);      
    max_conv_cl(:,:,i) = wiener2(smooth_cl,[3 3]);  
    max_conv_sm(:,:,i) = wiener2(smooth_sm,[3 3]);  
end
%% identify the latitudinal location of the ITCZ
mask = [50 0;130 0;130 25;50 25]; % area of influence by ITCZ over South Asia
[lon_mask,lat_mask,~] = data_mask3D(max_conv_ra(:,:,1),lon,lat,mask);
[LON_mask,LAT_mask] = meshgrid(lon_mask,lat_mask);
itcz_lat_ra = zeros(numel(lon_mask),L_ra);
itcz_lat_cl = zeros(numel(lon_mask),L_ra);
itcz_lat_sm = zeros(numel(lon_mask),L_ra);
for i = 1:L_ra
    % extract data
    [~,~,conv_ra_small] = data_mask3D(max_conv_ra(:,:,i),lon,lat,mask); 
    [~,~,conv_cl_small] = data_mask3D(max_conv_cl(:,:,i),lon,lat,mask); 
    [~,~,conv_sm_small] = data_mask3D(max_conv_sm(:,:,i),lon,lat,mask); 
    % find max convergence along latitude
    [~,lat_id_ra] = max(conv_ra_small,[],1); 
    [~,lat_id_cl] = max(conv_cl_small,[],1);
    [~,lat_id_sm] = max(conv_sm_small,[],1); 
    % find max latitude
    max_lat_ra = LAT_mask(lat_id_ra)';
    max_lat_cl = LAT_mask(lat_id_cl)';
    max_lat_sm = LAT_mask(lat_id_sm)';
    % smoothiing
    max_lat_ra = smooth(lon_mask,max_lat_ra,0.3,'loess');
    max_lat_cl = smooth(lon_mask,max_lat_cl,0.3,'loess');
    max_lat_sm = smooth(lon_mask,max_lat_sm,0.3,'loess');
    itcz_lat_ra(:,i) = max_lat_ra;
    itcz_lat_cl(:,i) = max_lat_cl;
    itcz_lat_sm(:,i) = max_lat_sm;
end
%% plot mean state
figure(1)
rivers = shaperead('worldrivers', 'UseGeoCoords', true);  
%plot mean
m_proj('miller','long',[50 130],'lat',[0 40]); 
[~,h] = m_contourf(LON,LAT,mean(max_conv_ra,3),256);      %mean state   
set(h,'edgecolor','none');
colormap(flipud(m_colmap('diverging',256)));
hold on
m_plot(lon_mask,mean(itcz_lat_ra,2),'y-','linewidth',2);
hold on
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
     'xtick',[50 70 90 110 130],'ytick',[0 10 20 30 40]);
h = colorbar;
h.Label.String = 'Convergence (s^-1)';
caxis([-1e-5,1e-5]);
title('JJA 850 hPa wind convergence over South Asia, mean (1421-2008 CE)');
%
%% plot temporal variability
itcz_mean_ra = median(itcz_lat_ra,1)';    
itcz_mean_cl = median(itcz_lat_cl,1)';    
itcz_mean_sm = median(itcz_lat_sm,1)';    
years = 1421:1:2008;
years = years';
id1 = years<1800;
id2 = years>=1800;
% t-test of trend 
% [H_sm,P_sm] = ttest2(itcz_mean_sm(id1),itcz_mean_sm(id2));
% [H_cl,P_cl] = ttest2(itcz_mean_cl(id1),itcz_mean_cl(id2));
% [H_ra,P_ra] = ttest2(itcz_mean_ra(id1),itcz_mean_ra(id2));
% itcz_ra1 = (itcz_mean_ra(id1) - mean(itcz_mean_ra(id1)))/std(itcz_mean_ra(id1));
% itcz_ra2 = (itcz_mean_ra(id2) - mean(itcz_mean_ra(id2)))/std(itcz_mean_ra(id2));
% itcz_ra = [itcz_ra1; itcz_ra2];
itcz_ra = (itcz_mean_ra-movmean(itcz_mean_ra,3))/std(itcz_mean_ra);
%
figure(2)
subplot(2,2,1); %plot ModE_RA
plot(years,itcz_mean_ra);
hold on
plot(years(id1),ones(numel(years(id1)),1).*mean(itcz_mean_ra(id1)));
hold on
plot(years(id2),ones(numel(years(id2)),1).*mean(itcz_mean_ra(id2)));
xlim([1490 2010]);
ylim([13 19]);
grid on;
xlabel('Year (CE)');
ylabel('Max latitude (N)');
subplot(2,2,2); %plot ModE_Sim
plot(years,itcz_mean_sm);
hold on
plot(years(id1),ones(numel(years(id1)),1).*mean(itcz_mean_sm(id1)));
hold on
plot(years(id2),ones(numel(years(id2)),1).*mean(itcz_mean_sm(id2)));
xlim([1490 2010]);
ylim([13 19]);
grid on;
xlabel('Year (CE)');
ylabel('Max latitude (N)');
subplot(2,2,3); %plot ModE_RAClim
plot(years,itcz_mean_cl);
hold on
plot(years(id1),ones(numel(years(id1)),1).*mean(itcz_mean_cl(id1)));
hold on
plot(years(id2),ones(numel(years(id2)),1).*mean(itcz_mean_cl(id2)));
xlim([1490 2010]);
ylim([13 19]);
grid on;
xlabel('Year (CE)');
ylabel('Max latitude (N)');
subplot(2,2,4); %plot ModE_RA z-score
plot(years,itcz_ra);
xlim([1490 2010]);
ylim([-3 3]);
grid on;
xlabel('Year (CE)');
ylabel('Max latitude (N)');
%
%% calculate strong ITCZ frequency
events = years(itcz_ra>0.2);
[yrs,freq] = flood_freq(events,1480,2020,31);
%bootstrap
sim = 1;
L = 10;
nbsamples = 5000;
Xb = bootstrap_series(freq,sim,L,nbsamples); %bootstrap
CI = prctile(Xb',[2.5 25 50 75 97.5]); %calcualte confidence interval
CI = CI';
CI(CI<0) = 0;
CI95 = [CI(:,1);flipud(CI(:,5));CI(1,1)]; 
CI50 = [CI(:,2);flipud(CI(:,4));CI(1,2)]; 
YRS = [yrs;flipud(yrs);yrs(1)];
% plot frequency
figure(4)
fill(YRS,CI95,'r');
hold on
fill(YRS,CI50,'b');
hold on
plot(yrs,movmean(freq,31),'g');
xlim([1480 2020]);
xlabel('Year (CE)');
ylabel('Frequency');
%
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
%% calculate annual max consecutive 3-month moisture flux and convergence
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
years_sm = 1420:1:2009;
id_sm = ismember(years_sm,years_sm(itcz_mean_sm>18)); 
id_ra = ismember(years,years(itcz_ra>1)); 
jja_fluxx_anom = jja_fluxx - mean(jja_fluxx,3);
jja_fluxy_anom = jja_fluxy - mean(jja_fluxy,3);
jja_convg_anom = jja_convg - mean(jja_convg,3);
composite_fluxx_anom = mean(jja_fluxx_anom(:,:,id_sm),3);
composite_fluxy_anom = mean(jja_fluxy_anom(:,:,id_sm),3);
composite_convg_anom = mean(jja_convg_anom(:,:,id_sm),3);
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
composite_ptot_anom = mean(jja_ptot_anom(:,:,id_ra),3);
%
figure(5)
load YRBasin.txt;
area = [111 28; 115 28; 115 31; 111 31; 111 28];
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
caxis([-8e-5,8e-5]);
title('Composite JJA mean moisture flux and convergence anomalies');
%
subplot(1,2,2)
m_proj('miller','long',[90 150],'lat',[15 45]);
hold on
[~,h] = m_contourf(LON,LAT,-composite_ptot_anom,256);
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
caxis([-100,100]);
title('Composite JJA mean precipitation anomalies');