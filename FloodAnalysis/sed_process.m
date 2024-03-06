function sed = sed_process(sed_data,age_model,low_pass,high_pass)
%function for pre-processing sedimentary data
%INPUT
%sed_data: name of sedimentary data file
%age_model: name of age-depth model
%low_pass: number of points for the movmin filter
%high_pass: number of points for the movmean filter
%OUTPUT
%sed: structure containing processed sedimentary record
%% reading sedimetnary data
h = fopen(sed_data);
sed_record = textscan(h,'%f %f','headerlines',1,'delimiter','\t');
fclose(h);
sed.depth = sed_record{1};
sed.EM = sed_record{2};
%% read age-depth model
h = fopen(age_model);
chron = textscan(h,'%f %f %f %f %f','headerlines',1,'delimiter','\t');
fclose(h);
model.depth = chron{1};
model.age_min = chron{2};
model.age_max = chron{3};
model.age_medn = chron{4};
model.age_mean = chron{5};
%% map depth to age
sed.age_min = interp1(model.depth,model.age_min,sed.depth,'linear','extrap');
sed.age_max = interp1(model.depth,model.age_max,sed.depth,'linear','extrap');
sed.age_medn = interp1(model.depth,model.age_medn,sed.depth,'linear','extrap');
sed.age_mean = interp1(model.depth,model.age_mean,sed.depth,'linear','extrap');
%% remove trend
sed.EM_trend = movmin(sed.EM,low_pass);
sed.EM_detrend = sed.EM - sed.EM_trend; 
%% remove local mean
sed.EM_mean = movmean(sed.EM_detrend,high_pass);
EM_demean = (sed.EM_detrend - sed.EM_mean);
%sed.EM_score = (sed.EM_detrend-sed.EM_mean)/std(sed.EM_detrend,1); %normalize
sed.EM_score = EM_demean;
%% get EM score by keeping only things above the local mean
sed.EM_score(sed.EM_score<0)=0;
sed.EM_score = sed.EM_score/max(sed.EM_score); %normalization
end