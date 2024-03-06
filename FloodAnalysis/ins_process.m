function ins = ins_process(instrumental_record,flood_stage,critical_h,record_end)
%function for processing instrumental data
%INPUT
%instrumental_record: name of instrumental data file
%flood_stage: local design flood stage
%critical_h: critical height above flood stage
%record_end: year the core was collected
%OUTPUT
%ins: structure constaining flood years and Qmax
%% read instrumetnal record
h = fopen(instrumental_record);
data = textscan(h,'%f %f %f','headerlines',1,'delimiter','\t');
fclose(h);
years = data{1};
Qmax = data{2};
Zmax = data{3};
id = years <= record_end;
years = years(id);
Zmax = Zmax(id);
%identify flood events from the instrumental record
%[pks,locs] = findpeaks(Zmax,double(years),'MinPeakProminence',0.5,'MinPeakDistance',1);
%refine flood years
z0 = flood_stage + critical_h;
%locs = locs(pks>=z0);
locs = years(Zmax>z0);
ins.flood_years = locs;
ins.years = years;
ins.Zmax = Zmax;
id = ismember(years,locs);
ins.flood_Qmax = Qmax(id);
ins.flood_Zmax = Zmax(id);
ins.flood_years = flipud(ins.flood_years);
ins.years = flipud(ins.years);
ins.Zmax = flipud(ins.Zmax);
ins.flood_Qmax = flipud(ins.flood_Qmax);
ins.flood_Zmax = flipud(ins.flood_Zmax);
end
