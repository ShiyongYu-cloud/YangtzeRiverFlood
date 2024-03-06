function [x_em,y_qm] = EM_Qmax(flood,ins)
%function for probabilistic matching EM score to Qmax
%INPUT
%flood: constructure containing sedimentary flood events
%ins: structure containing instrumental flood events
%OUTPUT
%x_em: EM scores corresponding to the sedimentary flood events
%q_qm: Qmax corresponding to the sedimentary flood events
%%
%extract data during the instrumental period
id = flood.age_mean>=min(ins.flood_years);
x.em_score = flood.EM_score(id);
x.age_min = flood.age_min(id);
x.age_max = flood.age_max(id);
x.age_mean = flood.age_mean(id);
% find Qmax from the instrumental record
sigmma = abs((x.age_max-x.age_min))/4;
y_Qmax = zeros(numel(x.em_score),1);
for i = 1:numel(x.em_score)
    prob = exp(-0.5*((x.age_mean(i)-ins.flood_years)./sigmma(i)).^2);
    [prob_max,id] = max(prob);
    c_low = exp(-2);
    c_upp = exp(0);
    if prob_max>c_low && prob_max<c_upp
        y_Qmax(i) = ins.flood_Qmax(id);
    else
        y_Qmax(i) = NaN;
    end
end
x_em = x.em_score;
y_qm = y_Qmax;
% remove NaN
id = isnan(y_qm);
x_em(id) = [];
y_qm(id) = [];
x_em = sort(x_em);
y_qm = sort(y_qm);
end