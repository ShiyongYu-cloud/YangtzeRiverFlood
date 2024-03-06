function [TPR,FPR,AUC] = my_roc(ins,sed)
%function for constructing a roc curve by considering age error
%INPUT
%ins: structure containing instrumental data
%sed: structure containing sedimentary data
%OUTPUT
%TPR: true positive rate
%FPR: false positive rate
%AUC: area under curve
%h: discrimination threshold
%%
% identify all flood events
N = 500; %number of points to be evaluated
h = linspace(0,1,N);
TPR = zeros(N,1);
FPR = zeros(N,1);
ACC = zeros(N,1);
for j = 1:N
    flood = find_flood(sed,h(j));
    % extract sedimentary data during the instrumental period
    x.mean = double(ins.years);
    %find age error
    depth = interp1(sed.age_mean,sed.depth,1950-x.mean);
    x.min = round(1950-interp1(sed.depth,sed.age_min,depth));
    x.max = round(1950-interp1(sed.depth,sed.age_max,depth));
    %swap min and max ages
    temp = x.min;
    x.min = x.max;
    x.max = temp;
    %age convergs at the core top
    x.min(1) = x.mean(1);
    x.max(1) = x.mean(1); 
    %label the sedimentary flood events
    x.predt = cell(length(x.mean),1);
    for i = 1:length(x.mean)
        if ismember(x.mean(i),flood.age_mean) 
            x.predt{i} = 'flood';
        else
            x.predt{i} = 'non_flood';    
        end
    end
    %label the validation
    x.valid = cell(length(x.mean),1);
    for i = 1:length(x.mean)
        sigm = 2;
        id = find(ins.years>=x.min(i)+sigm & ins.years<=x.max(i)-sigm & strcmpi(ins.labels,x.predt(i))==1, 1);
        if ~isempty(id)  %sum(id) >= 1
            if strcmpi(x.predt(i),'flood') ==1
                x.valid{i} = 'TP';
            elseif strcmpi(x.predt(i),'non_flood') ==1
                x.valid{i} = 'TN';
            end
        else
            if strcmpi(x.predt(i),'flood') ==1
                x.valid{i} = 'FP';
            elseif strcmpi(x.predt(i),'non_flood') ==1
                x.valid{i} = 'FN';
            end 
        end
    end
    id_TP = strcmpi(x.valid,'TP')==1;
    id_TN = strcmpi(x.valid,'TN')==1;
    id_FP = strcmpi(x.valid,'FP')==1;
    id_FN = strcmpi(x.valid,'FN')==1;
    tpr = sum(id_TP)/(sum(id_TP)+sum(id_FN)); %true positive rate 
    fpr = sum(id_FP)/(sum(id_FP)+sum(id_TN)); %false positive rate
    acc = (sum(id_TP)+sum(id_TN))/(sum(id_TP)+sum(id_TN)+sum(id_FP)+sum(id_FN));
    TPR(j) = tpr;
    FPR(j) = fpr;  
    ACC(j) = acc;
end
%pad 0 and 1
TPR = [1;TPR;0]; 
FPR = [1;FPR;0];
%calculate area under curve
AUC = 0;
for i = 2:N
    AUC = AUC + (FPR(i-1)-FPR(i))*(TPR(i-1)+TPR(i))/2;
end   
end