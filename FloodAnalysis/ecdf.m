function [X,F] = ecdf(x,Xmin,Xmax,M)
%% function for estimating empirical cdf
%INPUT
%x: array containing data for estimating cdf 
%Xmin: lower bound of the cdf
%Xmax: upper bound of the cdf
%M: number of data points at which cdf to be evaluated
%
%OUTPUToptions.
%X: data points at which cdf was estimated 
%F: cdfF at the data points
%% generate points at which cdf will be estimated
if Xmin > min(x)
    Xmin = min(x); %replace Xmin with min(x)
end
if Xmax < max(x)
    Xmax = max(x); %replace Xmax with max(x)
end
X = linspace(Xmin,Xmax,M);
X = X';
N = length(x);
F = zeros(M,1);
%% estimating CDF by definition
for i = 1:M
    p = 0;              % True Probability
    q = 0;              % False Probability
    for j = 1:N
        if x(j)<=X(i)   % Definition of CDF
            p = p + 1;
        else
            q = q + 1;
        end
    end
    F(i) = p/(p + q);   % Calulating Probability
end
%% remove duplicated values
[f,id,~] = unique(F,'stable');
x = X(id);
%% pad duplicated data through interpolation and extrapolation
id = f > 0 & f < 1;
x = [Xmin;x(id);Xmax];
f = [0;f;1];
[x,id,~] = unique(x,'stable');
f = f(id);
[f,id,~] = unique(f,'stable');
x = x(id);
F = interp1(x,f,X);
F = (F-F(1))/(F(end)-F(1)); %scale to [0 1]
end