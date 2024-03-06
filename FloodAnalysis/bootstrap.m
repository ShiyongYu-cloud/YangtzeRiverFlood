function X = bootstrap(values,age_scale,bsamples)
%% function for bootstrapping N sequences of random number
%INPUT
%values: vector containing data to be bootstrapped
%age_scale: scale of year to be reported (e.g. BCE/CE, BP, or B2K)
%bsamples: number of bootstraps
%OUTPUT
%X: matrix containing bootstrapped data sequences
%%
M = length(values);
y = 1:M;
y = y';
P = zeros(M,M);
for lambda = 1:M
    p = poisspdf(y,lambda);
    p = p/sum(p);
    P(:,lambda) = p;
end  
% P = zeros(M,M);
% for i = 1:M
%     p = [1:i-1 M:-1:i];
%     p = p';
%     p = p/sum(p);
%     P(:,i) = p;
% end 
X = zeros(M,bsamples);
for i = 1:M
    X(i,:) = durand(values,P(:,i),1,bsamples);
end
if strcmpi(age_scale,'BCE/CE') == 1
   X = sort(X,'ascend'); 
elseif strcmpi(age_scale,'BP') == 1
   X = 1950 - X;
   X = sort(X,'descend'); 
end    
return;
%%
function x = durand(values,prob,M,N)
%% function for drawing M*N random numbers from a list at given probability
%INPUT:
%values: a list of numbers to be drawn
%prob: corresponding probability value of the numbers
%M: number of rows of random numbers 
%N: number of columns of random numbers
%OUTPUT:
%x: M*N random numbers  
%%
values = values(:);
% prob = prob(:);
% if sum(prob)~=1
%    prob = prob/sum(prob);
% end
L = length(prob);
K = M*N;
psup = cumsum(prob);
pinf = [0; psup(1:end-1)];
Pinf = kron(ones(1,K),pinf(:));
Psup = kron(ones(1,K),psup(:));
u = rand(1,K);
U = kron(ones(L,1),u);
C = (U>Pinf) & (U<Psup);
V = kron(values(:),ones(1,K));
X = V.*C;
x = sum(X);
x = reshape(x,M,N);
return;