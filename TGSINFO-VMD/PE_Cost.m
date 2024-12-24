function [MinPE] = PE_Cost(alpha, K, signal)
tau=0;        
DC=1;         
init=1;       
tol=1e-6;     
[u, ~, ~] = VMD(signal, alpha, tau, K, DC, init, tol);
%% The minimum permutation entropy fitness value is calculated
SE_singal=u;
m=size(SE_singal,1);  
M = 3;                
T = 1;                
PE=zeros(1,m);
for i=1:m
    [PE(i),~] = PermutationEntropy(SE_singal(i,:),M,T);
end
MinPE=min(PE);
end