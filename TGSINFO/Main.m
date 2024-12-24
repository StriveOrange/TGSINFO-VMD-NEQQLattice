%% TGSINFO is an improved version based on the INFO algorithm.

clear 
close all
clc

nP=30;          % Number of Population

Func_name='F1'; % Name of the test function, range from F1-F23

MaxIt=500;      % Maximum number of iterations

% Load details of the selected benchmark function
[lb,ub,dim,fobj]=BenchmarkFunctions(Func_name);

[Best_fitness,BestPositions,Convergence_curve] = TGSINFO(nP,MaxIt,lb,ub,dim,fobj);

%% Draw objective space

figure(2),
hold on
semilogy(Convergence_curve,'Color','r','LineWidth',1);
title('Convergence curve')
xlabel('Iteration');
ylabel('Best fitness obtained so far');
axis tight
grid off
box on
legend('TGSINFO')


