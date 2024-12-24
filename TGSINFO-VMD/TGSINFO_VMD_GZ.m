clc
clear 
close all

%% Data loading
filename = 'Data_GuangZhou.xlsx';
Data = xlsread(filename);
rows_with_empty_data = any(isnan(Data), 2);     
Data_cleaned = Data(~rows_with_empty_data, :);  
data=flip(Data_cleaned,1);
Signal=data';
Signal=Signal(1:1522); 

Sequencelength=size(Signal,2);

%% Algorithm parameter settings
nP=5;          
MaxIt=5;      
lb=[100,2];                          %upper
ub=[1.5*Sequencelength,15];          %lower
dim=2;
fobj=@PE_Cost;

%% The TGSINFO algorithm optimizes VMD parameters
[Best_fitness,BestPositions,Convergence_curve] = TGSINFO(nP,MaxIt,lb,ub,dim,fobj,Signal);
BestaAlpha=BestPositions(1,1);
BestaK=BestPositions(1,2);

save('GZ_Best_fitness.mat', 'Best_fitness');
save('GZ_BestPositions.mat', 'BestPositions');
save('GZ_Convergence_curve.mat', 'Convergence_curve');

%% VMD decomposition based on optimal parameters
tau=0;        % Time step of the dual ascent (use 0 for noise relaxation), default value
DC=1;         % Set the first mode as DC (0 frequency), default value
init=1;       % Initialization can be 0, 1, or 2: 0 starts from zero; 1 starts with uniform distribution; 2 starts with random distribution
tol=1e-6;     % Tolerance for convergence criterion; typically about 1e-6

[u, u_hat, omega] = VMD(Signal, BestaAlpha, tau, BestaK, DC, init, tol);

%% Draws the results of the decomposition
fs=1;              %Sampling frequency, which is the time interval between two data points in a time series, samples every day here
Ts=1/fs;           %Sampling period
L=size(Signal,2);  
t=(1:L)*Ts;        

figure(1);
n = size(u, 1);
subplot(n+1, 1, 1);
plot(t, Signal);
ylabel('Original carbon price', 'fontsize', 8, 'fontname', 'Times New Roman');
colors = lines(n);  
for imfIndex = 1:n
    subplot(n+1, 1, imfIndex + 1);
    plot(t, u(imfIndex, :), 'color', colors(imfIndex, :), 'linewidth', 1);
    ylabel(['IMF' num2str(imfIndex)], 'fontsize', 8, 'fontname', 'Times New Roman');
end
xlabel('Time\it/day', 'fontsize', 14, 'fontname', 'Times New Roman');

%% Draw the convergence curve
figure(2);
hold on;
semilogy(Convergence_curve, ...
         'Color', [0, 0.5, 0.8], ...  
         'LineWidth', 1.5, ...  
         'Marker', 'd', ...  
         'MarkerSize', 7, ...  
         'MarkerEdgeColor', [0, 0.4, 0.7], ...  
         'MarkerFaceColor', [1, 0.8431, 0]);  
title('GZ Convergence Curve', 'FontSize', 12, 'FontName', 'Times New Roman');  
xlabel('Iteration', 'FontSize', 10, 'FontName', 'Times New Roman');  
ylabel('Best Fitness Obtained So Far', 'FontSize', 10, 'FontName', 'Times New Roman');  
axis tight;
grid off;
box on;
legend('TGSINFO-VMD', 'FontSize', 10, 'FontName', 'Times New Roman');  
