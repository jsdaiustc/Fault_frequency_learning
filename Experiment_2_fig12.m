clear;
close all;
addpath(genpath(fileparts(mfilename('fullpath'))));
data=load('data_for_fig9.mat');      
y=data.y;                            % load the dataset for Fig.9
Fs=data.Fs;                          % sampling frequence 
N=length(y);
t = (0 : N-1) / Fs;

%% Our method
y_envo= abs(hilbert(y))-mean(abs(hilbert(y))); % obtain the noisy fault impulse signal envelope
y_h=  hilbert(y_envo);                         % remove the redundant conjugate frequencies
f_sample=[0:2:1000];                           % set the grid
[res_x,res_sample] = fault_frequency_learning(y_h,f_sample,Fs);

%% GSSA
Q = 2;
r = 5;
J =10;
now = ComputeNow(N,Q,r,J,'radix2');
AH = @(Sig) tqwt_radix2(Sig, Q, r, J);
A = @(w) itqwt_radix2(w, Q, r , N);
lam = 1.0 * now;
rho = 1;
load Performance_Comparison_Combination_K_Index_Size5_Sigma6.mat
K1 = round(N*0.04); 
Method1.Name = 'WGL';
Method1.Initial_Size = 5;
Method1.SubName = 'MC';
Method1.gamma = 2;
Method1.window = 'gausswin';
z1 = IterGSS(y, A, AH, lam, rho, K1, Method1);
P_GSSA = real(A(z1));
y_GSSA_enve=abs(fft(abs(hilbert(P_GSSA)) -mean(abs(hilbert(P_GSSA)))  ))/(N/2);

 %% GSL
 [GSL_result] = GSL(y);
 our_GSL_enve=abs(fft(abs(hilbert(GSL_result)) -mean(abs(hilbert(GSL_result))) ))/(N/2);
 
%% AdaESPGL
Params.Fs            = Fs;     % The sampling frequency of the simulation signal
Params.N             = N;      % The length of the signal
Params.N1    = 4;              % The samples of one impulse
Params.M     = 4;              % The number of periods
Params.Fn_N  = 0;              % a vector which contains the period of each component (Fs / fc)
Params.mu    = 9.235e-4;       % The parameter related to sparsity within groups
Params.pen   = 'atan';         % The penalty function
Params.rho   = 1;              % The degree of nonconvex
Params.Nit   = 100;            % The number of iteration 
% Estimate noise
[C,L]=wavedec(y,5,'sym8');
c1=detcoef(C,L,1);
est_noise=median(abs(c1-median(c1)))/0.678;
Params.lam= 0.272*est_noise + 0.044; 
[AdaESPGL_result] = AdaESPGL(y, Params);
 y_AdaESPGL_enve=abs(fft(abs(hilbert(AdaESPGL_result)) -mean(abs(hilbert(AdaESPGL_result)))   ))/(N/2); 
 
%% BPD
N=Params.N ;
rho = 2;
Method.Name = 'L1';
k_sparsity=round(N*10/100);
BPD_result = IterGSS_modified(y, rho, k_sparsity, Method)';
y_BPD_enve=abs(fft(abs(hilbert(BPD_result)) -mean(abs(hilbert(BPD_result))) ))/(N/2);

%% %%%%%%%%%%%%%%% P-GSL
[P_GSL_result] = P_GSL(y, Fs);
our_PSBL_enve=abs(fft(abs(hilbert(P_GSL_result)) -mean(abs(hilbert(P_GSL_result))) ))/(N/2);


%% Plot figure
temp1=600; temp2=0.3;
F = ([1:N]-1)*Fs/N;
F2= F(1:2001);
figure(9);
subplot(3,2,1)
stem(res_sample,abs(res_x)/(2) ,'Marker','none')
axis([0 temp1 0  temp2])
ylabel('Amp.[m/s^2]')
title('Our method')


subplot(3,2,2)
plot(F2,  y_GSSA_enve(1:2001))
axis([0 temp1 0 temp2])
title('GSSA')

subplot(3,2,3)
plot(F2, our_GSL_enve(1:2001) );
axis([0 temp1 0 temp2])
ylabel('Amp.[m/s^2]')
title('GSL')

subplot(3,2,4)
plot(F2,  y_BPD_enve(1:2001))
axis([0 temp1 0 temp2])
title('BPD')

subplot(3,2,5)
plot(F2,  our_PSBL_enve(1:2001) )
axis([0 temp1 0 temp2])
xlabel('Frequency')
ylabel('Amp.[m/s^2]')
title('P-GSL')

subplot(3,2,6)
plot(F2,  y_AdaESPGL_enve(1:2001))
axis([0 temp1 0 temp2])
xlabel('Frequency')
title('AdaESPGL')
