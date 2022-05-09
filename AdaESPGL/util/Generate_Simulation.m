function [ Sig , t, Simulation] = Generate_Simulation(Params)
% This function creates the simulation containing impulses, harmonic
% and noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Params: a struct contains all the parameters
%%% Global parameters
%       Params.random_seed:   The random state
%       Params.Fs:            The sampling frequency of the simulation signal
%       Params.N:             The length of the signal
%       Params.Fn:            The fault characteristic frequency
%       Params.F:             The resonance frequency
%       Params.mixture_ratio: The mixing ratio of [impulses, harmonic, noise]
%%% noise type
%       Params.noise_type:    The noise type can be 'Gaussian' or 'Laplacian'       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig:                        The generated signal
% t:                          The time index
% Author : Zhibin Zhao
% Place  : Xi'an Jiaotong University
% Email  : zhibinzhao1993@gmail.com
% Date   : 2018.6
rng('default')
rng(Params.random_seed)
Fs = Params.Fs;
F = Params.F;
N = Params.N;
Fn = Params.Fn;
Simulation = zeros(N , 1);
K = floor(Fs / Fn);
Num = floor(N / K);
slip = round(2000*(2*rand(Num,1)-1)*0.001); 
Interval = 10;
n = (0 : Interval-1) / Fs;
n = n';
for i = 1 : Num
    U = randi(10);
    Transient = zeros(Interval , 1); 
    for j = 1 : U
        A = randn(1) + 1;
        w = sqrt(10)*randn(1) + F;
        beta =sqrt(1) * randn(1);
        Transient = Transient + A*sin(2*pi*w*n + beta) .* exp(-900*n);
    end 
    if i == 1
        Simulation((i-1)*K+1:(i-1)*K+Interval) = Transient;
    else
        Simulation((i-1)*K+1+slip(i):(i-1)*K+slip(i)+Interval) = Transient;
    end
end

shift=randi(100);
Simulation=[Simulation(shift+1:end);Simulation(1:shift)];

% Simulation(1 : fix(N*0.4)) = 0;

%% Generate the Noise
noise_type = Params.noise_type;
switch noise_type
    case 'Gaussian'
        Noise = randn(N , 1);
    case 'Laplacian'
        Noise = laprnd(N, 1);
    otherwise
        error('Unknown method.')
end
%% Generate the Combined Signal
mixture_ratio = Params.mixture_ratio;
Sig = mixture_ratio(1)*Simulation(:) + mixture_ratio(2)*Noise(:);
t = (0 : N-1) / Fs;
t = t(:);
end

