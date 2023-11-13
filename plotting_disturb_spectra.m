clear all%#ok
close all
clc

rng(481516);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the disturbance spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- operating frequencies ---
f_camp_ESO = 1e4;          % ESO and plant operating frequency [Hz]
T_camp_ESO = 1/f_camp_ESO; % ESO sampling period [s]
T_camp_sys = T_camp_ESO;   % plant sampling period [s]
% --- operating frequencies ---

% ---- simulation settings ----
simTime = 5; % time horizon: [0 simTime]
nsample = floor(simTime/T_camp_sys);
TimeSamples = linspace(0, simTime, nsample);

p = 7 ; % number of outputs
% ---- simulation settings ----

% --- additive output disturbances ----
disturbances = DT_disturb_create(p, nsample); 
% --- additive output disturbances ----

% ---- load the ADRMPC results for the case N=8 ----
load RESULTS_MPC_ADRC_10kHz_N8.mat ErrY Nstart
% ErrY   <--> output tracking error
% Nstart <--> avoiding to analyse the initial transient


nY = p;
figure('Name','Output Disturbances - no ARDMPC', ...
    'Units','normalized', 'Position',[0.02,0.10,0.90,0.80]);
for ijk = 1:nY
    hax(ijk) = subplot(nY,1, ijk);
    plot(TimeSamples, disturbances(ijk,:), 'LineWidth', 1.5);
    grid on; 
    xlabel('time [s]'); ylabel(['d_', num2str(ijk), ' [mm]']);
end
linkaxes(hax, 'x');
% ---
figure('Name','Output Disturbance Reduction using ADRMPC', ...
    'Units','normalized', 'Position',[0.02,0.10,0.90,0.80]);
for ijk = 1:nY
    hax(ijk) = subplot(nY,1, ijk);
    plot(TimeSamples(Nstart:end), disturbances(ijk,(Nstart:end)), 'LineWidth', 1.5);
    hold on
    plot(TimeSamples(Nstart:end), ErrY(ijk,:), 'LineWidth', 1.5);
    grid on; 
    xlabel('time [s]'); ylabel(['d_', num2str(ijk), ' [mm]']);
end
linkaxes(hax, 'x');
% ---

Fs = f_camp_ESO; % 10 kHz sampling frequency
MAXFreq = 500; % Hz - the max frequency component to appear in the graphics
m2Mu_m = 1e+6; % conversion factor: from meters to micrometers
n_fft = simTime*Fs;

figure('Name', 'Output Disturbance Spectra - no ADRMPC', ...
    'Units','normalized', 'Position',[0.03,0.10,0.94,0.80]);
for abc=1:nY
    dist1 = disturbances(abc,:);
    % Take fourier transform
    fftSignal = fft(dist1);
    % apply fftshift to put it in the form we are used to (see documentation)
    fftSignal = fftshift(fftSignal);
    MagSig = (m2Mu_m)*((abs(fftSignal)*(1/Fs)))/n_fft; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    f = (-n_fft/2:n_fft/2-1)*(Fs/n_fft);
    subplot(nY,1, abc);
    plot(f(floor(n_fft/2)+1:end), MagSig(floor(n_fft/2)+1:end),...
        'LineWidth',1.5);
    xlim([0, MAXFreq]);
    title('Power Spectrum of Disturbance $d_1(t)$', 'Interpreter','latex');
    xlabel('Frequency (Hz)', 'Interpreter','latex');
    ylabel('Magnitude [$\mu m$]', 'Interpreter','latex');
end

figure('Name', 'Output Disturbance Spectra Reduction',...
    'Units','normalized', 'Position',[0.03,0.10,0.94,0.80]);
for abc=1:nY
    dist1 = disturbances(abc,:);
    % Take fourier transform
    fftSignal = fft(dist1);
    % apply fftshift to put it in the form we are used to (see documentation)
    fftSignal = fftshift(fftSignal);
    MagSig = (m2Mu_m)*((abs(fftSignal)*(1/Fs)))/n_fft; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    f = (-n_fft/2:n_fft/2-1)*(Fs/n_fft);
    subplot(nY,1, abc);
    plot(f(floor(n_fft/2)+1:end), MagSig(floor(n_fft/2)+1:end),...
        'LineWidth',1.5);
    hold on
    fftY = fft(ErrY(abc,:));
    numY = numel(ErrY(abc,:));
    fftSignalY = fftshift(fftY);
    MagSigY = (m2Mu_m)*((abs(fftSignalY)*(1/Fs)))/numY; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    fY = (-numY/2:numY/2-1)*(Fs/numY);
    plot(fY(floor(numY/2)+1:end), MagSigY(floor(numY/2)+1:end),...
        'LineWidth',2.5);
    xlim([0, MAXFreq]);
    title('Magnitude Spectrum of Disturbance $d_1(t)$', 'Interpreter','latex');
    xlabel('Frequency (Hz)', 'Interpreter','latex');
    ylabel('Magnitude [$\mu m$]', 'Interpreter','latex');
end