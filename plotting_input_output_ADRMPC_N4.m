clear all%#ok
close all
clc

rng(481516);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting the inputs & outputs 
% of the system controlled by 
% ADRMPC with horizon N = 4
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
Tskip = 1; % plot data every Tskip samples
m = 15; % number of inputs
p = 7 ; % number of outputs

% ---- simulation settings ----
        
% --- ---
% loading control inputs & controlled outputs
load ADRC_MPC_7Yx15U_10khz_N4.mat U Y 
% the output references
load goldenORBITS.mat
% output traking errors, initial sample ID & disturbances 
load RESULTS_MPC_ADRC_10kHz_N4.mat ErrY Nstart disturbances 
% --- ---


Y_bounds = [-2, +2];
U_bounds = [-10, +10];
p = 7; % # of outputs
m = 15; % # of inputs

Yref = repmat(H_Orbits_REFS, [1, nsample]);
% the orbit references
Y_lowerB = repmat(min(Y_bounds),[1, nsample]);
Y_upperB = repmat(max(Y_bounds),[1, nsample]);

U_lowerB = repmat(min(U_bounds),[1, nsample]);
U_upperB = repmat(max(U_bounds),[1, nsample]);



%% plotting inputs & outputs
nY = p;
hfYU = figure(); %('Name','Controlled Outputs', ...
    % 'Units','normalized', 'Position',[0.02,0.10,0.50,0.50]);

hax(1) = subplot(2,1, 1);
    for ijk = 1:p
        plot(TimeSamples(1:Tskip:end), Y(ijk,(1:Tskip:end)), 'LineWidth', 1.0);
        hold on; 
    end
    for ijk = 1:p
        plot(TimeSamples(1:500*Tskip:end), Yref(ijk,(1:500*Tskip:end)), 'LineWidth', 0.5, ...
            'LineStyle','-', 'Marker','diamond','MarkerSize',2,...
            'MarkerFaceColor','k','MarkerEdgeColor', 'k');
        hold on; 
    end
%     plot(TimeSamples(1:Tskip:end), Y_lowerB(1:Tskip:end), 'LineWidth', 1.5, ...
%             'LineStyle',':', 'Color', 'r');
%     plot(TimeSamples(1:Tskip:end), Y_upperB(1:Tskip:end), 'LineWidth', 1.5, ...
%             'LineStyle',':', 'Color', 'r');
    grid on
    %xlabel('Time [s]'); 
    ylabel('Output [mm]');
    title('Outputs, References and Contraints');

    ylim([-0.55, +0.55]);
% ---
hax(2) = subplot(2,1, 2);

    for ijk = 1:m
        plot(TimeSamples(1:Tskip:end), U(ijk,(1:Tskip:end)), 'LineWidth', 1.0);
        hold on; 
    end
    plot(TimeSamples(1:Tskip:end), U_lowerB(1:Tskip:end), 'LineWidth', 1.5, ...
            'LineStyle',':', 'Color', 'r');
    plot(TimeSamples(1:Tskip:end), U_upperB(1:Tskip:end), 'LineWidth', 1.5, ...
            'LineStyle',':', 'Color', 'r');
    grid on
    xlabel('Time [s]'); ylabel('Input [A]');
    ylim([min(U_bounds)-2, max(U_bounds)+2]);
    title('Inputs and Contraints');

linkaxes(hax, 'x');

%% saving the figure as TIKZ

% Matlab2Tikz toolbox needed!

% set(hfYU, 'Units', 'pixels');
% cleanfigure('handle', hfYU);
% matlab2tikz('figurehandle', hfYU,'filename','ADRMPC_sim.tex',...
%     'standalone', true, 'width','10cm');
%% disturbances spectra
Fs = 1e4; % 10 kHz sampling frequency
MAXFreq = 550; % Hz - the max frequency component to appear in the graphics
m2Mu_m = 1e+6; % conversion factor: from meters to micrometers
% % figure; pspectrum(dist1, Fs)
n_fft = simTime*Fs;

hfFFTALL = figure('Units','normalized');%, 'Position',[0.03,0.10,0.94,0.80]);
for abc=1:p
    dist1 = disturbances(abc,:);
    % Take fourier transform
    fftSignal = fft(dist1);
    % apply fftshift to put it in the form we are used to (see documentation)
    fftSignal = fftshift(fftSignal);
    MagSig = (m2Mu_m)*((abs(fftSignal)*(1/Fs)))/n_fft; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    f = (-n_fft/2:n_fft/2-1)*(Fs/n_fft);
    subplot(p,1, abc);
    plot(f(floor(n_fft/2)+1:end), MagSig(floor(n_fft/2)+1:end),...
        'LineWidth',0.8);
    xlim([0, MAXFreq]);
    if (abc==1)
        title('Magnitude Spectrum of Disturbances $d_{i}(t), i=1, 2, \ldots 7$', 'Interpreter','latex');
    end
    if (abc==p)    
        xlabel('Frequency (Hz)', 'Interpreter','latex');
    end
    ylabel('[$\mu m$]', 'Interpreter','latex');
end

%% saving the figure as TIKZ

% Matlab2Tikz toolbox needed!

% set(hfFFTALL, 'Units', 'pixels');
% cleanfigure('handle', hfFFTALL);
% matlab2tikz( 'figurehandle', hfFFTALL,'filename','FFT_ALL.tex',...
%     'standalone', true, 'width','10cm', 'height', '14cm');
%% only the spectrum of selected disturb
selD = 1; % select the additive disturbance acting on output selID

hfspectr=figure('Units','normalized');
for abc=selD
    dist1 = disturbances(abc,:);
    % Take fourier transform
    fftSignal = fft(dist1);
    % apply fftshift to put it in the form we are used to (see documentation)
    fftSignal = fftshift(fftSignal);
    MagSig = (m2Mu_m)*((abs(fftSignal)*(1/Fs)))/n_fft; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    f = (-n_fft/2:n_fft/2-1)*(Fs/n_fft);
    plot(f(floor(n_fft/2)+1:end), MagSig(floor(n_fft/2)+1:end),...
        'LineWidth',1.5);
    xlim([0, MAXFreq]);
    title(['Spectrum of Disturbance $d_', num2str(abc), '(t)$'], 'Interpreter','latex');
    xlabel('Frequency (Hz)', 'Interpreter','latex');
    ylabel('Magnitude [$\mu m$]', 'Interpreter','latex');
end
%% saving the figure as TIKZ

% Matlab2Tikz toolbox needed!

% set(hfspectr, 'Units', 'pixels');
% cleanfigure('handle', hfspectr);
% matlab2tikz('figurehandle', hfspectr,'filename','FFT.tex',...
%     'standalone', true, 'width','10cm');

%% spectra comparison - pointing out the amplitude reduction    
% ---
selD = 1; % select output & corresponding additive disturbance

hfspectr2=figure('Units','normalized');
for abc=selD
    dist1 = disturbances(abc,:);
    % Take fourier transform
    fftSignal = fft(dist1);
    % apply fftshift to put it in the form we are used to (see documentation)
    fftSignal = fftshift(fftSignal);
    MagSig = (m2Mu_m)*((abs(fftSignal)*(1/Fs)))/n_fft; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    f = (-n_fft/2:n_fft/2-1)*(Fs/n_fft);
    plot(f(floor(n_fft/2)+1:end), MagSig(floor(n_fft/2)+1:end),...
        'LineWidth',0.75);
    hold on
    fftY = fft(ErrY(abc,:));
    numY = numel(ErrY(abc,:));
    fftSignalY = fftshift(fftY);
    MagSigY = (m2Mu_m)*((abs(fftSignalY)*(1/Fs)))/numY; 
    % ((abs(fftSignal)*(1/Fs)).^2)/n_fft;
    % Next, calculate the frequency axis, which is defined by the sampling rate
    fY = (-numY/2:numY/2-1)*(Fs/numY);
    plot(fY(floor(numY/2)+1:end), MagSigY(floor(numY/2)+1:end),...
        'LineWidth',1.25);
    xlim([0, MAXFreq]);
    title('Magnitude Spectrum', 'Interpreter','latex');
    xlabel('Frequency (Hz)', 'Interpreter','latex');
    ylabel('Magnitude [$\mu m$]', 'Interpreter','latex');
    legend(['Disturb $d_', num2str(abc), '(t)$'], ...
        ['Tracking Error $y_{', num2str(abc),...
        '}(t)-y_{ref ',num2str(abc),' }(t)$'],...
         'Interpreter','latex');
end

%% saving the figure as TIKZ

% Matlab2Tikz toolbox needed!

% set(hfspectr2, 'Units', 'pixels');
% cleanfigure('handle', hfspectr2);
% matlab2tikz('figurehandle',hfspectr2, 'filename','ADRMPC_dred.tex',...
%     'standalone', true, 'width','10cm');