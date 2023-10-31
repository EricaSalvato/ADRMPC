function [Ac, Bc, Cc, Dc, ...
           p, n, uID, yID, ...
           umin, umax, ymax, ymin] = CT_10kHz_Elettra_sys_create(m)
% 
% m <--> n. of inputs
% -------
% 
    load RmH95X82.mat RmH
    % loading the (horizontal) response matrix RmH
    [maxY, mMAX] = size(RmH);
    % mMAX <--> max number of inputs
    pMAX = round(mMAX/2)+1; 
    % pMAX <--> max admissible number of outputs

    % -------------------- simulated plant ----
    %

    %{
    % --- NOT USED ---
    % % autovalori per le dinamiche lente 
    % % (banda 10 Hz <-- email di G. Gaio 26/05/2020)
    f3dB_slow = 45; % Hz
    om3dB_slow = 2 * pi * f3dB_slow;
    slow_lambda = -om3dB_slow;
    % autovalori per le dinamiche veloci 
    % (banda 500Hz <-- email di G. Gaio 26/05/2020)
    f3dB_actuator = 70; % Hz --> the effective pass-band for the ELETTRA
                        % actually used beam position correction magnets
    % ref. --> 
    % M. Lonza, D. Bulfone, V. Forchi and G. Gaio,
    % "Commissioning of the Elettra fast orbit feedback system,"
    % 2007 IEEE Particle Accelerator Conference (PAC), 
    % Albuquerque, NM, USA, 2007, pp. 203-205
    % --- NOT USED ---
    % % FdT = 1/(1+s/om) 
    %}

    % --- actuators dynamics: 1st order dynamics ---
    %       FdT = 1/(1+s/om)
    %
    f3dB_actuator = 500; % Hz --> the pass-band for the beam 
                         % position correction magnets
    % ref. --> 
    % Gayadeen, S., Heron, M.T., & Rehm, G. (2016). 
    % A unified approach to the design of orbit feedback 
    % with fast and slow correctors. 
    % In Riches, Kathleen (Ed.). Proceedings of the 15th International 
    % Conference on Accelerator and Large Experimental 
    % Physics Control Systems ICALEPCS 2015, (p. 1225). Australia 
    % <--
    om3dB_act = 2 * pi * f3dB_actuator;
    act_lambda = -om3dB_act;

    % ---------------------------------------------------------------
    n = m; % # of state vars: as many as the inputs
    % each corrector magnet can be described by a 1st order dynamics
    % ---------------------------------------------------------------

    %{
    % NOT USED
    % ----------------------------------------------------------
    % Questi sono gli indici dei magneti correttori lenti sulla
    % macchina reale
    
    slow_corr82 = [1 6 7 8 13 14 15 20 21 22 27 28 29 34 35 36 40 41 42 43 ...
        47 48 49 50 54 55 56 57 61 62 63 64 68 69 70 71 75 76 77 78 81 82];
    index4selection = find(slow_corr82 <= m, 1, 'last');
    slow_corr = slow_corr82(1:index4selection);
    p=length(slow_corr);
    % NOT USED
    %}

    % -------------------------------
    % available selection approaches:
    selectionACT = 1; % see the code below

    switch selectionACT 

    % a) a list of "selected" outputs (pMAX outputs)
        case 1    
            selectY = [1, 6, 7, 8, 13, 14, 15, 20, 21, 22, 27, 28, 29,...
                       34, 35, 36, 40, 41, 42, 43, 47, 48, 49, 50, 54, ...
                       55, 56, 57, 61, 62, 63, 64, 68, 69, 70, 71, 75,...
                       76 77 78 81 82];
    % ---
    % b) randomly choice of pMAX outputs from the available maxY in RmH
        case 2
            selectY = randi(maxY, [pMAX, 1]);
    % ---
    % c) a set of pMAX outputs, starting from the 1st one
        case 3   
            selectY = (1:pMAX);
    % ---
    % d) a set of pMAX outputs, backwards from the last one
        case 4
            selectY = (maxY-pMAX-1:maxY);
        otherwise
            selectY = [];
    end % switch-case 

    index4selection = find(selectY <= m, 1, 'last');
    selectY = selectY(1:index4selection);
    p=length(selectY);

    % --- select the output candidates ---
    ORtest = RmH(selectY, :); % it was (probably wrong) RmH(:, selectY);
    
    % --- now select the most effective inputs, given the outputs ---
    addpath(genpath('rga'));
    %     Relative Gain Array
    %     For nonsquare systems, the general rga can also be used 
    %     to select inputs if the number of inputs is larger than 
    %     the number of outputs. The second output in this case 
    %     gives the input effectiveness, E.
    %  
    %     [R,E] = rga(G)
    %  
    %     Where E is a column vector corresponding to the input. Inputs
    %     with largest input effectiveness are recommented to be
    %     selected for the control system.
    %
    [~, EE] = rga(ORtest);
    th = sort(EE); % order the input effectiveness, 
                   % from the least effective (the 1st position in EE)
                   % to the most effective input (the last position in EE)
    % th_value = th(end-p+1);    % <-- WRONG!!
    % KK = find(EE >= th_value); % <-- WRONG!!
    threshold_value = th(end-m+1); % <-- select the m most effective inputs 
    KK = find(EE >= threshold_value);% <-- the indices of the 
                                     %     most effective inputs
    ORtest = RmH(selectY, KK);
    % ORtest = RmH(KK, 1:m);
    ORMat = ORtest;
    uID = KK;
    yID = selectY;
    rmpath(genpath('rga'));
    % ----------------------------------------------------

    % -------------------------------------------------
    % state eqs. - continuous time system: corrector magnets+ELETTRA+BPM
    Adiag = ones(n,1) * act_lambda; % the 1st order dynamics of each
                                    % corrector magnet
    % -------------------------------------------------
    % MIMO continuous time system
    Ac = diag(Adiag);
    Bc = -Ac;
    Cc = ORMat;
    Dc = zeros(p, m);
    % MIMO continuous time system
    % -------------------------------------------------

    % --- Plant constraints ---
    Ucon = 10; %50; % [A] max input current for each corrector magnet
    
    umin = -Ucon*ones(m,1); % input constraints
    umax = +Ucon*ones(m,1); % input constraints

    ymax = +2; % output constraints [mm]
    ymin = -2; % output constraints [mm]
    % --- Plant constraints ---
end