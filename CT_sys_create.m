function [Ac, Bc, Cc, Dc, ...
          n, ...
           umin, umax, ymax, ymin] = CT_sys_create(m, p)
% 
% m <--> n. of inputs
% p <--> n. of outputs
% -------
% Ac, Bc, Cc, Dc <--> continuous-time system matrices
% p <--> number of outputs
% n <--> # of state variables (equal to # of inputs)
% umin, umax <--> range of admissible input values
% ymin, ymax <--> range of admissible  values for the outputs

    load ORmatrix.mat ORMat
    % loading the response matrix ORmat
    [pMAX, mMAX] = size(ORMat);
    % mMAX <--> max number of inputs
    % pMAX <--> max admissible number of outputs
    p = min(p, pMAX);
    m = min(m, mMAX);

    % -------------- simulated plant ---------------
    

    % ------ Plant constraints ------
    Ucon = 10; % [A] max input current for each corrector magnet
    Ycon = 2; % output constraints [mm]
    % ------ Plant constraints ------

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

    % --------------------------------------
    % select inputs & outputs according to 
    % the values of m and p
    actualORmatrix = ORMat(1:p, 1:m);
    % --------------------------------------

    % -------------------------------------------------
    % state eqs. - continuous time system: 
    %              corrector magnets+synchrotron+BPM
    Adiag = ones(n,1) * act_lambda; % the 1st order dynamics of each
                                    % corrector magnet
    % -------------------------------------------------
    % MIMO continuous time system
    Ac = diag(Adiag);
    Bc = -Ac;
    Cc = actualORmatrix;
    Dc = zeros(p, m);
    % MIMO continuous time system
    % -------------------------------------------------

    % --- Plant constraints ---    
    umin = -Ucon*ones(m,1); % input constraints
    umax = +Ucon*ones(m,1); % input constraints

    ymax = +Ycon; % output constraints [mm]
    ymin = -Ycon; % output constraints [mm]
    % --- Plant constraints ---
end