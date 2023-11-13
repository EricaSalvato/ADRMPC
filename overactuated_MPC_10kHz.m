clear all%#ok
close all
clc

N_values = [4,8,12,16,20]; % Control and prediction horizon

for n_index=1:length(N_values)
    clearvars -except N_values n_index
    rng(481516);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NB both MPC & ADRC running @10kHz!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLANT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % --- operating frequencies ---
    f_camp_ESO = 1e4;          % ESO and plant operating frequency [Hz]
    T_camp_ESO = 1/f_camp_ESO; % ESO sampling period [s]
    Tsampling = T_camp_ESO;   % plant sampling period [s]
    % --- operating frequencies ---

    % ---- simulation settings ----
    simTime = 5; % time horizon: [0 simTime]
    nsample = floor(simTime/Tsampling);
    TimeSamples = linspace(0, simTime, nsample);
    % ---- simulation settings ----
    
    % ------- continuous-time model -------
     m = 15; % number of inputs
     p = 7 ; % number of outputs
    [Ac, Bc, Cc, Dc, n,... 
        umin, umax, ymax, ymin] = CT_sys_create(m, p); 
    % Ac, Bc, Cc, Dc <--> continuous-time system matrices
    % p <--> number of outputs
    % n <--> # of state variables (equal to # of inputs)
    % umin, umax <--> range of admissible input values
    % ymin, ymax <--> range of admissible  values for the outputs
    %
    sysC = ss(Ac, Bc, Cc, Dc);

    F = Cc;

    numY = p;
    numU = m;
    numX = n;

    % ------ continuous-time model ------

    
    % ---- sampled-time model with no disturbance ----
    E2BOCS_td = c2d(sysC, Tsampling, 'zoh' );
    Ad=E2BOCS_td.A;
    Bd=E2BOCS_td.B;
    Cd = E2BOCS_td.C;
    Dd = E2BOCS_td.D;
    % -------------------------------------------
    % -- compute the steady-state state, input & output ---
    mE = [Ad-eye(size(Ad)), Bd; ...
           Cd, Dd];
    mF = [mE, [zeros(n,p); eye(p)]];
    
    assert(rank(mE)==rank(mF));
    % if OK then the steady state equation
    %
    %     [A-I  B] [ x_s ]   [   0  ]
    %     [      ] [     ] = [      ]
    %     [ C   D] [ u_s ]   [ y_ref]
    %
    % has a solution FOR ANY SET-POINT y_ref
    %
    % ref: D. Limon, I. Alvarado, T. Alamo, E.F. Camacho,
    % "MPC for tracking piecewise constant references for 
    % constrained linear systems," Automatica, Volume 44, 
    % Issue 9, 2008, Pages 2382-2387
    %
    
    
    % -------------- reachability ----------------
    [ABAR,BBAR,CBAR,T,K] = ctrbf(Ad,Bd,Cd);
    assert(sum(K)==size(Ad,1)); 
    % (Ad, Bd) stabilizable <--> see 
    % G. Pannocchia, "Offset-free tracking MPC: A tutorial review 
    % and comparison of different formulations," 
    % 2015 European Control Conference (ECC), 
    % Linz, Austria, 2015, pp. 527-532, 
    % -------------------------------------------
    detectability_rank = rank(obsv(Ad, Cd));
     % Note: the obsv rank is < n 
     % BUT the system is asympt. stable, so it's
     % detectable!
     % see Pannocchia 2015
    check_rank = rank([Ad-eye(size(Ad)), Bd;...
          Cd, zeros(numY, numU)]);
    if (check_rank == numX+numY)
        disp('Assumption 1 in (Pannocchia, 2015) fulfilled!');
    else
        disp('The assumptions in (Pannocchia, 2015) are NOT fulfilled.')
    end
    
    %------- input & output constraints ----------
    maxU = umax(1); % [A]
    maxX = maxU;
    minU = umin(1); % [A]
    maxY = ymax(1); % [mm]
    minY = ymin(1); % [mm]
    
    % -----   reference orbits -----
    % golden orbit [mm]
    load("goldenORBITS.mat")
    Y_REF = H_Orbits_REFS;  % output references

    U_ref_effective = pinv(F)*Y_REF;
    X_ref_effective = U_ref_effective;
    % check the steady state solution computed exploiting pinv(F)
    errSS = mE * [X_ref_effective;U_ref_effective];
    checkXU = errSS(1:n,1);
    checkY = errSS(n+1:end,1);
    %     norm(checkXU)
    %     norm(checkY-Y_REF)
    %     max(U_ref_vere)
    
    % -- 2nd approach --
    % using linsolve --> mE * XX = YY
    YY = [zeros(n,1); Y_REF];
    XX = linsolve(mE, YY);
    dummyRES = mE*XX;
    Yss = dummyRES(n+1:end,1);
    XUs = dummyRES(1:n,1);
    %     norm(XUs)
    %     norm(Yss-Y_REF)
    X_ref_effective = XX(1:n,1);
    U_ref_effective = XX(n+1:end,1);
    %     max(U_ref_vere)
    OrbitH_REF = Y_REF;
    % Note: better approach: 
    % smaller error norm, inputs far from the
    % constraints
    
    % ---- sampled time model with augmented outputs ---
    % -------- the system looks like: 
    %
    % dx = x - xref
    % du = u - uref
    % dy = y - yref
    %
    % with y --> [y y_aug] 
    %
    %               \dot{dx} = A dx+ B du
    %                     dy = C dx
    %
    
    % ---- MPC weights ----
    R_weights = 1e-4; % Control input cost weights
    R = R_weights * eye(numU, numU);
    Q_weights = 1e3;  % State cost weights
    Q = Q_weights * eye(numY,numY);
    
    Wy = diag(10+rand(n,1));
    Qp = Wy'*Wy; % Wy = chol(Qy)
    % ---- MPC weights ----

    % --- the augmented system ---
    Daug = zeros(numY+numU, numU);
    
    E2BOCS_td_augOut = ss(Ad, Bd, [Cd; Wy], Daug, Tsampling);
    Ad_aug = E2BOCS_td_augOut.A;
    Bd_aug = E2BOCS_td_augOut.B;
    Cd_aug = E2BOCS_td_augOut.C;
    % -------------------------------------------
    
    
    
    %---------- MPC ----------
    simulation_time = 5; % 5 sec
    simulation_steps = simulation_time/Tsampling;
    N = N_values(n_index);
    Pred_H = N;
    Control_H = Pred_H;
    
    % --- using the MPC toolbox -----
    Weights.ManipulatedVariables = diag(R)';
    Weights.ManipulatedVariablesRate = [];
    Weights.OutputVariables = [diag(Q)', zeros(1,numU)];
    
    for i=1:numU
        ManipulatedVariablesCharacteristics(i).Min = minU-U_ref_effective(i);
        ManipulatedVariablesCharacteristics(i).Max = maxU-U_ref_effective(i);
    end
    
    for i=1:numY
        OutputVariablesCharacteristics(i).Min = minY-OrbitH_REF(i);
        OutputVariablesCharacteristics(i).Max = maxY-OrbitH_REF(i);
        OutputVariablesCharacteristics(i).MinECR = 0;
        OutputVariablesCharacteristics(i).MaxECR = 0;
    end

    % we need yaug just to insert the terminal cost
    setmpcsignals(E2BOCS_td_augOut, 'MO', (1:numY), 'UO', (numY+1:numY+numU)); 
    MPCobj = mpc(E2BOCS_td_augOut, Tsampling, Pred_H, Control_H, Weights, ...
        ManipulatedVariablesCharacteristics, OutputVariablesCharacteristics);
    
    % random inital states
    dummy_Y = 2*rand(numY,1)-1; % mm
    U_init = pinv(F)*dummy_Y; % corresponding inputs
    dummy_Y_acceptable4range = maxU/max(U_init)*dummy_Y; 
    % check the constraints and clip, if needed
    X_init = pinv(F)*dummy_Y_acceptable4range; % the initial state vector
    
    DeltaX_init = X_init - X_ref_effective;
    
    
    
    % let's generate the augmented outputs references
    % Note: these refs are zero as the aug. outputs are deltay
    DeltaY_ref_amp_aug = zeros(numY+numU,1);
    DeltaY_ref_amp_aug_signal = DeltaY_ref_amp_aug.*ones(numY+numU,...
        simulation_steps);
    
    
    % introducing the terminal cost
    TCWeight = 10000; % 1000
    Y = struct('Weight',[zeros(1,numY), TCWeight*ones(1,numU)]);
    U.Weights =[];
    setterminal(MPCobj,Y,U,Pred_H);
    
    % --- additive output disturbances ----
    disturbance_signals = DT_disturb_create(p, nsample); 
    % --- additive output disturbances ----    
            
    xref = X_ref_effective;
    uref = U_ref_effective;
    yref = OrbitH_REF;

    DeltaX = [];
    DeltaY = [];
    DeltaU = [];
    X = [];
    Y = [];
    U = [];
               
    deltax = DeltaX_init;
    deltay = Cd_aug*deltax;   
    x = deltax + xref;
    y = deltay(1:numY) + yref;
    xmpc = mpcstate(MPCobj);
    xmpc.Plant = deltax;

    beat_count =0;beatN=1000;
    tic
    for t = 1:simulation_steps

        beat_count = beat_count+1;
        if (beat_count == beatN)
            fprintf(1, '%d ...\n\r', t); drawnow;
            beat_count =0;
        end
    
        % Store the plant state.
        DeltaX = [DeltaX; deltax];
        X = [X; x];
        DeltaY = [DeltaY; deltay'];
        Y = [Y; y'];
          
        % Compute the MPC control action.
        [deltau, Info] = mpcmove(MPCobj,xmpc,deltay,DeltaY_ref_amp_aug);
        u = deltau + uref;
        DeltaU = [DeltaU; deltau'];
        U = [U; u'];
    
    
        % --- Update the state and the output of the plant of the plant ---
        % ------------  in the real plant we do not use it  ---------------
        % -----------  we  will receive the real value of y  --------------
        x = Ad*x + Bd*u;
        
        y = Cd*x+disturbance_signals(:,t); %il disturbo è campionato 10kHz
            
        % -----------------------------------------------------------------
        deltay = [y-yref; x-xref];
    end % for 
    
    elapsed_time_period = toc;
    
    comutation_time_controlinput = elapsed_time_period/simulation_steps;

    % plot results
    yref_sig = yref.*ones(numY,simulation_steps);   

    filename = ['MPC_7Yx15U_10khz_N',num2str(N),'.mat'];
    save(filename, 'TimeSamples', 'X', 'Y', 'U', 'yref', ...
        'elapsed_time_period')
    
    Nstart = 1000;
    ErrY = DeltaY(Nstart:end,1:p);
    errYSTD = std(ErrY,0,1)';
    errYVAR = var(ErrY,0,1)';
    maxErrY = max(ErrY,[],1)';
    disturbSTD = std(disturbance_signals,0,2);
    disturbVAR = var(disturbance_signals,0,2);
    maxDisturb = max(disturbance_signals,[],2);
    STD_reduction = abs(errYSTD ./ disturbSTD);
    VAR_reduction = abs(errYVAR ./ disturbVAR);
    maxReduction = abs(maxErrY ./ maxDisturb);
    errYnorm1 = sum(abs(DeltaY(Nstart:end,:)),1);
    
    filename = ['RESULTS_MPC_10kHz_N',num2str(N),'.mat'];
    save(filename, 'Nstart', 'ErrY', 'errYSTD', 'errYVAR', 'maxErrY','disturbSTD',...
        'disturbVAR', 'maxDisturb', 'STD_reduction', 'VAR_reduction', ...
        'maxReduction', 'errYnorm1','disturbance_signals')
   
end
