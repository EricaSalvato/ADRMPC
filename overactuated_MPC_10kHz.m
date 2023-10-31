clear all%#ok
close all
clc

N_values = [4,8,12,16,20];%36; % Control and prediction horizon

for n_index=1:length(N_values)
    clearvars -except N_values n_index
    rng(481516);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NB both MPC & ADRC running @10kHz!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% PLANT
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    m = 15; % number of inputs
    f_camp = 1e4; % ESO and plant operating frequency [Hz]
    simTime = 5; % time horizon: [0 simTime]
            
    % operating frequencies
    T_camp_sys = 1/f_camp; % sampling period [s]
    Tsampling = T_camp_sys;
    % simulation settings
    nsample = floor(simTime/T_camp_sys);
    TimeSamples = linspace(0, simTime, nsample);
    
    % --- continuous-time model ---
    [Ac, Bc, Cc, Dc, p, n, U_indices, ...
        Y_indices, umin, umax, ymax, ymin] = CT_10kHz_Elettra_sys_create(m); 
    % Ac, Bc, Cc, Dc <--> continuous time system matrices
    % p <--> number of outputs
    % n <--> # of state variables
    % U_indices <--> 
    % Y_indices <-->
    % umin, umax <-->
    % ymin, ymax <-->
    
    analisi_svd = false;
    
    if ~analisi_svd 
        F = Cc;
    else
        S = svd(Cc); % i valori singolari in un vettore
        % analisi con SVD <-- 7 valori singolari NON nulli
        %
        % S =
        %       4.151462447642301
        %       2.090673771233182
        %       0.563358832018362
        %       0.464391717916251
        %       0.282887875974029
        %       0.160685076694029
        %       0.027072355063732
        %
        [U_m, Fsvd, V] = svd(Cc);
        howManySVal = numel(S);
        FkSVD = cell(howManySVal,1); 
        sVal_ijk_contribMAT = cell(howManySVal,1); 
        for ijk = 1:howManySVal
            U_ColVect = U_m(:,ijk);
            sVal = S(ijk);
            VT_RowVect = V(:,ijk)';
            sVal_ijk_contribMAT{ijk} = sVal * U_ColVect * VT_RowVect;
        end % for ijk
        
        FkSVD{1} = sVal_ijk_contribMAT{1};
        for abc=2:howManySVal
            FkSVD{abc} = FkSVD{abc-1} + sVal_ijk_contribMAT{abc};
        end % for abc
        Cc = FkSVD{7};
    end
    
    sysC = ss(Ac, Bc, Cc, Dc);
    
    numY = p;
    numU = m;
    numX = n;
    % --- continuous-time model ---
    
    % --- sampled-time model ---
    % sysD = c2d(sysC, T_camp_sys);
    % Ad = sysD.A;
    % Bd = sysD.B;
    % Cd = sysD.C;
    % Dd = sysD.D;
    % --- sampled-time model ---
    
    
    % ----------------- modello a tempo campionato senza disturbo--------------
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
    
    
    % -------------- controllabilita' a tempo campionato ----------------
    [ABAR,BBAR,CBAR,T,K] = ctrbf(Ad,Bd,Cd);
    assert(sum(K)==size(Ad,1)); 
    % (Ad, Bd) stabilizable <--> see 
    % G. Pannocchia, "Offset-free tracking MPC: A tutorial review 
    % and comparison of different formulations," 
    % 2015 European Control Conference (ECC), 
    % Linz, Austria, 2015, pp. 527-532, 
    % -------------------------------------------
    detectability_rank = rank(obsv(Ad, Cd));
     % Note: the obsv rank is < n BUT the system is asympt. stable, so it's
     % detectable!
     % see Pannocchia 2015
    check_rank = rank([Ad-eye(size(Ad)), Bd;...
          Cd, zeros(numY, numU)]);
    if (check_rank == numX+numY)
        disp('Assumption 1 in (Pannocchia, 2015) fulfilled!');
    else
        disp('The assumptions in (Pannocchia, 2015) are NOT fulfilled.')
    end
    
    %------- vincoli input output ----------
    maxU = umax(1); % [A]
    maxX = maxU;
    minU = umin(1);
    maxY = ymax(1);   % [mm]
    minY = ymin(1);
    
    % -----   Orbite reali da ottenere -----
    real_ref = true;
    if real_ref 
        File4Data = 'BPM_DATA_FILES/BPM_TREND_INFO_DATA.mat';
        load(File4Data);
        OrbitH_vere = zeros(numel(bpm_h_off_trendInfo),1);
        for i=1:numel(bpm_h_off_trendInfo)
            OrbitH_vere(i,1) = bpm_h_off_trendInfo(i).offset; % in micron!!!
        end
        
        OrbitH_vere = OrbitH_vere(Y_indices).*0.001; %Da micron a mm
        U_ref_vere = pinv(F)*OrbitH_vere;
        Y_REF = OrbitH_vere;
    else
        % OrbitH_vere = 2*rand(numY,1)-1;
        OrbitH_vere = 0.5*randn(p,1); % golden orbit [mm]
        Y_REF = OrbitH_vere;
        U_ref_vere = pinv(F)*OrbitH_vere;
        OrbitH_vere = (maxU/max(U_ref_vere)).*OrbitH_vere;
        U_ref_vere = pinv(F)*OrbitH_vere;
    end
    X_ref_vere = U_ref_vere;
    
    
    % check the steady state solution computed exploiting pinv(F)
    
    errSS = mE * [X_ref_vere;U_ref_vere];% - [zeros(n,1); Y_REF]
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
    X_ref_vere = XX(1:n,1);
    U_ref_vere = XX(n+1:end,1);
%     max(U_ref_vere)
    OrbitH_vere = Y_REF;
    
    % ---- modello a tempo campionato con uscite aumentate per terminal cost---
    % -------- stiamo considerando in questo caso un sistema del tipo ---------
    %
    % dx = x - xref
    % du = u - uref
    % dy = y - yref
    %
    % con y = [y y_aug] come nel manuale <-- ?? quale??
    %
    %               \dot{dx} = A dx+ B du
    %                     dy = C dx
    %
    % l'MPC lavora con questo.
    WeightQ = 1000;
    % Q = WeightQ * eye(numY,numY);
    % R = 0.0001 * eye(numU, numU);
    
    R_weights = 1e-4; %Control input cost weights
    R = R_weights * eye(numU, numU);
    Q_weights = 1e3;  % State cost weights
    Q = Q_weights * eye(numY,numY);
    
    Wy = diag(10+rand(n,1));
    Qp = Wy'*Wy; % Wy = chol(Qy)
    
    Daug = zeros(numY+numU, numU);
    
    E2BOCS_td_augOut = ss(Ad, Bd, [Cd; Wy], Daug, Tsampling   );
    Ad_aug = E2BOCS_td_augOut.A;
    Bd_aug = E2BOCS_td_augOut.B;
    Cd_aug = E2BOCS_td_augOut.C;
    % -------------------------------------------
    
    
    
    %---------- MPC ----------
    simulation_time = 5; %5 sec
    simulation_steps = simulation_time/Tsampling;
    N = N_values(n_index);
    Pred_H = N;
    Control_H = Pred_H;
    
    Weights.ManipulatedVariables = diag(R)';
    Weights.ManipulatedVariablesRate = [];
    Weights.OutputVariables = [diag(Q)', zeros(1,numU)];
    
    %ManipulatedVariablesCharacteristics(1:numU) = struct('Min', minU, 'Max',...
    % maxU, 'Units', "A"); % di default sono hard
    for i=1:numU
        ManipulatedVariablesCharacteristics(i).Min = minU-U_ref_vere(i);
        ManipulatedVariablesCharacteristics(i).Max = maxU-U_ref_vere(i);
    end
    
    %OutputVariablesCharacteristics(1:numY+numU) =  struct('Min', minY,...
    % 'Max', maxY, 'MinECR', 0, 'MaxECR', 0, 'Units', "mm");
    for i=1:numY
        OutputVariablesCharacteristics(i).Min = minY-OrbitH_vere(i);
        OutputVariablesCharacteristics(i).Max = maxY-OrbitH_vere(i);
        OutputVariablesCharacteristics(i).MinECR = 0;
        OutputVariablesCharacteristics(i).MaxECR = 0;
    end
    % settiamo yaug non osservabile poichè è un gioco che facciamo solo per 
    % inserire il terminal cost
    setmpcsignals(E2BOCS_td_augOut, 'MO', (1:numY), 'UO', (numY+1:numY+numU)); 
    MPCobj = mpc(E2BOCS_td_augOut, Tsampling, Pred_H, Control_H, Weights, ...
        ManipulatedVariablesCharacteristics, OutputVariablesCharacteristics);
    
    
    % genera stati inziali casuali.
    Y_di_servizio = 2*rand(numY,1)-1; % mm
    U_init = pinv(F)*Y_di_servizio;
    Y_di_servizio_acceptable4range = maxU/max(U_init)*Y_di_servizio;
    X_init = pinv(F)*Y_di_servizio_acceptable4range;
    
    DeltaX_init = X_init - X_ref_vere;
    
    
    
    % genera le uscite di riferimento aumentate che vanno a zero in quando sono
    % le deltay 
    DeltaY_ref_amp_aug = zeros(numY+numU,1);
    DeltaY_ref_amp_aug_signal = DeltaY_ref_amp_aug.*ones(numY+numU,...
        simulation_steps);
    
    
    % aggiungo il terminal cost
    TCWeight = 10000; % 1000
    Y = struct('Weight',[zeros(1,numY), TCWeight*ones(1,numU)]);
    U.Weights =[];
    setterminal(MPCobj,Y,U,Pred_H);
    
    real_disturbance = 1; 
    if real_disturbance == 0
        % Creo i disturbi come la somma di tre sinusoidi di ampiezza tra -1/2e 1/2
        % e di frequenza compresa tra i 10 e i 30Hz
        minF = 10;
        maxF = 30;
        ampd =1/2;
        df = maxF - minF;
        t = (0:Tsampling:Tsampling*simulation_steps-Tsampling)';
        disturbance_signal = (rand()-ampd)*sin(2*pi*(df*rand()+minF)*t) + ...
                    (rand()-ampd)*sin(2*pi*(df*rand()+minF)*t) + ...
                    (rand()-ampd)*sin(2*pi*(df*rand()+minF)*t);
        disturbance_signals = repmat(disturbance_signal',numY,1);
    else
        File4Data = 'ELETTRA2_BPM_OFF_ON_detrend_2021_21_29.mat';
        load(File4Data);
        disturbance_data = bpm_h_off_detrend(:,Y_indices)';
        % disturbance_data = cell2mat(BPM_H_disturbances);
        % disturbance_signals = reshape(disturbance_data,[numY,length(BPM_H_disturbances{1, 1})]);
        disturbance_signals = disturbance_data .* 0.001; 
        % converto micron in mm
    end
    
    simulate_run_time = 1; %0 se si vuole lavorare solo in simulazione
    
    if simulate_run_time == 0
        % simulazione completa sistema con i delta (riferimenti a zero)
        simOpt = mpcsimopt(MPCobj);
        simOpt.PlantInitialState = DeltaX_init;
        [Ympc,~,Umpc,~,~,Options,status] = sim(MPCobj,simulation_steps,...
            DeltaY_ref_amp_aug_signal',simOpt);
        
%         h1 = figure();
%         for plotId = 1 : numY
%             subplot(12, 8, plotId)
%             plot(Ympc(:,plotId)')
%             hold on 
%             plot(DeltaY_ref_amp_aug_signal(plotId,:))
%             grid on
%             name = sprintf("DeltaY_%d",plotId);
%             title(name);
%         end
%         namefig=sprintf("DeltaY_PH%d_CH%d", Pred_H, Control_H);
%         savefig(h1,namefig)
%         
%         h2 = figure();
%         for plotId = 1 : numU
%             subplot(12, 7, plotId)
%             plot(Umpc(:,plotId)')
%             grid on
%             name = sprintf("DeltaU_{%d}",plotId);
%             title(name);
%         end
%         namefig=sprintf("DeltaU_PH%d_CH%d", Pred_H, Control_H);
%         savefig(h2,namefig)
    
    else  
        
        xref = X_ref_vere;
        uref = U_ref_vere;
        yref = OrbitH_vere;
    
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
    %         Cost = Info.Cost;
            u = deltau + uref;
            DeltaU = [DeltaU; deltau'];
            U = [U; u'];
    
    
            % --- Update the state and the output of the plant of the plant ---
            % ------------  in the real plant we do not use it  ---------------
            % -----------  we  will receive the real value of y  --------------
            x = Ad*x + Bd*u;
            
            y = Cd*x+disturbance_signals(:,t); %il disturbo è campionato 10kHz
            
            % -----------------------------------------------------------------
            
    %         % Update the state of the model
    %         deltax = Ad_aug*deltax + Bd_aug*deltau;
    
            % deltay = Cd_aug*deltax;
            deltay = [y-yref; x-xref];
    %         deltay = [y-yref; deltax];
        end
        fine = toc;
        
        durata_comp_controlinput = fine/simulation_steps;
    
        % plot results
        yref_sig = yref.*ones(numY,simulation_steps);   
%         h3 = figure();
%             
%         max_error = zeros(1,numY);
%         RMSE = zeros(1,numY);
%         for plotId = 1 : numY
%             subplot(1, numY, plotId)
%             plot(Y(:,plotId)', 'Linewidth',1.5)
%             hold on 
%             plot(yref_sig(plotId,:),'Linewidth',2.5)
%             plot(minY*ones(1,simulation_steps),'Linewidth',1,'Color','r')
%             plot(maxY*ones(1,simulation_steps),'Linewidth',1,'Color','r')
%             grid on
%             ylim([minY-1 maxY+1])
%             name = sprintf("Correctors500Hz_Orbits_{%d}",plotId);
%             title(name,'Interpreter','none');
%             max_error(plotId) = max((Y(simulation_steps-2*Pred_H:simulation_steps,plotId)'-yref_sig(plotId, simulation_steps-2*Pred_H:simulation_steps))); 
%             RMSE(plotId) = sqrt(mean(Y(simulation_steps-2*Pred_H:simulation_steps,plotId)'-yref_sig(plotId, simulation_steps-2*Pred_H:simulation_steps)).^2);
%         end
%         namefig=sprintf("Correctors500Hz_Orbits_PH%d_CH%d_withDisturbancesOnY_TCWeights%d_YWeights%d", Pred_H, Control_H, TCWeight, WeightQ);
%         savefig(h3,namefig)
%         save(namefig,'Y')
%         h4 = figure();
%         for plotId = 1 : numY
%             subplot(1, numY, plotId)
%             plot(DeltaY(:,plotId)','Linewidth',1.5)
%             hold on 
%             plot(DeltaY_ref_amp_aug_signal(plotId,:),'Linewidth',2.5)
%             plot((minY-OrbitH_vere(plotId))*ones(1,simulation_steps),'Linewidth',1,'Color','r')
%             plot((maxY-OrbitH_vere(plotId))*ones(1,simulation_steps),'Linewidth',1,'Color','r')
%             grid on
%             ylim([minY-OrbitH_vere(plotId)-1 maxY-OrbitH_vere(plotId)+1])
%             name = sprintf("Correctors500Hz_DeltaY_{%d}",plotId);
%             title(name,'Interpreter','none');
%         end
%         namefig=sprintf("Correctors500Hz_DeltaY_PH%d_CH%d_simulationStepByStep_withDisturbancesOnY_TCWeights%d_YWeight%d", Pred_H,...
%             Control_H, TCWeight, WeightQ);
%         savefig(h4,namefig)
%     
%         h5 = figure();
%         for plotId = 1 : numU
%             subplot(5, 3, plotId)
%             plot(U(:,plotId)','Linewidth',1.5)
%             hold on
%             plot(minU*ones(1,simulation_steps),'Linewidth',1,'Color','r')
%             plot(maxU*ones(1,simulation_steps),'Linewidth',1,'Color','r')
%             grid on
%             name = sprintf("Correctors500Hz_Control_{%d}",plotId);
%             ylim([minU-1 maxU+1])
%             title(name,'Interpreter','none');
%         end
%         namefig=sprintf("Correctors500Hz_Control_PH%d_CH%d_withDisturbancesOnY_TCWeights%d_YWeights%d", Pred_H, Control_H, TCWeight, WeightQ);
%         savefig(h5,namefig)
%         save(namefig,'U')
    end
    filename = ['MPC_7Yx15U_10khz_N',num2str(N),'.mat'];
    save(filename, 'TimeSamples', 'X', 'Y', 'U', 'yref', ...
        'fine', 'U_indices', 'Y_indices')
    
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
