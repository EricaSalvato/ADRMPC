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
    f_camp_ESO = 1e4; % ESO and plant operating frequency [Hz]
    simTime = 5; % time horizon: [0 simTime]
            
    % operating frequencies
    T_camp_ESO = 1/f_camp_ESO; % sampling period [s]
    T_camp_sys = T_camp_ESO; 
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
    sysC = ss(Ac, Bc, Cc, Dc);
    % --- continuous-time model ---
    
    % --- sampled-time model ---
    sysD = c2d(sysC, T_camp_sys);
    Ad = sysD.A;
    Bd = sysD.B;
    Cd = sysD.C;
    Dd = sysD.D;
    % --- sampled-time model ---
    
    % --- C*B pseudo-inverse ---
    % CBinverse useful also if by using it 
    % we respect constraints on the plant input (avoid to use the LP) 
    CBinv_prov = pinv(Cc*Bc);
    % --- C*B pseudo-inverse ---
    
    % --- disturbances - from REAL plant ---
    disturbances = DT_Elettra_dist_create(Y_indices, nsample); 
    errELETTRA = DT_Elettra_orbitERROR(Y_indices, nsample);
    % --- disturbances - from REAL plant ---
     
    % figure('Name','output disturbance');
    % plot(TimeSamples, Fdisturbance_mm, 'LineWidth', 1.5);
    % grid on; xlabel('time [s]'); ylabel('output disturbance [mm]');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESO 4 ADRC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold = 41;
            
    [Lp, Lc, Aeso, Bes, ...
        BesD, index_ESO_stable] = ESO_create(threshold, f_camp_ESO,...
                                                T_camp_ESO, p);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SoMPC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq_red_factor = 1; %10% frequency reduction factor of MPC with respect
    R_weights = 10; %Control input cost weights
    Q_weights = 10000;  % State cost weights

    N = N_values(n_index);
    fprintf(1, 'Horizon_value: %d\n', N);
    % --- continuous time system viewed by MPC ---
    [A_MPC, B_MPC, C_MPC] = TC_MPC_create(p);
    [~,c_B] = size(B_MPC);
    % --- continuous time cost matrices ---
    R = R_weights*eye(c_B);
    Q = Q_weights*eye(size(A_MPC)); 
    % --- Euler Auxiliary system ---
    % [A_EAS, B_EAS, C_EAS, tau, Q_EAS, R_EAS] = EAS_create(freq_red_factor,f_camp_ESO, A_MPC, B_MPC, C_MPC, Q, R);
    [A_EAS, B_EAS, C_EAS, tau, Q_EAS, R_EAS] = ...
        EAS10kHz_create(f_camp_ESO, A_MPC, B_MPC, C_MPC, Q, R);
    
    
    [r_B_EAS, c_B_EAS] = size(B_EAS);
    % --- N step reachability test ---
    [Reach_matrix, Obs_matrix] = reach_obsMatrices(N, A_EAS, B_EAS, C_EAS);
    assert(rank(Reach_matrix)==length(A_EAS),'Not reachable')
    assert(rank(Obs_matrix)==length(A_EAS),'Not observable')
    % --- terminal weigth matrix 4 stability ensurance ---
    % [P_EAS,~,~] = idare(A_EAS,B_EAS,1/2*Q_EAS,1/2*R_EAS,[],[]);  % Discrete time riccati equation Matlab 2021
    [P_EAS,~,~] = dare(A_EAS,B_EAS,1/2*Q_EAS,1/2*R_EAS,[],[]); % Discrete time riccati equation Matlab 2017
    % --- MPC ---
    [H, F, M, hatS, hatT, Xm_MPC, XM_MPC] = MPC_create(N, Reach_matrix, A_EAS, B_EAS, R_EAS, Q_EAS, P_EAS, ymax);
    % --- Quadprog settings 4 MPC ---
    H_quadprog = H;
    A_quadprog_no_uoConstr = [hatS;-hatS];
    Aeq_quadprog = hatS(end-length(B_EAS)+1:end,:);
    
    % --- storages for data savings
    U0 = []; % MPC storage
    Jcomp = 0; % MPC Cost value
    
    
    % ---- disturbance settings ----
    %Yref = 0.5*randn(p,1); % golden orbit [mm]
    load("goldenORBITS.mat")
    Yref = OrbitH_vere;
    YrefSEQ=[];%#ok
    YrefSEQ = repmat(Yref,1, nsample+N);
    
    %% Controlled system
    % ---- initialize ------
    x0=rand(n,1);
    
    xk = x0;
    z1hat = []; z2hat = [];
    for yn=1:p
        zhat0{yn} = zeros(2,1);%#ok
        z2hat0{yn} = zhat0{yn}(2);%#ok
        z1hat0{yn} = zhat0{yn}(1);%#ok
        
        z2hatk{yn} = z2hat0{yn};%#ok
        z1hatk{yn} = z1hat0{yn};%#ok
        zhatk{yn} = zhat0{yn};%#ok
        z1hat = [z1hat;z1hatk{yn}];%#ok
        z2hat = [z2hat;z2hatk{yn}];%#ok
    end
    
    X = [ ];
    Y = [];
    U = [];
    Zhat = [];
    H_old=[];
    u0_range_old = [];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Simulation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% ---- flag ----
    z2_bounds_updated_LP = 0;
    z2_bounds_updated_MPC = 0;
    
    %%% ---- linprog settings ----
    CB = Cc*Bc;
    %%% ---- {1,2}-inverse of CB computing ---
    CB_inv = left_pseudo_inverse12(CB);
    assert(isequal(any(CB*CB_inv-eye(p)>1.0e-12*ones(p)),boolean(zeros(1,p))),'Condition for linprog formulation not satisfied')
    
    %%% ---- for MPC input constraints ----
    [rCBinv, cCBinv] = size(CB_inv);
    CB_inv_rep = [CB_inv;zeros(N*rCBinv-rCBinv,cCBinv)];
    [rCBrep,cCBrep] = size(CB_inv_rep);
    for i=1:N-1
        new_col = [zeros(i*rCBinv,cCBinv);CB_inv;zeros(rCBrep-i*rCBinv-rCBinv,cCBinv)];
        CB_inv_rep = [CB_inv_rep, new_col];%#ok
    end
    %%% ---- (I-CB_inv*CB) SVD ----
    [Usvd, Ssvd, Vsvd] = svd(eye(m)-CB_inv*CB);
    Ssvd(Ssvd<0.00001) = 0; %********* 0.00001 è a capocchia possiamo trovare un valore idoneo?? ********* 
    Yfactor4linprog = Usvd*Ssvd*Vsvd';
    assert(rank(Yfactor4linprog) == m-p,'Condition for y_linprog not ok')
    [row_Yfactor4linprog,col_Yfactor4linprog] = size(Yfactor4linprog);
    
    %%% ---- ESO bounds init ----
    z2hat_max=1/10*ones(length(z2hat),1);
    z2hat_min=-1/10*ones(length(z2hat),1);
    %%% ---- Loop START ----
    fprintf(1, 'starting simulation loop:\n\r %d simulation steps \n\r' ,nsample )
    tic;
    beat_count = 0;
    beatN = max(floor(1/T_camp_ESO), ceil(1/T_camp_ESO))/100;
    
    for i=1:nsample
        beat_count = beat_count+1;
        if (beat_count == beatN)
            fprintf(1, '%d ...\n\r', i); drawnow;
            beat_count =0;
        end
        
        yk = Cd*xk+disturbances(:,i);
    
        if i > index_ESO_stable  %we wait ESO stabilization before starting the MPC
            %%% ---- z2 bounds updating check ----
            z2hat_max = max(z2hat_max,z2hat);
            z2hat_min = min(z2hat_min,z2hat);  
    
            if any((z2hat_max - z2hat)==0) || any((z2hat_min - z2hat)==0)
                z2_bounds_updated_MPC = 1;
                z2_bounds_updated_LP = 1;
                z2hat_mean = (z2hat_min+z2hat_max)/2;
                delta = (z2hat_max-z2hat_min)/2; 
                Mat_comb = dec2bin(0:2^p-1)' - '0';
                Mat_comb(Mat_comb==0) = -1 ;
                Mat_comb = delta.*Mat_comb;
                b_quadprog_upper=[];
                b_quadprog_lower=[];
                for iiii=1:rCBinv
                    b_quadprog_u_constraints_pos = min(umax(iiii)+CB_inv(iiii,:)*z2hat_mean+CB_inv(iiii,:)*Mat_comb);
                    b_quadprog_u_constraints_neg = min(-umin(iiii)-CB_inv(iiii,:)*z2hat_mean-CB_inv(iiii,:)*Mat_comb);
    
                    b_quadprog_upper = [b_quadprog_upper; b_quadprog_u_constraints_pos];%#ok 
                    b_quadprog_lower = [b_quadprog_lower; b_quadprog_u_constraints_neg];%#ok
                end
                A_quadprog = [A_quadprog_no_uoConstr;CB_inv_rep;-CB_inv_rep];
            end
            if rem(i-1,freq_red_factor) == 0 % the MPC works freq_red_factor times slower than ESO
                % ---- MPC start
                mpc_updated = 1;
                x0_MPC = z1hat;
                dummy = YrefSEQ(:,i:i+N-1);
    
                f_quadprog = F*x0_MPC-M*dummy(:); 
                b_quadprog_no_uoConstr = [XM_MPC-hatT*x0_MPC;...
                     -(Xm_MPC-hatT*x0_MPC)];
    
                if z2_bounds_updated_MPC
                    z2_bounds_updated_MPC = 0;
                    b_quadprog = [b_quadprog_no_uoConstr;repmat(b_quadprog_upper,N,1);repmat(b_quadprog_lower,N,1)];
                end
                beq_quadprog = YrefSEQ(:,i) - A_EAS^N*x0_MPC;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%   MANCA IL CHECK DELL'ORIZZONTE MPC  %%%%%%%%%%%%%%%
    %             N_feasibility=max(max(abs(ymin-YrefSEQ(:,i)),abs(ymax-YrefSEQ(:,i))))/min(min(abs(u0min),abs(u0max))); ----> ho CB^{(1)}*u_0<b ma come determino u0min e u0max???           
    %             error=sprintf('N_feasibility must be > %d',N_feasibility);
    %             assert(N>N_feasibility, error)
    %             
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % ---- MPC control sequence
                options_quadprog = optimoptions('quadprog','Display', 'off','MaxIterations',100);% ********
                [u_star,fval,exitflag] = quadprog(H_quadprog,f_quadprog,A_quadprog,b_quadprog,...
                    Aeq_quadprog,beq_quadprog,[],[],[],options_quadprog);
                assert(exitflag==1||exitflag==0,'MPC does not converge to a solution.')
                if exitflag==0
                    disp('exitflag=0')
                end
                % ---- only first input
                u0 = u_star(1:length(B_EAS),1); 
                Jcomp = Jcomp+1/2*x0_MPC'*Q_EAS*x0_MPC+1/2*u0'*R_EAS*u0;
                % ---- MPC end
                U0 = [U0, u0];%#ok
            end
            UK = (u0-z2hat);
       
            if z2_bounds_updated_LP || mpc_updated %Solve LP if Z2 bounds changed or MPC worked
                mpc_updated = 0;
                z2_bounds_updated_LP = 0;
                
                %%% ---- LP for u0-z2 -> u mapping with u constraints satisfaction
                f_linprog = [zeros(m,1);1];
                
                options_linprog = optimoptions('linprog','Display','off','MaxIterations',4000);
                
                A_linprog = [Yfactor4linprog,zeros(row_Yfactor4linprog,1);...
                    -Yfactor4linprog,zeros(row_Yfactor4linprog,1);...
                    eye(row_Yfactor4linprog),-ones(col_Yfactor4linprog,1);...
                    -eye(row_Yfactor4linprog),-ones(col_Yfactor4linprog,1)];
    
                b_linprog_upper = b_quadprog_upper-CB_inv*u0;
                b_linprog_lower = b_quadprog_lower-CB_inv*u0;
    
                b_linprog = [b_linprog_upper; b_linprog_lower;zeros(2*col_Yfactor4linprog,1)];
                [y_linprog,fval_linprog,exitflag_linprog] = linprog(f_linprog,A_linprog,b_linprog,[],[],[],[],options_linprog);
                assert(exitflag_linprog==1||exitflag_linprog==0,'linprog fault')
                if exitflag_linprog==0
                    disp('linprog exitflag=0')
                end
    
            end
            uk = CB_inv*UK+Yfactor4linprog*y_linprog(1:end-1);
            
            if any(abs(uk)>umax)
                keyboard;
            end
    
            
        else
            u0 = zeros(p,1);
            UK = zeros(p,1);
            uk= pinv(Cc*Bc)*UK; % in evoluzione libera???
            if any(abs(uk)>umax)
                keyboard; %#ok<*KEYBOARDFUN>
            end
        end
    
        %%% ---- one-step plant evolution
        xk_new = Ad*xk+Bd*uk;
        
        %%% ---- variables update
        for yn=1:p
            AesoMAT = Aeso{yn};
            BesMAT = Bes{yn};
            zhatkV = zhatk{yn};
            LcY = Lc{yn};
            LpY = Lp{yn};
            zhatk_new{yn} = AesoMAT*zhatkV + BesD{yn}*UK(yn) + LpY*yk(yn);%#ok
            
        end
        
        X = [X,xk];%#ok
        Y = [Y, yk];%#ok
        U = [U, uk];%#ok
        
        Zhat = [Zhat, [z1hat;z2hat]];%#ok
        for yn=1:p
            z1hatk{yn} = zhatk_new{yn}(1);
            z2hatk{yn} = zhatk_new{yn}(2);
            zhatk{yn} = zhatk_new{yn};
        end
        
        z1hat = []; z2hat = [];
        for yn=1:p
            z1hat = [z1hat;z1hatk{yn}];%#ok
            z2hat = [z2hat;z2hatk{yn}];%#ok
        end
        
        xk = xk_new;
    
    end
elapsed_time = toc;
fprintf(1, 'elapsed time: %g [s]\n', elapsed_time);
% -------------------------------------------------

% figure('Name','obs. state vars')
% ax(1) = subplot(2,1,1);
% plot(TimeSamples(1:i),Zhat(1:p ,:), 'LineWidth',1.5);
% grid on; xlabel('time [s]');ylabel('y estimate z_1');
% 
% ax(2) = subplot(2,1,2);
% plot(TimeSamples(1:i),Zhat(p+1:end,:), 'LineWidth',1.5);
% grid on; xlabel('time [s]');ylabel('disturb estimate z_2');
% 
% figure('Name','output vs reference')
% plot(TimeSamples(1:i),Y, 'LineWidth',1.5);
% grid on; xlabel('time [s]'); hold on
% plot(TimeSamples(1:i),Yref.*ones(size(Y)), 'LineWidth',2.5);
% 
% figure('Name','control input')
% plot(TimeSamples(1:i),U, 'LineWidth',1.5);
% grid on; xlabel('time [s]');
% 
errY = Y - Yref.*ones(size(Y));
% errYSTD = std(errY,0,1);
% % ---- per confronto carico anche i dati di macchina -------
% errSTD_ELETTRA = std(errELETTRA,0,2);
% 
% figure;
% subplot(2,2,1);plot(TimeSamples(1:i),errY);grid on;xlabel('samples');ylabel('a.u.');title('Sensor Errors (Orbita)');
% subplot(2,2,2);plot(errYSTD(1:i));grid on;xlabel('samples');ylabel('a.u.');title('Feedback error (RMS)');
% 
% subplot(2,2,3);%plot(TimeSamples(1:i),errELETTRA(1:i));grid on;xlabel('samples');ylabel('a.u.');title('ELETTRA Sensor Errors (Orbita)');
% subplot(2,2,4);%plot(errSTD_ELETTRA(1:i));grid on;xlabel('samples');ylabel('a.u.');title('ELETTRA Feedback error (RMS)');
% 
% figure('Name','control input')
% for input=1:m
%     subplot(4,ceil(m/4),input)
%     plot(TimeSamples(1:i),U(input,:), 'LineWidth',1.5);
%     hold on
%     plot(TimeSamples(1:i),umin(input)*ones(i,1), 'LineWidth',1.5);
%     plot(TimeSamples(1:i),umax(input)*ones(i,1), 'LineWidth',1.5);
% end
% 
% figure('Name','disturbance')
% for input=1:7
%     subplot(4,ceil(m/4),input)
%     plot(TimeSamples(1:i),disturbances(input,:), 'LineWidth',1.5);
% end

    filename = ['ADRC_MPC_7Yx15U_10khz_N',num2str(N),'.mat'];
    save(filename, 'TimeSamples', 'X', 'Y', 'U', 'U0', 'YrefSEQ', 'Zhat', ...
        'T_camp_ESO', 'freq_red_factor', 'elapsed_time', ...
        'U_indices', 'Y_indices')
    
    Nstart = 1000;
    ErrY = errY(:,Nstart:end);
    errYSTD = std(ErrY,0,2);
    errYVAR = var(ErrY,0,2);
    maxErrY = max(ErrY,[],2);
    disturbSTD = std(disturbances,0,2);
    disturbVAR = var(disturbances,0,2);
    maxDisturb = max(disturbances,[],2);
    STD_reduction = abs(errYSTD ./ disturbSTD);
    VAR_reduction = abs(errYVAR ./ disturbVAR);
    maxReduction = abs(maxErrY ./ maxDisturb);
    errYnorm1 = sum(abs(errY(:,Nstart:end)),2);
    
    filename = ['RESULTS_MPC_ADRC_10kHz_N',num2str(N),'.mat'];
    save(filename, 'Nstart', 'ErrY', 'errYSTD', 'errYVAR', 'maxErrY','disturbSTD',...
        'disturbVAR', 'maxDisturb', 'STD_reduction', 'VAR_reduction', ...
        'maxReduction', 'errYnorm1','disturbances')
end
