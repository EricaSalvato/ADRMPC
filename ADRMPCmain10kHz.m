clear all%#ok
close all
clc
N_values = [4,8,12,16,20]; % Control and prediction horizon

% exploring different possible control and prediction horizon
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
    T_camp_sys = T_camp_ESO;   % plant sampling period [s]
    % --- operating frequencies ---

    % ---- simulation settings ----
    simTime = 5; % time horizon: [0 simTime]
    nsample = floor(simTime/T_camp_sys);
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
    % ------ continuous-time model ------
    
    % --- sampled-time model ---
    sysD = c2d(sysC, T_camp_sys);
    Ad = sysD.A;
    Bd = sysD.B;
    Cd = sysD.C;
    Dd = sysD.D;
    % --- sampled-time model ---
    
    % --- C*B pseudo-inverse ---
    % CBinverse useful also if by using it 
    % we respect constraints on the plant input
    % (avoid to use the LP) 
    CBinv_prov = pinv(Cc*Bc);
    % --- C*B pseudo-inverse ---
    
    % --- additive output disturbances ----
    disturbances = DT_disturb_create(p, nsample); 
    % --- additive output disturbances ----
     
    %{
    figure('Name','Output Disturbance');
    plot(TimeSamples, disturbances, 'LineWidth', 1.5);
    grid on; xlabel('Time [s]'); 
    ylabel('Output Disturbance [mm]');
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% ESO 4 ADRC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold = 41;
            
    [Lp, Lc, Aeso, Bes, ...
        BesD, index_ESO_stable] = ESO_create(threshold, f_camp_ESO, p);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% SoMPC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    freq_red_factor = 1; 
    % frequency reduction factor of MPC 
    % with respect to the ESO sampling
    % frequency

    % -------- cost weights ------------------
    R_weights = 10; %Control input cost weights
    Q_weights = 10000;  % State cost weights
    % -------- cost matrices ------------------

    N = N_values(n_index);
    fprintf(1, 'Horizon_value: %d\n', N);
   
    % --- continuous time system viewed by MPC ---
    [A_MPC, B_MPC, C_MPC] = TC_MPC_create(p);
    [~,c_B] = size(B_MPC);
    
    % --- continuous time cost matrices ---
    R = R_weights*eye(c_B);
    Q = Q_weights*eye(size(A_MPC)); 
   
    % --- discretizing the MPC-related system ---    
%     [A_EAS, B_EAS, C_EAS, tau, Q_EAS, R_EAS] = ...
%         EAS10kHz_create(f_camp_ESO, A_MPC, B_MPC, C_MPC, Q, R);
    [A_DT, B_DT, C_DT, tau, Q_DT, R_DT] = ...
        DS10kHz_create(f_camp_ESO, A_MPC, B_MPC, C_MPC, Q, R);
    
    
    [r_B_EAS, c_B_EAS] = size(B_DT);
    
    % --- N step reachability test ---
    [Reach_matrix, Obs_matrix] = reach_obsMatrices(N, A_DT, B_DT, C_DT);
    assert(rank(Reach_matrix)==length(A_DT),'Not reachable')
    assert(rank(Obs_matrix)==length(A_DT),'Not observable')
    
    % --- terminal weigth matrix 4 stability ensurance ---
    [P_DT,~,~] = idare(A_DT,B_DT,1/2*Q_DT,1/2*R_DT,[],[]);  
    % Discrete time riccati equation Matlab R2021
    
    %[P_DT,~,~] = dare(A_DT,B_DT,1/2*Q_DT,1/2*R_DT,[],[]); 
    % Discrete time Riccati equation Matlab R2017
    
    % --- create the MPC-related matrices ---
    [H, F, M, hatS, hatT, ...
        Xm_MPC, XM_MPC] = MPC_create(N, Reach_matrix,...
                                        A_DT, B_DT, R_DT, ...
                                        Q_DT, P_DT, ymax);
    % --- Quadprog settings 4 MPC ---
    H_quadprog = H;
    A_quadprog_no_uoConstr = [hatS;-hatS];
    Aeq_quadprog = hatS(end-length(B_DT)+1:end,:);
    
    % --- storages for data savings
    U0 = []; % MPC storage
    Jcomp = 0; % MPC Cost value
    
    
    % ---- reference settings ----
    % golden orbit [mm]
    load("goldenORBITS.mat")
    Yref = H_Orbits_REFS;
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
    Ssvd(Ssvd<0.00001) = 0; % setting by trials 
    Yfactor4linprog = Usvd*Ssvd*Vsvd';
    assert(rank(Yfactor4linprog) == m-p,'Condition for y_linprog not ok')
    [row_Yfactor4linprog,col_Yfactor4linprog] = size(Yfactor4linprog);
    
    %%% ---- ESO bounds init ----
    z2hat_max=1/10*ones(length(z2hat),1);
    z2hat_min=-1/10*ones(length(z2hat),1);
    %%% ---- Loop START ----
    fprintf(1, 'Starting simulation loop:\n\r %d simulation steps \n\r' ,nsample )
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
                beq_quadprog = YrefSEQ(:,i) - A_DT^N*x0_MPC;
                
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
                u0 = u_star(1:length(B_DT),1); 
                Jcomp = Jcomp+1/2*x0_MPC'*Q_DT*x0_MPC+1/2*u0'*R_DT*u0;
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
            uk= pinv(Cc*Bc)*UK; 
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

    filename = ['ADRC_MPC_7Yx15U_10khz_N',num2str(N),'.mat'];
    save(filename, 'TimeSamples', 'X', 'Y', 'U', 'U0', 'YrefSEQ', 'Zhat', ...
        'T_camp_ESO', 'freq_red_factor', 'elapsed_time')
    
    Nstart = 1000; % avoiding to analyse the initial transient
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
