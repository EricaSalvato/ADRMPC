function [Lp, Lc, Aeso, Bes, BesD,...
    index_ESO_stable] = ESO_create(threshold,f_camp_ESO, p)
% -------------------
% cfr. [1] Han, "From PID to Active Disturbance Rejection Control," 
% in IEEE Transactions on Industrial Electronics, vol. 56, no. 3, 
% pp. 900-906, March 2009, doi: 10.1109/ TIE.2008.2011621.
% -------------------

% state equation of the system
%
% \dot{x} = Ax + Bu
%      y  = Cx + F(t)
%
%
% \dot{y} = CBu + W(t) <-- W(t) dynamics to wash out and 
%                          first derivative of the output additive 
%                          disturbance F(t)
% CBu = U = [U1, U2, U3]^T
%
% then with ESO converging and deleting G(t) we obtain
%
% \dot{y1} = U1
% \dot{y2} = U2
% \dot{y3} = U3
% -----
% the ESO structure
%
% \dot{Z1} = Z2 + U
% \dot{Z2} = W(t)
%        Y = Z1
%
% Y =[y1, y2, y3]^T
%

T_camp_ESO = 1/f_camp_ESO;
for yn=1:p
        Aes{yn} = [0, 1;...
            0, 0];%#ok
    
    Bes{yn} = [1; ...
            0];%#ok
    
    Ces{yn} = [1, 0];%#ok

    sysESO = ss(Aes{yn}, Bes{yn}, Ces{yn}, 0);
    % c2d using ZOH --> see [2]
    sysESO_td = c2d(sysESO, T_camp_ESO, 'zoh');
    AesD{yn} = sysESO_td.A;%#ok
    BesD{yn} = sysESO_td.B;%#ok
    CesD{yn} = sysESO_td.C;%#ok

    % --> current discrete estimator scheme <--
    % [2] R. Miklosovic, A. Radke and Zhiqiang Gao, 
    % "Discrete implementation and generalization 
    %  of the extended state observer," 
    % 2006 American Control Conference, 2006
    %
    om_eso_tc = 5000*2*pi; % rad/s <-> 5 kHz full band ESO 
                      % 
    om_beta = exp(-om_eso_tc*T_camp_ESO); % [2], Eq. (14) 
    % as in [2], Eqs (13) & (14)


    Lp{yn} = [2*(1-om_beta); ...
               ((om_beta^2)-1+2*(1-om_beta))/T_camp_ESO];%#ok
    Lc{yn} = AesD{yn} \ Lp{yn}; %#ok % inv(AesD{yn})*Lp; % [2], Eq. (9)
    
    Aeso{yn} = AesD{yn} -Lp{yn} *CesD{yn}; %#ok
    Beso{yn} = [BesD{yn} , Lp{yn} ];%#ok
    Ceso{yn} = [eye(size(AesD{yn}))-Lc{yn}*CesD{yn} ];%#ok
    Deso{yn} = [-Lc{yn}*zeros(1) , Lc{yn}];%#ok
    
    % index from which we suppose ESO is stable (MPC does not work until ESO is not stable)
    index_ESO_stable = ceil(threshold/(-log(exp(-om_eso_tc/f_camp_ESO))));

end

