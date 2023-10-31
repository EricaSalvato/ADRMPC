function [A_EAS, B_EAS, C_EAS, tau, Q_EAS, R_EAS] = EAS10kHz_create(f_camp_ESO, A_MPC, B_MPC, C_MPC, Q, R)

% undersampling
T_camp_MPC = 1/f_camp_ESO;
tau = T_camp_MPC;

%{
% Euler auxiliary system representation
A_EAS = eye(size(A_MPC))+tau*A_MPC;
B_EAS = tau*B_MPC;
C_EAS = tau*C_MPC;
%}

D_MPC = zeros(size(C_MPC,1),size(B_MPC,2));
sysC_MPC = ss(A_MPC, B_MPC, C_MPC, D_MPC);
sysD_MPC = c2d(sysC_MPC, T_camp_MPC);
A_EAS = sysD_MPC.A;
B_EAS = sysD_MPC.B;
C_EAS = sysD_MPC.C;
% Cost matrices D.T.
Q_EAS = tau*Q;
R_EAS = tau*R;
end