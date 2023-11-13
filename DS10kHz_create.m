function [A_DS, B_DS, C_DS, tau, Q_DS, R_DS] = DS10kHz_create(f_camp_ESO, A_MPC, B_MPC, C_MPC, Q, R)
% DS10kHz_create() implements a discrete-time representation of the
% continuous-time system (A_MPC, B_MPC, C_MPC, 0), with sampling 
% frequency f_camp_ESO and using the step-invariant transform 
% for discretizing. 
% Moreover, given the continuous-time cost matrices Q and R, the function 
% DS10kHz_create() computes the corresponding cost matrices of the
% discretized cost function
%
% IN:
% f_camp_ESO          <--> sampling frequency [Hz]
% A_MPC, B_MPC, C_MPC <--> matrices of the continuous-time 
%                          dynamic system to control using         
%                          MPC approach 
% Q, R                <--> matrices of the continuous-time cost
%                          to optimize: 
%                           integral[x'(t)Qx(t) +u'(t)Ru(t)]
%
% OUT:
% A_DS, B_DS, C_DS    <--> matrices of the discrete-time 
%                          dynamic system to control using         
%                          MPC approach
% tau                 <--> the MPC sampling time
%
% Q_DS, R_DS          <--> matrices of the discrete-time cost
%                          to optimize

    % MPC sampling period
    T_camp_MPC = 1/f_camp_ESO;
    tau = T_camp_MPC;
    
    % discretizing the system
    D_MPC = zeros(size(C_MPC,1),size(B_MPC,2));
    sysC_MPC = ss(A_MPC, B_MPC, C_MPC, D_MPC);
    sysD_MPC = c2d(sysC_MPC, T_camp_MPC);
    A_DS = sysD_MPC.A;
    B_DS = sysD_MPC.B;
    C_DS = sysD_MPC.C;
    
    % Cost matrices D.T.
    Q_DS = tau*Q;
    R_DS = tau*R;

end