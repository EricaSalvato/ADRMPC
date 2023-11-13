function [A_MPC, B_MPC, C_MPC] = TC_MPC_create(p)
% create the system to be controlled by the MPC
% p discrete-time integrators, completely decoupled each other

    A_MPC = zeros(p,p);
    B_MPC = eye(p);
    C_MPC = eye(p);

end