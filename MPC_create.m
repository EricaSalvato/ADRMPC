function [H, F, M, hatS, hatT, Xm_MPC, XM_MPC] = MPC_create(N, Reach_matrix, A_EAS, B_EAS, R_EAS, Q_EAS, P_EAS, ymax)
% http://cse.lab.imtlucca.it/~bemporad/teaching/mpc/imt/1-linear_mpc.pdf
% The procedure is the same of Magni Scattolini: Advanced and multivariable
% control with the introduction of the cost that lead to follow a ref
[r_B_EAS, c_B_EAS] = size(B_EAS);
hatS = Reach_matrix;
for i=1:N-1
    new_col = [zeros(i*r_B_EAS,c_B_EAS);...
        Reach_matrix(1:(end-i*r_B_EAS),:)];
    hatS = [hatS, new_col];%#ok
end

hatT= A_EAS;
for i=2:N
    hatT= [hatT;A_EAS^i];%#ok
end


hatR = [R_EAS;zeros(N*length(R_EAS)-length(R_EAS),length(R_EAS))];
[rhatR,~] = size(hatR);
for i=1:N-1
    new_col = [zeros(i*length(R_EAS),length(R_EAS));R_EAS;zeros(rhatR-i*length(R_EAS)-length(R_EAS),length(R_EAS))];
    hatR = [hatR, new_col];%#ok
end


hatQ = [Q_EAS;zeros((N-1)*length(Q_EAS),length(Q_EAS))];
[rhatQ,~] = size(hatQ);
for i=1:N-1
    if i ~= N-1
        new_col = [zeros(i*length(Q_EAS),length(Q_EAS)); Q_EAS;...
            zeros(rhatQ-i*length(Q_EAS)-length(Q_EAS),length(Q_EAS))];
    else
        new_col = [zeros(i*length(Q_EAS),length(Q_EAS)); P_EAS];
    end
    hatQ = [hatQ, new_col];%#ok
end       

H = 2*(hatR + hatS'*hatQ*hatS);
F = 2*hatS'*hatQ'*hatT;
M = 2*(hatS'*hatQ');

%--- U constraints
ampX_MPC = ymax;
xmin = -ampX_MPC*ones(length(A_EAS),1);
xmax = ampX_MPC*ones(length(A_EAS),1);
Xm_MPC = repmat(xmin,N,1);
XM_MPC = repmat(xmax,N,1);


