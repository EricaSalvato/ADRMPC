function [Reach_matrix, Obs_matrix] = reach_obsMatrices(N, A_DT, B_DT, C_DT)
% Verify the reachability and the observability properties
    Reach_matrix = B_DT;
    for i=1:N-1
        Reach_matrix = [Reach_matrix;A_DT^i*B_DT]; %#ok
    end
    
    Obs_matrix = C_DT;
    for i=1:N-1
        Obs_matrix = [Obs_matrix;C_DT*A_DT^i]; %#ok
    end
end 