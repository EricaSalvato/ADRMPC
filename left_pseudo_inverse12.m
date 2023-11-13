function A_pINV_12 = left_pseudo_inverse12(Matrix)
    % compute the pseudoinverse {1,2} matrix of "Matrix"
    % Ref:
    % A. Ben-Israel and T. N. Greville, 
    % Generalized inverses: theory and applications,
    % Springer Science & Business Media, 2003, 
    % pag. 258 - eq. (5.19)

    A=Matrix;
    [Ua,~,Va] = svd(A);
    % A = Ua*Sa*Va'
    %
    % Ua'*A*Va = Sa
    P = Ua';
    Q = Va;
    pA = P*A*Q;
    r = rank(pA);
    A11=pA(1:r,1:r);
    [m,n]=size(A);
    dA = diag(A11);
    dA1=1./dA; % inv(A11) = diag(dA1)
    A_pINV_12 = Q*[diag(dA1) zeros(r,m-r);...
                    zeros(n-r,m)]*P;

end