function [D, Q, iter] = nambati_final_p2s(A, tol)
%NAMBATI_FINAL_P2 function is based on pure QR algorithm which returns D, Q and Iter using Wilkson shift.
%
%   Inputs :
%       A - Real and symmetric matrix.
%       tol  - Tolerance for convergence.
%
%   Outputs :
%       D - Matrix with Diagnal eigen values as elements.
%       Q - Eigenvector.
%       iter - Number of iteration it took to converge.

    T=A;
    iter=0;
    U=eye(size(T));
    [m,n]=size(A);
    while true
        if n >= 2
            A = A(n-1:n, n-1:n);  % Select the bottom-right 2x2 submatrix
        else
            error('Original matrix is too small to create a 2x2 matrix.');
        end
        trace_A = trace(A);
        det_A = det(A);
        temp=sqrt(trace_A^2 - 4 * det_A);
        u1 = (trace_A+temp) / 2;
        u2 = (trace_A-temp) / 2;
        if abs(u1-A(2,2))<abs(u2-A(2,2))
            shift=u1;
        else
            shift=u2;
        end
        [Q,R]=qr(T-shift*eye(size(T)));
        T=R*Q+shift*eye(size(T));
        U=U*Q;
        iter=iter+1;
        if (norm(diag(T,-1)/norm(diag(T)))<tol)
            D=T;
            break
        end
    end % nambati_final_p2