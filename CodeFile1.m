function [D, Q, iter] = nambati_final_p1(A, tol)
%NAMBATI_FINAL_P1 function is based on pure QR algorithm which returns D, Q and Iter.
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
    while 1
        [Q,R]=qr(T);
        T=R*Q;
        U=U*Q;
        iter=iter+1;
        disp((norm(diag(T,-1))))
        disp((norm(diag(T))))
        if (norm(diag(T,-1)/norm(diag(T)))<tol)
            D=T;
            break
        end
    end
end % nambati_final_p1
