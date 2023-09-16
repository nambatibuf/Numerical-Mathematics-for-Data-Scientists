function[l, v] = nambati_pp11(A, v0, tol, itermax)
%NAMBATI_PP11 function implements the Rayleigh-Quotient iteration for eigenvalue and eigenvector calculation.
%
%   Inputs :
%       A - Real and square matrix.
%       v0  - Initial vector.
%       tol  - Tolerance for convergence.
%       itermax  - Maximum number of iterations.
%
%   Outputs :
%       l - Eigenvalue.
%       v - Eigenvector.

    i=0;
    v=v0/norm(v0);
    l=v.' * A * v;
    [m,~]=size(A);
    P=eye(m);
    while 1
        v0=(A-l*P)\v;
        v=v0/norm(v0);
        ln=v.' * A * v;
        if i>=itermax
            error("Maximum Iteration reach");
        elseif abs((ln-l)/ln)<=tol
            l=ln;
            break;
        else
            l=ln;
        end
        i=i+1;
    end % nambati_pp11