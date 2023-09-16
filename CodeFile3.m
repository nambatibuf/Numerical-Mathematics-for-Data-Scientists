function [mu, var] = nambati_final_p3(eigFunc, m, tol, N)
%NAMBATI_FINAL_P3 function which returns variance and mean of number of iteration for random matrix.
%
%   Inputs :
%       eigFunc - eigFunc is a function handle.
%       m  - Size of the square matrix.
%       tol  - Tolerance.
%       N  - Number of samples.
%
%   Outputs :
%       mu - mean.
%       var - variance.

    rng(0);
    iter=[];
    counter=1;
    i=1;
    while i<=N
        rand_matrix=rand(m);
        if ~issymmetric(rand_matrix)
            rand_matrix=(rand_matrix+rand_matrix')/2;
            [D,Q,cur_iter]=eigFunc(rand_matrix,tol);
            iter(counter)=cur_iter;
            counter=counter+1;
            i=i+1;
        end
    end
    mean_v=mean(iter);
    numerator=sum((iter-mean_v).^2);
    var=numerator/(N-1);
    mu=mean_v;
end % nambati_final_p3