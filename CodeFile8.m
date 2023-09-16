function[L, U, P] = nambati_sbalraj2_pp8(A)
%NAMBATI_SBALRAJ2_PP8 function that accepts a matrix A of size n cross times n  and computes the LU decomposition of A
%
%   Inputs :
%       A - A matrix A of size n cross times n.
%
%   Outputs :
%       L - Lower triangular matrix.
%       U - Upper traiangular matrix.
%       P - Lower triangular matrix.

    [m,n]=size(A);
    P=eye(n);
    L=eye(n);
    U=A;
    
    for i=1:n-1
        [~, index]=max(abs(U(i:n,i)));
        index=index+i-1;
    
        tol = 1e-10;
        if abs(U(index,i))>tol
            P([i,index],:)=P([index,i],:);
            L([i,index],1:i-1)=L([index,i],1:i-1);
            U([i,index],i:n)=U([index,i],i:n);
    
            for j = i+1:n
                factor=U(j,i)/U(i,i);
                L(j,i)=factor;
                U(j,i:n)=U(j,i:n)-(factor*U(i,i:n));
            end
        end
    end
    end % nambati_sbalraj2_pp8