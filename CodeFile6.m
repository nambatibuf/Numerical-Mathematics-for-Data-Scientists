function [R] = sbalraj2_nambati_pp6(A)
%SBALRAJ2_PP6 Returns the Reduced Row Echelon Form of matrix 'A' using 
%             Gaussian Elimination with Pivoting in 'R'
%   A - Input Matrix
%   R - RREF of Matrix A
    
    m = size(A, 1); % Total rows in input matrix
    n = size(A, 2); % Total columns in input matrix

    p = 1;          %Initializing Matrix Iterators to 1st element
    q = 1;

    while p<=m && q<=n
        if p == 1            %Finding the row with largest value in 1st -
            f = A(:, q);     %-column
            w = max(f);
            e = 0;
            for s = 1:max(size(f))
                if f(s) == w
                    e = s;
                    break
                end    
                e = e + 1;
            end  
            k = e;
            A([k p], :) = A([p k], :); % Swapping largest column value rows
        end
         if abs(A(p, q)) > 10^(-16)   % Checking whether 0 or not
            h = A(p, q);              % storing 1st non-zero element
            for t = q : n
                A(p,t) = A(p,t)*(1/h);% Creating Pivots
            end   
            for k=[1:p-1 p+1:m]       % Making other values in that column- 
              h = A(k, q);            % -zero
              for y = q:n
                  A(k,y) = A(k,y) - A(p,y)*h ;
              end    
            end
            p = p + 1;               %Next row
            q = q + 1;               %Next Column
         else
             q = q + 1;              %If encountered Zero/Next Column

         end    
    end
    R = A;                           %Passing Transformed Matrix to Output

end

