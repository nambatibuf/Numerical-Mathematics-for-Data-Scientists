function [length] = sbalraj2_nambati_pp5(A, i, j)
%SBALRAJ2_PP5 The function computes and returns the minimum path length 
%             between point 'i' and 'j' within the graph using Adjacency
%             Matrix. The function returns an error for path length greater
%             than 15
%   Input Parameters
%   A = Adjacency Matrix
%   i = Initial or departure point
%   j = Final or Goal point

    loop_count = 0;  %Initializing Loop Variable
    T = A;           % Storing A in Temp variable
    while true
        loop_count = loop_count + 1;   %Counter + 1
        if loop_count > 15             % Check for length > 15
            msg = 'Path length must be less than 15';
            error(msg);                % Passing above error
        end   
    
        A = T^loop_count;             % As per corollary of path length 
        if A(i,j) ~= 0                
            length = loop_count;
            break
    
        end    
    end
end

