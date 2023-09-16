function [x, iter] = sbalraj2_nambati_pp9(g, x0, eps, delta, itermax)
%SBALRAJ2_NAMBATI_PP9 implements Newton-Raphson's method to evaluate the root of function g(x)
%
%   Inputs :
%       g - function handle g has the following [f, fx] = g(x) where f is the function and fx is the derivative            
%       x0 - Initial guess 
%       eps - Convergence criteria |xn+1-xn|<eps 
%       delta - Divergence criteria |xn+1-xn|>delta.
%       itermax - Maximum allowed iterations 
% 
%   Outputs :
%       x - Calculated root        
%       iter - The number of iterations required to obtain the root

    iter = 0;
    while true
        if g(x0) == 0
            break
        end
        [func , func_d] = g(x0);
        x_next = x0-(func/func_d);

        if func_d == 0
            error("Newton Raphson Failed to converge to root.") 
        end

        if abs(x_next - x0) < eps   
            x0 = x_next;
            break
        end
        if abs(x_next-x0)>delta   
            error("Method diverged.")
        end

        iter = iter + 1;
        x0=x_next;

        if iter > itermax  
            error("Maximum iteration limit exceeded.")
        end

    end   
    x=x0;
    
end %sbalraj2_nambati_pp9

