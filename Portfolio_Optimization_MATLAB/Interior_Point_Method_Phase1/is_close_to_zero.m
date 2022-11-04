%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function boolean = is_close_to_zero(x)
    
    
    
    %% Description
    %
    %  This is a function that checks if the elements of a vector x
    %  (or matrix x) are all close to zero.
    
    
    
    %% Input
    %
    %  x -> Vector x or matrix x
    
    
    
    %% Output
    %
    %  boolean -> 1 (true)  if x is close to zero.
    %             0 (false) if x is not close to zero.
    
    
    
    %% Function's body
    
    if (sum(abs(x) < 10^(-5)*ones(size(x)),'all') == numel(x))
        boolean = 1;
    else
        boolean = 0;
    end
    
    
    
end