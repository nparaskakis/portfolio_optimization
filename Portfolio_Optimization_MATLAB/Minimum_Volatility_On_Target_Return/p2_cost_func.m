%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function f = p2_cost_func(w,Sigma)
    
    
    
    %% Description
    %
    %  This is a function that computes the value of the cost function,
    %  that corresponds to the problem of exercise 12b, at a specific
    %  point x.
    
    
    
    %% Input
    %
    %  x -> Vector of size nx1.
    %  d -> Vector of size nx1.
    
    
    
    %% Output
    %
    %  f -> The value of the cost function computed at point x.
    
    
    
    %% Function's body
    
    f = (1/2)*w'*Sigma*w;
    
    
    
end