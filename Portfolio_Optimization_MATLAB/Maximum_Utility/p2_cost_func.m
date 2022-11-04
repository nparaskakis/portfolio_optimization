%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function f = p2_cost_func(w,mu,Sigma,risk_aversion)
    
    
    
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
    
    f = - mu'*w + (1/2)*risk_aversion*w'*Sigma*w;
    
    
    
end