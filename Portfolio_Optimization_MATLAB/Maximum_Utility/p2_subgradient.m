%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function df = p2_subgradient(w,mu,Sigma,risk_aversion)
    
    
    
    %% Description
    %
    %  This is a function that computes the subgradient associated
    %  with projected subgradient algorithm for exercise 4b.
    
    
    
    %% Input
    %
    %  x -> The point x at which the value of the subgradient will be computed
    %  A -> Matrix of size mxn
    %  b -> Vector of size mx1
    
    
    
    %% Output
    %
    %  df -> The value of the subgradient computed at point x.
    
    
    
    %% Function's body
    
    df = - mu + risk_aversion*Sigma*w;
    
    
    
end