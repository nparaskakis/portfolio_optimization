%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function step = dynamic_step(~,df_x,~,k,Lf)
    
    
    
    %% Description
    %
    %  This is a function that computes the dynamic step, which is
    %  needed in projected subgradient algorithm.
    
    
    
    %% Input
    %
    %  df_x  -> The value of the subgradient of the
    %           cost function at a specific point.
    %  k     -> Number of iteration.
    %  Lf    -> Lipschitz constant of f over set C.
    
    
    
    %% Output
    %
    %  step -> The dynamic step at a specific point.
    
    
    
    %% Function's body
    
    if (is_close_to_zero(df_x) == 1)
        step = 1/Lf;
    else
        step = 1/((norm(df_x,2))*(sqrt(k+1)));
    end
    
    
    
end