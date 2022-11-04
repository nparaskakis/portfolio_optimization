%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function step = polyak_step(f_x,df_x,f_opt,~,~)
    
    
    
    %% Description
    %
    %  This is a function that computes the polyak step, which is
    %  needed in projected subgradient algorithm.
    
    
    
    %% Input
    %
    %  f_x   -> The value of the cost function at point x.
    %  df_x  -> The value of the subgradient of
    %           the cost function at point x.
    %  f_opt -> The value of the cost function at the optimal point.
    
    
    
    %% Output
    %
    %  step -> The polyak step at a specific point.
    
    
    
    %% Function's body
    
    if (is_close_to_zero(df_x) == 1)
        step = 1;
    else
        step = (f_x - f_opt)/(norm(df_x,2)^2);
    end
    
    
    
end