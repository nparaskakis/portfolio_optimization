%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function [x_opt,f_opt,x_best,f_best,iter] = projected_subgradient_method(s)
    
    
    
    %% Description
    %
    %  This is a function that implements the projected
    %  subgradient algorithm to solve an optimization problem.
    
    
    
    %% Input
    %
    %  Structure s with the following member variables:
    %  print_val   -> Print log in console (1 for true, 0 for false).
    %  c           -> Tolerance on finding the optimal value.
    %  Lf          -> Lipschitz constant for the cost function f.
    %  f_opt_cvx   -> The value of the cost function at the optimal point.
    %  x_initial   -> Initial point (vector) x.
    %  cost_func   -> Problem's cost function f.
    %  subgradient -> Subgradient function of the cost function f.
    %  projection  -> Projection function onto constraint set C.
    %  stepsize    -> Stepsize function for the algorithm (Polyak/Dynamic).
    
    
    
    %% Output
    %
    %  x_opt  -> The optimal point.
    %  f_opt  -> The value of the cost function f at the optimal point.
    %  x_best -> A matrix with vectors x_k that correspond to f_k_best.
    %  f_best -> A vector with the values of f_k_best.
    %  iter   -> The number of iterations that the algorithm needed.
    
    
    
    %% Function's body
    
    print_val = s.print_val;
    
    c = s.c;
    
    Lf = s.Lf;
    
    f_opt_cvx = s.f_opt_cvx;
    
    x_initial = s.x_initial;
    x_k = x_initial;
    
    cost_func = s.cost_func;
    f = @(x) cost_func(x);
    
    subgradient = s.subgradient;
    df = @(x) subgradient(x);
    
    projection = s.projection;
    proj = @(x) projection(x);
    
    stepsize = s.stepsize;
    step = @(f_x_k,df_x_k,f_opt_cvx,k,Lf) stepsize(f_x_k,df_x_k,f_opt_cvx,k,Lf);
    
    x_k_best = x_initial;
    f_k_best = f(x_k);
    
    x_best = [];
    f_best = [];
    
    k = 0;
    
    if (print_val == 1)
        disp('\nProcess of projected subgradient method:');
    end
    
    while(1)
        
        f_x_k = f(x_k);
        
        df_x_k = df(x_k);
        
        if (f_x_k < f_k_best)
            f_k_best = f_x_k;
            x_k_best = x_k;
        end
        
        x_best = [x_best x_k_best]; %#ok<AGROW>
        
        f_best = [f_best f_k_best]; %#ok<AGROW>
        
        if (print_val == 1)
            fprintf('\nk = %d\n',k);
            fprintf('x_k = %s\n',print_vector_inline(x_k));
            fprintf('f(x_k) = %.4f\n',f_x_k);
            fprintf('df(x_k) = %s\n',print_vector_inline(df_x_k));
        end
        
        if (f_x_k <= c*f_opt_cvx)
            break;
        end
        
        t_k = step(f_x_k,df_x_k,f_opt_cvx,k,Lf);
        
        x_k_plus_1 = proj(x_k - t_k*df_x_k);
        
%         if (is_close_to_zero(norm((x_k-x_k_plus_1)/x_k,2)))
%             break;
%         end
        
        x_k = x_k_plus_1;
        
        k = k+1;
        
    end
    
    x_opt = x_k;
    f_opt = f_x_k;
    iter = k;
    
    
    
end