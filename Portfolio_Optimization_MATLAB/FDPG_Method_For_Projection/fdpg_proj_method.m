%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function [x_opt,f_opt,x_all,f_all,iter] = fdpg_proj_method(s)
    
    
    
    %% Description
    %
    %  This is a function that implements the Fast Dual Proximal
    %  Gradient (FDPG) method to solve an optimization problem.
    
    
    
    %% Input
    %
    %  Structure s with the following member variables:
    %  L           -> Lipschitz constant.
    %  c           -> Tolerance on finding the optimal point.
    %  f_opt_cvx   -> The value of the cost function at the optimal point.
    %  cost_func   -> Cost function.
    %  constraint  -> Function that checks if the
    %                 constraints are satisfied.
    %  argmax      -> Functiomn that computes the argmax.
    %  proj_1      -> Function that computes the projection_1.
    %  proj_2      -> Function that computes the projection_2.
    %  transform_A -> Function that computes the A().
    %  y_initial   -> Initial point for the algorithm.
    
    
    
    %% Output
    %
    %  x_opt -> Optimal point.
    %  f_opt -> Value of cost function at the optimal point.
    %  x_all -> A matrix with all vectors x_k.
    %  f_all -> A vector with the values of the
    %           cost function f at all points x_k.
    %  iter  -> The number of iterations that the algorithm needed.
    
    
    
    %% Function's body
    
    L = s.L;
    
    constraint = s.constraint;
    constraint = @(x) constraint(x);
    
    cost_func = s.cost_func;
    f = @(x) cost_func(x);
    
    argmax = s.argmax;
    argmax = @(y) argmax(y);
    
    proj_1 = s.proj_1;
    proj_1 = @(x) proj_1(x);
    
    proj_2 = s.proj_2;
    proj_2 = @(x) proj_2(x);
    
    transform_A = s.transform_A;
    transform_A = @(x) transform_A(x);
    
    y_initial = s.y_initial;
    y_k = y_initial;
    
    w_initial = y_initial;
    w_k = w_initial;
    
    t_initial = 1;
    t_k = t_initial;
    
    x_all = [];
    
    f_all = [];
    
    k = 0;
    
    while(1)
        
        x_k = argmax(y_k);
        
        f_x_k = f(x_k);
        
        f_all = [f_all f_x_k]; %#ok<AGROW>
        
        x_all = [x_all x_k]; %#ok<AGROW>
        
%         if ((f_x_k <= c*f_opt_cvx) && (constraint(x_k) == 1))
%             break;
%         end
        
        u_k = argmax(w_k);
        
        y_k_plus_1(:,1) = w_k(:,1) - (1/L)*transform_A(u_k) + (1/L)*proj_1(transform_A(u_k) - L*w_k(:,1));
        
        y_k_plus_1(:,2) = w_k(:,2) - (1/L)*transform_A(u_k) + (1/L)*proj_2(transform_A(u_k) - L*w_k(:,2));

        x_k_plus_1 = argmax(y_k_plus_1);
        
        if (is_close_to_zero(norm((x_k - x_k_plus_1)/x_k,2)) && (constraint(x_k) == 1))
           break;
        end
        
%         if (is_close_to_zero(norm((x_k - x_k_plus_1)/x_k,2)))
%            break;
%         end
        
%         if ((constraint(x_k) == 1))
%            break;
%         end
        
        t_k_plus_1 = (1 + sqrt(1 + 4*(t_k^2)))/2;
        
        w_k_plus_1 = y_k_plus_1 + ((t_k - 1)/t_k_plus_1)*(y_k_plus_1 - y_k);
        
        y_k = y_k_plus_1;
        
        w_k = w_k_plus_1;
        
        t_k = t_k_plus_1;
        
        k = k+1;
        
    end

    x_opt = x_k;
    f_opt = f_x_k;
    iter = k;
    
    
    
end