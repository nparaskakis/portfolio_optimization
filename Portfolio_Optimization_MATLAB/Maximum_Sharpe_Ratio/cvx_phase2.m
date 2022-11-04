%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function [w_opt,f_opt] = cvx_phase2(s)
    
    
    
    %% Description
    %
    %  This is a function that solves the problem
    %  that corresponds to exercise 12b using CVX.
    
    
    
    %% Input
    %
    %  Structure s with member variables:
    %  n         -> Dimension n.
    %  A         -> Matrix of size mxn.
    %  b         -> Vector of size mx1.
    %  v1        -> Left boundary of interval of x.
    %  v2        -> Right boundary of interval of x.
    %  cost_func -> Cost function of problem.
    
    
    
    %% Output
    %
    %  x_opt -> The optimal point.
    %  f_opt -> The value of the cost function at the optimal point.
    
    
    
    %% Function's body
    
    n = s.n; %#ok<NASGU>
    mu = s.mu;
    C_hat = s.C_hat;
    
    cost_func = s.cost_func;
    f = @(y) cost_func(y);
    
    % Set warnings off
    warning('off');
    
    % Use CVX to minimize f(x)
    cvx_begin quiet;
        variable y(n);
        minimize(f(y));
        subject to
            mu'*y == 1;
            y >= 0; %#ok<VUNUS>
            C_hat*y >= 0; %#ok<VUNUS>
    cvx_end;
    
    % Result of problem solution via CVX
    y_opt = y;
    f_opt = cvx_optval;
    
    w_opt = y_opt/sum(y_opt);
    
    
end