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
    
    n = s.n;
    
    Sigma = s.Sigma;
    var = s.var;
    
    v1 = s.v1;
    v2 = s.v2;
    
    cost_func = s.cost_func;
    f = @(w) cost_func(w);
    
    % Set warnings off
    warning('off');
    
    % Use CVX to minimize f(x)
    cvx_begin quiet;
        variable w(n);
        minimize(f(w));
        subject to
            ones(n,1)'*w == 1;
            w'*Sigma*w <= var; %#ok<VUNUS>
            v1 <= w <= v2; %#ok<VUNUS,CHAIN>
    cvx_end;
    
    % Result of problem solution via CVX
    w_opt = w;
    f_opt = cvx_optval;
    
    
    
end