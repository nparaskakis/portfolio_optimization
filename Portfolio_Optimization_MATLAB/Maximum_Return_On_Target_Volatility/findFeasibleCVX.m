function w_0 = findFeasibleCVX(n,v1,v2,Sigma,var)
    
    % Use CVX to compute a feasible initial point x_0
    cvx_begin quiet
    variable w(n,1)
    minimize 1
    subject to
        ones(n,1)'*w == 1;
        w'*Sigma*w <= var; %#ok<VUNUS>
        v1 <= w <= v2; %#ok<VUNUS,CHAIN>
    cvx_end
    
    w_0 = w;
    
end