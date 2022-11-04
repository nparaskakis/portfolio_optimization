function y_0_cvx = findFeasibleCVX(n,A,b,C_hat) %#ok<INUSL,STOUT>

    % Use CVX to compute a feasible initial point x_0
    cvx_begin quiet
    variable s
    variable y_0_cvx(n,1)
    minimize 1
    subject to
%         eye(n)*y_0_cvx == pinv(A)*b;
%         - y_0_cvx - s <= 0; %#ok<VUNUS>
%         - C_hat*y_0_cvx - s <= 0; %#ok<VUNUS>
        A*y_0_cvx == b;
        y_0_cvx >= 0; %#ok<VUNUS>
        C_hat*y_0_cvx >= 0; %#ok<VUNUS>
    cvx_end
    
    
    
end