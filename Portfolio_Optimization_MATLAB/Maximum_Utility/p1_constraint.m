%%  Advanced Topics in Convex Optimization
%%  Nikolaos Paraskakis - A.M.: 2018030027



function boolean = p1_constraint(x,A,b,v1,v2)
    
    
    
    %% Description
    %
    %  This is a function that checks if a point x satsfies the
    %  constraints that correspond to the problem of exercise 12b.
    
    
    
    %% Input
    %
    %  x  -> Vector of size nx1.
    %  A  -> Matrix of size mxn.
    %  b  -> Vector of size mx1.
    %  v1 -> Left boundary of interval of x.
    %  v2 -> Right boundary of interval of x.
    
    
    
    %% Output
    %
    %  boolean -> 1: if the constraints are satisfied.
    %             0: if the constraints are not satisfied.
    
    
    
    %% Function's body
    
    for i = (1:1:length(x))
        if ((x(i) <= v1(i)) && (is_close_to_zero(v1(i) - x(i)) == 1))
            x(i) = v1(i);
        end
    end
    
    for i = (1:1:length(x))
        if ((x(i) >= v2(i)) && (is_close_to_zero(v2(i) - x(i)) == 1))
            x(i) = v2(i);
        end
    end
    
    if (is_close_to_zero(A*x-b) &&...
            sum(v1 <= x) == length(x) &&...
            sum(v2 >= x) == length(x))
        boolean = 1;
    else
        boolean = 0;
    end
    
    
    
end