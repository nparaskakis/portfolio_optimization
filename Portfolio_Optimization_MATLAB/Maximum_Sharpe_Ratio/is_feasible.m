function boolean = is_feasible(y,A,b,C_hat)
    
    tmp1 = A*y-b;
    
    for i = (1:1:length(tmp1))
        if (is_close_to_zero(tmp1(i)))
            tmp1(i) = 0;
        end
    end
    
    for i = (1:1:length(y))
        if (is_close_to_zero(y(i)) == 1)
            y(i) = 0;
        end
    end
    
    tmp2 = C_hat*y;
    
    for i = (1:1:length(y))
        if (is_close_to_zero(tmp2(i)) == 1)
            tmp2(i) = 0;
        end
    end
    
    if (is_close_to_zero(tmp1) &&...
        (sum(y >= 0) == length(y)) &&...
        (sum(tmp2 >= 0) == length(y)))
        boolean = 1;
    else
        boolean = 0;
    end

end