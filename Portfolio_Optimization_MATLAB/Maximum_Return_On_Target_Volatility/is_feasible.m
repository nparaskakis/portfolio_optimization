function boolean = is_feasible(w,A,b,Sigma,var)
    
    tmp1 = A*w-b;
    
    for i = (1:1:length(tmp1))
        if (is_close_to_zero(tmp1(i)))
            tmp1(i) = 0;
        end
    end
    
    for i = (1:1:length(w))
        if (is_close_to_zero(w(i)) == 1)
            w(i) = 0;
        end
    end
    
    for i = (1:1:length(w))
        if (is_close_to_zero(1 - w(i)) == 1)
            w(i) = 1;
        end
    end
    
    tmp2 = w'*Sigma*w - var;
    
    for i = (1:1:length(w))
        if (is_close_to_zero(tmp2) == 1)
            tmp2(i) = 0;
        end
    end
    
    if (is_close_to_zero(tmp1) &&...
        (sum(w >= 0) == length(w)) &&...
        (sum(w <= 1) == length(w)) &&...
        (sum(tmp2 < 0) == 1))
        boolean = 1;
    else
        boolean = 0;
    end

end