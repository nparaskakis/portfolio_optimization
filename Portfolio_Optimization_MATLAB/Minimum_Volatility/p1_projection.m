function w_proj = p1_projection(d,n,v1,v2)
    
    m = 2;
    L = m+1;
    
    A = ones(1,n);
    b = 1;
    
    y_initial = rand(n,2);
    
    s = struct('y_initial', y_initial, 'L',L, 'cost_func',@(w) p1_cost_func(w,d), 'argmax',@(w) p1_argmax(w,d), 'transform_A',@(u) p1_transform_A(u), 'proj_1',@(w) p1_proj_1(w,A,b), 'proj_2',@(w) p1_proj_2(w,v1,v2), 'constraint',@(w) p1_constraint(w,A,b,v1,v2));
    
    w_proj = fdpg_proj_method(s);
    
end