function hessian = hessian_phi_nbm_phase1(y_aug,C_hat)
    
    y = y_aug(1:end-1);
    s = y_aug(end);
    
    n = length(y);
    
    hessian = zeros(n+1,n+1);
    
    for i = (1:1:n)
        
        e_i = zeros(n+1,1);
        e_i(i) = 1;
        
        e_n_plus_1 = zeros(n+1,1);
        e_n_plus_1(end) = 1;
        
        c_hat_i = C_hat(i,:)';
        
        hessian = hessian + (1/(-y(i)-s)^2)*(-e_i-e_n_plus_1)*(-e_i-e_n_plus_1)' + (1/(-c_hat_i'*y-s)^2)*(-[c_hat_i ; 0] - e_n_plus_1)*(-[c_hat_i ; 0] - e_n_plus_1)';
    
    end
        
end