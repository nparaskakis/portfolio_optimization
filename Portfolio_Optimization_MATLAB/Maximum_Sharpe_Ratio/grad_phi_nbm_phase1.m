function grad = grad_phi_nbm_phase1(y_aug,C_hat)
    
    y = y_aug(1:end-1);
    s = y_aug(end);
    
    n = length(y);
    
    grad = zeros(n+1,1);
    
    for i = (1:1:n)
        
        e_i = zeros(n+1,1);
        e_i(i) = 1;
        
        e_n_plus_1 = zeros(n+1,1);
        e_n_plus_1(end) = 1;
        
        c_hat_i = C_hat(i,:)';
        
        grad = grad + (1/(y(i)+s))*(-e_i-e_n_plus_1) + (1/(c_hat_i'*y+s))*(-[c_hat_i ; 0] - e_n_plus_1);
        
    end
        
end