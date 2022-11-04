function grad = grad_phi_nbm_phase1(w_aug,Sigma,v)
    
    w = w_aug(1:end-1);
    s = w_aug(end);
    
    n = length(w);
    
    grad = zeros(n+1,1);
    
    for i = (1:1:n)
        
        e_i = zeros(n+1,1);
        e_i(i) = 1;
        
        e_n_plus_1 = zeros(n+1,1);
        e_n_plus_1(end) = 1;
        
        grad = grad + (1/(w(i)+s))*(-e_i-e_n_plus_1) + (1/(s+1-w(i)))*(e_i-e_n_plus_1);
        
    end
    
    grad = grad + (1/(s+v-w'*Sigma*w))*([2*Sigma*w ; -1]);
    
end