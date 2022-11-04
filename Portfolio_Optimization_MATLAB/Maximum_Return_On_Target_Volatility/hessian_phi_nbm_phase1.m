function hessian = hessian_phi_nbm_phase1(w_aug,Sigma,v)
    
    w = w_aug(1:end-1);
    s = w_aug(end);
    
    n = length(w);
    
    hessian = zeros(n+1,n+1);
    
    for i = (1:1:n)
        
        e_i = zeros(n+1,1);
        e_i(i) = 1;
        
        e_n_plus_1 = zeros(n+1,1);
        e_n_plus_1(end) = 1;
        
        hessian = hessian + (1/(-w(i)-s)^2)*(-e_i-e_n_plus_1)*(-e_i-e_n_plus_1)' + (1/(w(i)-1-s)^2)*(e_i-e_n_plus_1)*(e_i-e_n_plus_1)';
    
    end
    
    hessian = hessian + (1/(w'*Sigma*w-v-s)^2)*([2*Sigma*w ; -1])*([2*Sigma*w ; -1])' + (1/(s+v-w'*Sigma*w))*([2*Sigma' zeros(7,1) ; zeros(1,8)]);
    
end