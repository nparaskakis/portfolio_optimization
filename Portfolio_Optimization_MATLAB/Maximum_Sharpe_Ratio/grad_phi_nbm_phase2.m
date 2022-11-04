function grad_phi = grad_phi_nbm_phase2(y,C_hat)

    n = size(y,1);
    
    grad_phi = zeros(n,1);
    
    for i = (1:1:n)
        
        e_i = zeros(n,1);
        e_i(i) = 1;
        
        c_i = C_hat(i,:)';
        
        grad_phi = grad_phi - e_i/y(i) - c_i/(c_i'*y);
    
    end
        
end