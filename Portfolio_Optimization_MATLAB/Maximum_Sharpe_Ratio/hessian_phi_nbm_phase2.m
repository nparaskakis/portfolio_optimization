function hessian_phi = hessian_phi_nbm_phase2(y,C_hat)

    n = size(y,1);
    
    hessian_phi = zeros(n,1);

    for i = (1:1:n)
        
        e_i = zeros(n,1) ;
        e_i(i) = 1;
        
        c_i = C_hat(i,:)';
        
        hessian_phi = hessian_phi + (e_i * e_i.') / (y(i))^2 + (c_i * c_i.') / (c_i'*y)^2;
    
    end
    
end