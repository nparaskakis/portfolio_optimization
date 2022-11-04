function grad_phi = grad_phi_nbm_phase2(w,Sigma,v)

    n = size(w,1);
    
    grad_phi = zeros(n,1);
    
    for i = (1:1:n)
        e_i = zeros(n,1);
        e_i(i) = 1;
        grad_phi = grad_phi - e_i/w(i) - e_i/(w(i)-1);
    end
    
    grad_phi = grad_phi - (1/(w'*Sigma*w-v))*2*Sigma*w;
    
end