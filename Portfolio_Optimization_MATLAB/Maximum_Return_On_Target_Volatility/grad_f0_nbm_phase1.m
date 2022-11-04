function grad = grad_f0_nbm_phase1(w_aug)
    
    n = length(w_aug)-1;
    
    grad = zeros(n+1,1);
    
    grad(end) = 1;
    
end