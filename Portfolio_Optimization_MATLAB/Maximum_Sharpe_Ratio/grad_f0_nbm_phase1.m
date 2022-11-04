function grad = grad_f0_nbm_phase1(y_aug)
    
    n = length(y_aug)-1;
    
    grad = zeros(n+1,1);
    
    grad(end) = 1;
    
end