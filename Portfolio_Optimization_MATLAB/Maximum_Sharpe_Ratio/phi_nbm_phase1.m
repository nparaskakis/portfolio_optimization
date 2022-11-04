function phi = phi_nbm_phase1(y_aug,C_hat)
    
    y = y_aug(1:end-1);
    s = y_aug(end);
    
    phi = - sum(log(y + s*ones(size(y)))) - sum(log(C_hat*y + s*ones(size(y))));
    
end