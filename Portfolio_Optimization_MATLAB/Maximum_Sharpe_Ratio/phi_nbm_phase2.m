function phi = phi_nbm_phase2(y,C_hat)
    
    phi = - sum(log(y)) - sum(log(C_hat*y));
    
end