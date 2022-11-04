function phi = phi_nbm_phase2(w,Sigma,v)

    phi = - sum(log(w)) - sum(log(1-w)) - log(v - w'*Sigma*w) ;
    
end