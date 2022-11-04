function phi = phi_nbm_phase1(w_aug,Sigma,v)
    
    w = w_aug(1:end-1);
    s = w_aug(end);
    
    phi = - sum(log(w + s*ones(size(w)))) - sum(log(- w + (s+1)*ones(size(w)))) - log(s+v-w'*Sigma*w);
    
end