function hessian_phi = hessian_phi_nbm_phase2(w,Sigma,v)

    n = size(w,1);
    
    hessian_phi = zeros(n,1);

    for i = (1:1:n)
        e_i = zeros(n,1) ;
        e_i(i) = 1;
        hessian_phi = hessian_phi + (e_i * e_i.') / w(i)^2 + (e_i * e_i.') / (w(i)-1)^2;
    end
    
    hessian_phi = hessian_phi + (1/(w'*Sigma*w - v)^2)*4*Sigma*(w*w')*Sigma' + (1/(v-w'*Sigma*w))*2*Sigma';

end