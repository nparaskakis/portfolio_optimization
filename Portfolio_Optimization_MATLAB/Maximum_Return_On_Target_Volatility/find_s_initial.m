function s = find_s_initial(w,Sigma,var)

    s = max([-w ; w-1 ; w'*Sigma*w - var]) + 1;

end