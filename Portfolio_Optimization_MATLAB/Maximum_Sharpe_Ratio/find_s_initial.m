function s = find_s_initial(y,C_hat)

    s = max([-y ; -C_hat*y]) + 1;

end