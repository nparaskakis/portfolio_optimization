function [x_opt,f_opt,x_all,f_all,iter] = interior_point_method(s)
    
    
    % Cost function of primal problem

    f0 = s.f0;
    f0 = @(w) f0(w);

    grad_f0 = s.grad_f0;
    grad_f0 = @(w) grad_f0(w);

    hessian_f0 = s.hessian_f0;
    hessian_f0 = @(w) hessian_f0(w);


    % Barrier function and its derivatives

    phi = s.phi;
    phi = @(w) phi(w);

    grad_phi = s.grad_phi;
    grad_phi = @(w) grad_phi(w);

    hessian_phi = s.hessian_phi;
    hessian_phi = @(w) hessian_phi(w);


    % Cost function of barrier method and its derivatives
    
    f_t = @(w,t) t*f0(w) + phi(w);
    
    grad_f_t = @(w,t) t*grad_f0(w) + grad_phi(w);
    
    hessian_f_t = @(w,t) t*hessian_f0(w) + hessian_phi(w);
    
    x_initial = s.x_initial;
    
    % Point x_k_bm (solution estimate) at k-th iteration of Barrier method
    x_k_bm = [];
    
    % To be filled with all x_k_bm through the iterations of the Barrier method
    x_bm_all = [x_k_bm]; %#ok<NBRAK>
    
    % To be filled with all f_x_k_bm through the iterations of the Barrier method
    f_x_bm_all = [];
    
    % Point x_k_nm (solution estimate) at k-th iteration of the Newton algorithm
    % Starting point, x_0_cvx (feasible point)
    x_k_nm = x_initial;

    % To be filled with all x_k_nm through the iterations of the Newton algorithm
    x_nm = [x_k_nm]; %#ok<NBRAK>

    % To be filled with all f_x_k_nm through the iterations of the Newton algorithm
    f_x_nm_all = [];

    % Duality gap (to be filled through the iterations of the Barrier method)
    duality_gap_k = [];

    % Parameter a of backtracking line search must be between 0 and 0.5
    alpha = s.alpha;

    % Parameter b of backtracking line search must be between 0 and 1
    beta = s.beta;

    % Parameter mu of Barrier method
    mu = s.mu;

    % Step scaler, t_bm, of Barrier method
    t_bm = s.t_bm;
    
    e = s.e;
    
    A = s.A;
    
    b = s.b;
    
    p = s.p;
    
    m = s.m;
    
    n = s.n;
    
    C_hat = s.C_hat;
    % Barrier method iterations counter
    iter_bm = 0;

    
    while(1)

        iter_nm = 0;

        % Newton iterations
        while(1)

            % Value of function f computed at x_k_nm
            f_w_k_t_nm = f_t(x_k_nm,t_bm);

            % Store f_x_k_t_nm to f_x_t_nm
            f_x_nm_all = [f_x_nm_all f_w_k_t_nm]; %#ok<AGROW>

            % Gradient of f computed at x_k_nm
            grad_f_x_k_t_nm = grad_f_t(x_k_nm,t_bm);

            % Hessian of f computed at x_k_nm
            hessian_f_x_k_t_nm = hessian_f_t(x_k_nm,t_bm);

            % Finding the Newton step, delta_x_k_nm
            block_matrix = [hessian_f_x_k_t_nm A.' ; A zeros(p,p)];
            steps = block_matrix \ [-grad_f_x_k_t_nm ; zeros(p,1)];
            delta_x_k_nm = steps(1:n);
            
            % Newton decrement squared
            lambda_squared = delta_x_k_nm.'*hessian_f_x_k_t_nm*delta_x_k_nm;

            % Newton algorithm termination criterion
            if ((lambda_squared/2) <= e)
                break;
            end

            % Step scaler, t_nm, of Newton algorithm
            t_nm = 1;

            % Check feasibility
            while ( (~is_close_to_zero(A*(x_k_nm + t_nm*delta_x_k_nm) - b)) ||...
                    (sum(C_hat*(x_k_nm + t_nm*delta_x_k_nm) >= 0) ~= length(x_k_nm)) ||...
                    (sum((x_k_nm + t_nm*delta_x_k_nm) >= 0) ~= length(x_k_nm)) )
                % Updated value of step scaler, t_nm, of Newton algorithm
                t_nm = beta*t_nm;
            end

            % Backtracking line search
            while (f_t(x_k_nm+t_nm*delta_x_k_nm,t_bm) > f_w_k_t_nm + alpha*t_nm*grad_f_x_k_t_nm.'*delta_x_k_nm)
                % Updated value of step scaler, t_nm, of Newton algorithm
                t_nm = beta*t_nm;
            end

            % Next point x_(k+1)_nm of Newton algorithm
            x_k_nm = x_k_nm + t_nm*delta_x_k_nm;

            % Store x_k_nm to x_nm
            x_nm = [x_nm x_k_nm]; %#ok<AGROW>

            % Increase Newton algorithm iterations counter
            iter_nm = iter_nm + 1;

        end

        % Duality gap
        duality_gap_k = [duality_gap_k m/t_bm]; %#ok<AGROW>

        % Next point x_(k+1)_bm of Barrier method
        x_k_bm = x_k_nm;

        % Store x_k_bm to x_bm
        x_bm_all = [x_bm_all x_k_bm]; %#ok<AGROW>

        % Value of function f computed at x_k_bm
        f_x_k_bm = f_t(x_k_bm,t_bm);

        % Store f_x_k_bm to f_x_bm
        f_x_bm_all = [f_x_bm_all f_x_k_bm]; %#ok<AGROW>

        % Barrier method termination criterion
        if ((m/t_bm) < e)
            break;
        end

        % Updated value of step scaler, t_bm, of Barrier method
        t_bm = mu * t_bm;

        % Increase Barrier method iterations counter
        iter_bm = iter_bm + 1;

    end
    
    x_opt = x_k_bm;
    f_opt = f0(x_k_bm);
    x_all = x_bm_all;
    f_all = f_x_bm_all;
    iter = iter_bm;

end