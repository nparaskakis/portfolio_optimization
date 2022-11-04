function x_0_barrier = findFeasibleBarrier(A,b,n,p,m,e,Sigma,var)
        
    % Cost function of primal problem
    f0 = @(x) f0_nbm_phase1(x);
    grad_f0 = @(x) grad_f0_nbm_phase1(x);
    hessian_f0 = @(x) hessian_f0_nbm_phase1(x);
    
    
    % Logarithmic barrier function and its derivatives
    phi = @(x) phi_nbm_phase1(x,Sigma,var);
    grad_phi = @(x) grad_phi_nbm_phase1(x,Sigma,var);
    hessian_phi = @(x) hessian_phi_nbm_phase1(x,Sigma,var);
    
    
    % Cost function of barrier method and its derivatives
    f_t = @(x,t) t*f0(x) + phi(x);
    grad_f_t = @(x,t) t*grad_f0(x) + grad_phi(x);
    hessian_f_t = @(x,t) t*hessian_f0(x) + hessian_phi(x);
    
    
    % Point x_k_bm (solution estimate) at k-th iteration of Barrier method
    x_k_bm = [];

    % To be filled with all x_k_bm through the iterations of the Barrier method
    x_bm = [x_k_bm]; %#ok<NBRAK>

    % To be filled with all f_x_k_bm through the iterations of the Barrier method
    f_x_bm = [];


    % Point x_k_nm (solution estimate) at k-th iteration of the Newton algorithm
    % Starting point, x_0_cvx (feasible point)
    x_k_nm = pinv(A) * b;
    x_k_nm(end+1) = max(-x_k_nm) + 1;
    
    % To be filled with all x_k_nm through the iterations of the Newton algorithm
    x_nm = [x_k_nm]; %#ok<NBRAK>

    % To be filled with all f_x_k_nm through the iterations of the Newton algorithm
    f_x_t_nm = [];


    % Duality gap (to be filled through the iterations of the Barrier method)
    duality_gap_k = [];


    % Parameter a of backtracking line search must be between 0 and 0.5
    alpha = 0.05;

    % Parameter b of backtracking line search must be between 0 and 1
    beta = 0.5;


    % Parameter mu of Barrier method
    mu = 1.5;

    % Step scaler, t_bm, of Barrier method
    t_bm = rand;
    
    
    % Matrix A_phase1
    A_phase1 = [A zeros(p,1)];
    
    
    % Barrier method iterations counter
    iter_bm = 1;
    
    
    % Barrier iterations
    while(1)
        
        % Newton algorithm iterations counter
        iter_nm = 0;
        
        % Newton iterations
        while(1)

            % Value of function f computed at x_k_nm
            f_x_k_t_nm = f_t(x_k_nm,t_bm);
            
            % Store f_x_k_t_nm to f_x_t_nm
            f_x_t_nm = [f_x_t_nm f_x_k_t_nm]; %#ok<AGROW>

            % Gradient of f computed at x_k_nm
            grad_f_x_k_t_nm = grad_f_t(x_k_nm,t_bm);

            % Hessian of f computed at x_k_nm
            hessian_f_x_k_t_nm = hessian_f_t(x_k_nm,t_bm);
            
            % Finding the Newton step, delta_x_k_nm
            block_matrix = [hessian_f_x_k_t_nm A_phase1.' ; A_phase1 zeros(p,p)];
            steps = block_matrix \ [-grad_f_x_k_t_nm ; zeros(p,1)];
            delta_x_k_nm = steps(1:n+1);

            % Newton decrement squared
            lambda_squared = delta_x_k_nm.'*hessian_f_x_k_t_nm*delta_x_k_nm;

            % Newton algorithm termination criterion
            if ((lambda_squared/2) <= e)
                break;
            end
            
            % Step scaler, t_nm, of Newton algorithm
            t_nm = 1;
            
            % Check feasibility
            while (~is_close_to_zero(A*(x_k_nm(1:end-1) + t_nm*delta_x_k_nm(1:end-1)) - b))
                % Updated value of step scaler, t_nm, of Newton algorithm
                t_nm = beta*t_nm;
            end
            
            % Check feasibility
            %x_test = x_k_nm + t_nm*delta_x_k_nm;
            %while ( sum( - x_test(1:end-1) - x_test(end) <= 0) ~= length(x_test) - 1 )
                % Updated value of step scaler, t_nm, of Newton algorithm
                %t_nm = beta*t_nm;
                %x_test = x_k_nm + t_nm*delta_x_k_nm;
            %end
            
            % Backtracking line search
            while (f_t(x_k_nm+t_nm*delta_x_k_nm,t_bm) > f_x_k_t_nm + alpha*t_nm*grad_f_x_k_t_nm.'*delta_x_k_nm)
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
        x_bm = [x_bm x_k_bm]; %#ok<AGROW>
        
        % Value of function f computed at x_k_bm
        f_x_k_bm = f_t(x_k_bm,t_bm);
        
        % Store f_x_k_bm to f_x_bm
        f_x_bm = [f_x_bm f_x_k_bm]; %#ok<AGROW>
    
        % Barrier method termination criterion
        if ((m/t_bm) < e)
            break
        end
        
        % Updated value of step scaler, t_bm, of Barrier method
        t_bm = mu * t_bm;
        
        % Increase Barrier method iterations counter
        iter_bm = iter_bm + 1;

    end
    
    % Feasible initial point found via Barrier method
    x_0_barrier = x_k_bm(1:end-1);
    
end