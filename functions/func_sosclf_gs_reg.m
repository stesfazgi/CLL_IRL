function [Vc_star, B1c_star, B2c_star] = func_sosclf_gs_reg(x, f, g, opt, data_pos, data_vel)
% SOSCLF algorithm with gradientsquare constraint with regularization term

    degV = opt.degV;
    degB = opt.degB;
    gammaV = opt.gammaV;
    gammaB = opt.gammaB;
    n_iter = opt.n_iter;
    init_scale = opt.init_scale;

    sossettings = sdpsettings('solver', 'mosek', 'verbose', 0);

    % Track total algorithm time
    t_start = tic;

    % Create symbolic expressions
    [V, V_coeff, ~] = polynomial(x, degV, 2);
    [VR, VR_coeff, ~] = polynomial(x, degV, 2);
    
    dVdx = jacobian(V, x);
    dVRdx = jacobian(VR, x);
        
    %%%%%%%%%%%% Initial guesses
    V_init = zeros(size(V_coeff));
    V_init(1) = init_scale;
    V_init(3) = init_scale;
    
    assign(V_coeff, V_init);
    
    % Create objective for GS iteration
    [cnstrOBJ, epivar] = func_create_accel_cone(x, f, g, 0.5, dVdx, data_pos, data_vel);
    
    % Create objective for regularization
    epivar_reg = sdpvar(1, 1);
    cnstrREG = cone(V_coeff, epivar_reg);

    % Create proxy constraint
    gs_cnstr = dVdx*dVRdx';
    
    % Create underbounding polynomial
    bnd_poly = 1e-2 * sum(x.^2);
    
    gs_losses = zeros(n_iter, 1);
    all_V_coeffs = zeros(size(V_coeff, 1), n_iter);
    for i=1:n_iter
        
        gs_curr = replace(gs_cnstr, VR_coeff, value(V_coeff));
        
        resV = solvesos([sos(V), sos(gs_curr - bnd_poly), cnstrOBJ, cnstrREG], (epivar + gammaV*epivar_reg), sossettings, V_coeff);
        
        gs_losses(i) = value(epivar);
        all_V_coeffs(:, i) = value(V_coeff);
    end
    
    % Pick solution with lowest loss
    [~, best_idx] = min(gs_losses);
    Vc_star = all_V_coeffs(:, best_idx);

    % Create beta vector
    n_states = size(data_pos, 1);
    B = sdpvar(n_states, 1);
    B_coeffs = sdpvar(0, 0);
    for i=1:n_states
        [Bx, Bx_coeffs, ~] = polynomial(x, degB);
        B(i) = Bx;
        B_coeffs = [B_coeffs; Bx_coeffs];
    end
    B_init_coeffs = zeros(size(Bx_coeffs, 1), 1);
    B_init_coeffs(1) = 0.5;
    B_init_coeffs = repmat(B_init_coeffs, n_states, 1);

    
    % Create objective for B fitting
    dVdx_curr = replace(dVdx, V_coeff, Vc_star);
    [cnstrOBJ, epivar] = func_create_accel_cone(x, f, g, diag(B), dVdx_curr, data_pos, data_vel);

    % Create objective for regularization
    cnstrREGB = cone(B_coeffs, epivar_reg);
    
    % Create second lyapunov constraint
    f_cl_curr = f - g*(diag(B)*g'*dVdx_curr');
    l2_cnstr = -dVdx_curr*f_cl_curr;
    
    resB = solvesos([sos(B(1) - bnd_poly), sos(B(2) - bnd_poly), sos(l2_cnstr - bnd_poly), cnstrOBJ, cnstrREGB], (epivar + gammaB*epivar_reg), sossettings, B_coeffs);


    Bc_star = value(B_coeffs);
    B1c_star = Bc_star(1:size(Bx_coeffs, 1));
    B2c_star = Bc_star(size(Bx_coeffs, 1)+1:end);

    
    t_total = round(toc(t_start));
    fprintf("Took %i s \n", t_total);
    
end

