function [f_poly, f_coeffs, savepath] = func_fit_poly_sys(x, degF, ref_dyn, ref_pos, tol, add_f, scale_f, filename)

    should_save = exist('filename', 'var');

    % Halt execution if any coeffs are unreasonably large
    max_tol = 1e2;

    % add_f is used to specify an additional dynamic on top of the
    % data-based dynamic to fit
    if ~exist('add_f', 'var')
        add_f = zeros(size(x));
    end


    % Get number of states and data
    n_states = size(x, 1);
    n_data = size(ref_pos, 2);
    
    % Create placeholder for symbolic polynomial system
    f_poly = sdpvar(n_states, 1);
    
    % Create symbolic expression for polynomial system dynamics (1 state)
    [fx, fx_coeff, ~] = polynomial(x, degF, 1);
    
    % Get gradient (symobilc row of regressor matrix A)
    dFxdc = jacobian(fx, fx_coeff);
    
    % Get number of polynomial coefficients
    n_coeff = size(fx_coeff, 1);
    
    % Create matrix to hold all fitted coefficients
    f_coeffs = zeros(n_states, n_coeff);
    
    
    % Fit to each state-dynamic independently
    for j=1:n_states
        
        % Target vector
        y = ref_dyn(j, :)';
        % Regressor matrix
        A = zeros(n_data, n_coeff);
        
        % Evaluate polynomial to fill regressor matrix
        for i=1:n_data
            A(i, :) = replace(dFxdc, x, ref_pos(:, i));
            y(i) = y(i) + replace(add_f(j), x, ref_pos(:, i));
        end
        
        % Create objective (=least squares)
        obj = norm(y - A*fx_coeff, 2);
        
        % Call solver, disable logging
        res = optimize([], obj, sdpsettings('verbose', 0));
        
        % Output warning if something went wrong
        if res.problem ~= 0
            fprintf('Encountered problem while fitting dynamics for state %i: %i', j, res.problem);
        end
        
        % Scale all coefficients and set coeffs under the tolerance to hard zero.
        fx_result = value(fx_coeff) * scale_f;
        fx_result = fx_result.*(abs(fx_result) > tol);
        
        % Update fitted dynamics for current state
        f_poly(j) = replace(fx, fx_coeff, fx_result);
        
        % Update matrix of coefficients
        f_coeffs(j, :) = fx_result;
        
    end
    
    n_large_coeffs = sum(f_coeffs(:) > max_tol);
    if n_large_coeffs > 0
        error('Encountered %i coefficients larger %i! Stopping execution', n_large_coeffs, max_tol);
    end
    
    if should_save
        save(filename, 'f_coeffs', 'degF', 'n_states');
        savepath = filename;
    end
    
    
end
