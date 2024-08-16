function [f_poly, f_coeffs] = func_restore_poly_sys(x, savepath)

    load(savepath, 'degF', 'f_coeffs', 'n_states');
    
    f_poly = sdpvar(n_states, 1);
    
    [fx, fx_coeff, ~] = polynomial(x, degF, 1);
    
    % Get number of polynomial coefficients
    n_coeff = size(fx_coeff, 1);
    
    
    for i=1:n_states
        f_poly(i) = replace(fx, fx_coeff, f_coeffs(i, :));
    end
    
end

