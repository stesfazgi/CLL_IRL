function fn_handle = func_sdpvar2fnhandle(f, x)
% func_sdpvar2fnhandle Convert a sdpvar expression to a function handle
%
%       IMPORTANT:  The input sdpvar has to only depend on the state x
%       IMPORTANT2: This only works with polynomial expressions
%
%   Input parameters:
%       f:          MxK SDPVAR representing a polynomial expression in x
%       x:          Nx1 SDPVAR representing the free variables in f
%
%   Output parameters:
%       fn_handle:  Handle of f that can be called with fn_handle(x)
%                   The argument of fn_handle is a vector with the size of x
    

    % create vector from input function and split into coefficients and
    % monomials
    try
        [fc, fm] = coefficients(f(:), x);
    catch
        % splitting fails if some entries in f are of degree 0 and some are
        % higher
        f_vec_tmp = f(:);
        % flag entries with degree 0
        deg_flags = zeros(size(f_vec_tmp));
        for i=1:size(f_vec_tmp, 1)
            deg_flags(i) = ~degree(f_vec_tmp(i));
        end
        
        % Temporarily raise degree for problematic entries
        f_vec_tmp = f_vec_tmp + deg_flags*x(1);
        
        % Rerun splitting
        [coeffs, fm] = coefficients(f_vec_tmp, x);
        
        % Revert raised degrees
        correction_coeffs = sparse(size(coeffs, 1), size(coeffs, 2));
        correction_coeffs(:, 1) = ones(size(correction_coeffs, 1), 1);
        fc = coeffs.*(~deg_flags) + correction_coeffs.*((deg_flags).*f(:));

    end
    
    % Manipulate monomial expression into compatible string
    sym_f_string = sdisplay(fm);
    sym_f_string = strrep(sym_f_string, "(", "");
    sym_f_string = strrep(sym_f_string, ")", "");
    
    % Create symbolic expression from string
    sym_f = str2sym(sym_f_string);

    % Output format is different if f maps to a scalar
    if size(f, 1)==1
        sym_f = fc'*sym_f;
    else
        sym_f = fc*sym_f;
    end
    
    
    % Find symbolic variables in new expression
    x_sym_new = symvar(sym_f).';

    % Create function handle from string
    fn_handle = matlabFunction(reshape(sym_f, size(f)), 'Vars', {x_sym_new});

    % Adjust indexing if input f is not dependend on all states 
    if size(x_sym_new, 1)~=size(x, 1)
        % Find all depended state
        fn_handle = @(y) fn_handle(y(degree(f, x) > 0));
    end
    
end

