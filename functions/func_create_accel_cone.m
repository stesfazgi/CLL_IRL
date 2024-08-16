function [cone_cnstr, epivar, err_dataset] = func_create_accel_cone(x, f, g, B, dVdx, target_pos, target_vel)
%   func_create_accel_cone Accelerates SOC-Constraint creation via fn-handles
%
%   Input parameters:
%       x:              Nx1 SDPVAR representing the statevector
%       f:              Nx1 SDPVAR representing the open loop dynamics
%       g:              NxM SDPVAR representing the input matrix
%       B:              Mx1 SDPVAR representing the diagonal of the beta matrix
%                       (B=0.5 if a cone for V should be created)
%       dVdx:           1xN SDPVAR representing the gradient of V
%       target_pos:     NxK double containing K states for evaluation
%       target_vel:     NxK double containing K velocities for evaluation
%
%   Output parameters:
%       cone_cnstr:     SOC constraint
%       epivar:         SDPVAR representing the epigraph variable
%       err_dataset:    NKx1 SDPVAR containing the error expression
%                       (only used for debugging)

    t_start = tic;

    v_opt = false;

    % Figure out what we are optimizing over
    if isa(B, "double")
        v_opt = true;
    end

    [n_states, n_data] = size(target_pos);

    % Create handles from sdpvars
    f_h = func_sdpvar2fnhandle(f, x);
    g_h = func_sdpvar2fnhandle(g, x);

    % Evaluate f and g
    f_eval = zeros(n_states, n_data);
    g_eval = cell(n_data, 1);

    for i=1:n_data
        f_eval(:, i) = f_h(target_pos(:, i));
        g_eval{i} = g_h(target_pos(:, i));
    end
    
    err_dataset = sdpvar(n_states, n_data);
    

    if ~v_opt
        % We are optimizing over beta
        
        % Evaluate gradient
        dVdxT_h = func_sdpvar2fnhandle(dVdx', x);
        
        dVdxT_eval = zeros(n_states, n_data);
        for i=1:n_data
            dVdxT_eval(:, i) = dVdxT_h(target_pos(:, i));
        
            err_dataset(:, i) = target_vel(:, i) - ( f_eval(:, i) -g_eval{i}*replace(B, x, target_pos(:, i))*g_eval{i}'*dVdxT_eval(:, i) );
        end
        
    else
        % Optimization over V, beta is a scalar
        
        for i=1:n_data
            err_dataset(:, i) = target_vel(:, i) - ( f_eval(:, i) -g_eval{i}*B*g_eval{i}'*replace(dVdx', x, target_pos(:, i)));
        end
        
    end

    % Create epigraph variable
    epivar = sdpvar(1, 1);

    % Create cone constraint
    cone_cnstr = cone(err_dataset(:), epivar);

    t_end = toc(t_start);

end

