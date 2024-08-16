clear;
clc;
yalmip('clear');
addpath('..\functions');
addpath('..\LASA\DataSet\');

%%%%%%%%%% Evaluation specific parameters
% Dataset
n_samples = 100;
n_demos = 7;

% Simulation/Plotting
plot_padding = 0.5;
n_interpolation = 100;
n_sim_steps = 2000;

% Output directory
algoname = 'sosclf';
output_dir = '..\eval_output\eval_sosclf_v8_b16_rV1e-4_rB1e-4\';

if ~exist(output_dir, "dir")
    mkdir(output_dir);
    fprintf("Making new directory: %s \n", output_dir);
else
    warning("Run already exists. Continuing...");
end

%%%%%%%%%% Algorithm specific parameters
degV = 8;
degB = 16;
n_iter = 5;
regV = 1e-4;
regB = 1e-4;

% eval_shape_ids = [1 3 19 23 26 30];
eval_shape_ids = linspace(1, 30, 30);

for esid=1:length(eval_shape_ids)
    
    shape_id = eval_shape_ids(esid);

    yalmip('clear');
    fprintf("Evaluating Shape %i \n", shape_id);
    
    %%%%%%%%% Perpare dataset
    % Load dataset
    [data_pos, data_vel, shapename, dt] = plot_shape(shape_id, n_demos, n_samples, false, [0 0]);
    % Normalize dataset
    scale_factor = 1/max(abs(data_pos(:)));
    data_pos = data_pos*scale_factor;
    data_vel = data_vel*scale_factor;
    data_start_idx = linspace(n_samples, n_samples*n_demos, n_demos) - (n_samples - 1);

    output_file = fullfile(output_dir, strcat('plot_data_', algoname, '_', shapename, '.mat'));
    if exist(output_file, "file")
        warning("Shape %s already exist. Skipping...", shapename);
        continue;
    end

    %%%%%%%%% Plot settings
    axis_bounds = zeros(2, 2);
    axis_bounds(1, :) = [min(data_pos(1, :))-plot_padding max(data_pos(1, :))+plot_padding];
    axis_bounds(2, :) = [min(data_pos(2, :))-plot_padding max(data_pos(2, :))+plot_padding];

    %%%%%%%%% Load unstable system dynamics
    x1 = sdpvar(1, 1);
    x2 = sdpvar(1, 1);
    % State vector
    x = [x1; x2];

    savepath = fullfile("..\eval_sys_unstable", strcat(string(shape_id), ".mat"));
    if ~exist(savepath, "file")
        error("System not found!")
    else
        [f, f_coeffs] = func_restore_poly_sys(x, savepath);
        fprintf("Using old system.\n");
    end

    % Input matrix of system
    g = eye(2);

    % Create options struct
    opt.degV = degV;
    opt.degB = degB;
    opt.gammaV = regV;
    opt.gammaB = regB;
    opt.n_iter = n_iter;
    opt.init_scale = 1;

    % Run algorithm
    t_start = tic;
    [Vc_star, B1c_star, B2c_star] = func_sosclf_gs_reg(x, f, g, opt, data_pos, data_vel);
    t_algo = toc(t_start);


    % Create expressions
    [V, V_coeff, ~] = polynomial(x, degV, 2);
    [B1, B1_coeff, ~] = polynomial(x, degB);
    [B2, B2_coeff, ~] = polynomial(x, degB);

    % Recreate closed loop system expression
    dVdx = jacobian(V, x);
    B = [B1, 0; 0 B2];
    f_cl = f - B*g*g'*dVdx';
    f_cl_star = replace(f_cl, [V_coeff; B1_coeff; B2_coeff], [Vc_star; B1c_star; B2c_star]);
    V_star = replace(V, V_coeff, Vc_star);

    f_cl_star_h = func_sdpvar2fnhandle(f_cl_star, x);
    V_star_h = func_sdpvar2fnhandle(V_star, x);

    % Create reproductions
    reproductions = {};
    initial_states = data_pos(:, data_start_idx);

    % Simulate CL system for each demo
    for i=1:n_demos
        % Start at initial state
        curr_x = initial_states(:, i);
        curr_traj = [curr_x];

        % Simulate for n_samples steps
        for j=1:n_sim_steps
            % Get system dynamics at current state
%             dx = replace(f_cl_star, x, curr_x);
            dx = f_cl_star_h(curr_x);
            % Discrete integration + update state
            curr_x = curr_x + dx*dt;
            % Save state in trajectory
            curr_traj = [curr_traj curr_x];
        end

        reproductions{end + 1} = curr_traj;
    end

    % Create data for CL dynamics and V contours
    x1_range = linspace(axis_bounds(1, 1), axis_bounds(1, 2), n_interpolation);
    x2_range = linspace(axis_bounds(2, 1), axis_bounds(2, 2), n_interpolation);
    v_surf = zeros(n_interpolation);
    x1_dot_surf = zeros(n_interpolation);
    x2_dot_surf = zeros(n_interpolation);

    for i=1:n_interpolation
        for j=1:n_interpolation

            x_curr = [x1_range(i); x2_range(j)];
%             x_dot_curr = replace(f_cl_star, x, x_curr);
%             v_surf(j, i) = replace(V_star, x, x_curr);
            x_dot_curr = f_cl_star_h(x_curr);
            v_surf(j, i) = V_star_h(x_curr);


            x1_dot_surf(j, i) = x_dot_curr(1);
            x2_dot_surf(j, i) = x_dot_curr(2);

        end
    end
    
    target_data = {};
    for i=1:n_demos
        if i==n_demos
            target_data{end+1} = data_pos(:, data_start_idx(i):end);
        else
            target_data{end+1} = data_pos(:, data_start_idx(i):data_start_idx(i+1)-1);
        end
        
    end

    algo_params.degV = degV;
    algo_params.degB = degB;
    algo_params.n_iter = n_iter;

    data_params.n_samples = n_samples;
    data_params.n_demos = n_demos;

    save(output_file, 'shapename', 'algoname', 'algo_params', 'data_params', ...
        'v_surf', 'x1_dot_surf', 'x2_dot_surf', 'reproductions', ...
        'x1_range', 'x2_range', 't_algo', 'target_data');
    
end


beep;

