function [finished] = func_run_sosclf_gs_reg_from_config(cfg)
%FUNC_RUN_SOSCLF_GS_REG_FROM_CONFIG Run a single evaluation cycle based on
%a given configuration, used for gridsearching

    if ~exist(cfg.root, "dir")
        mkdir(cfg.root);
    end

    for gs_shape_index=1:length(cfg.shapeids)
        
        % Extract settings from config
        plot_padding = cfg.plot_padding;
        n_interpolation = cfg.n_interpolation;
        n_sim_steps = cfg.n_sim_steps;

        shape_id = cfg.shapeids(gs_shape_index);
        degV = cfg.degV;
        degB = cfg.degB;
        n_demos = cfg.n_demos;
        n_samples = cfg.n_samples;
        n_iter = cfg.n_iter;
        gammaV = cfg.regV;
        gammaB = cfg.regB;

        algoname = 'sosclf';

        yalmip("clear");


        %%%%%%%%% Perpare dataset
        % Load dataset
        [data_pos, data_vel, shapename, dt] = plot_shape(shape_id, n_demos, n_samples, false, [0 0]);
        % Normalize dataset
        scale_factor = 1/max(abs(data_pos(:)));
        data_pos = data_pos*scale_factor;
        data_vel = data_vel*scale_factor;
        data_start_idx = linspace(n_samples, n_samples*n_demos, n_demos) - (n_samples - 1);

        output_file = fullfile(cfg.root, strcat('plot_data_', algoname, '_', shapename, '.mat'));
    
        if exist(output_file, "file")
            fprintf("Shape already completed. Skipping...")
            continue;
        end
    
        %%%%%%%%% Plot settings
        axis_bounds = zeros(2, 2);
        axis_bounds(1, :) = [min(data_pos(1, :))-plot_padding max(data_pos(1, :))+plot_padding];
        axis_bounds(2, :) = [min(data_pos(2, :))-plot_padding max(data_pos(2, :))+plot_padding];
    
        %%%%%%%%% Load system dynamics
        x1 = sdpvar(1, 1);
        x2 = sdpvar(1, 1);
        % State vector
        x = [x1; x2];
    
        if cfg.eigendyns
            savepath = fullfile("..\eval_sys_unstable", strcat(string(shape_id), ".mat"));
            if ~exist(savepath, "file")
                error("System not found!")
            else
                [f, f_coeffs] = func_restore_poly_sys(x, savepath);
                fprintf("Using old system.\n");
            end
        else
            f = [0; 0];
        end
    
        % Input matrix of system
        g = eye(2);

        % Create options struct
        opt.degV = degV;
        opt.degB = degB;
        opt.gammaV = gammaV;
        opt.gammaB = gammaB;
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
                target_data{end+1} = data_pos(:, data_start_idx(i):data_start_idx(i+1));
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

    finished = true;

    save(fullfile(cfg.root, "config.mat"), "cfg");

end

