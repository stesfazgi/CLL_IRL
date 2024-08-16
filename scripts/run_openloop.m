cd('..');
addpath('functions\');

% Dataset
n_demos = 7;
n_samples = 1000;

% Plotting
plot_padding = 0.5;
n_interpolation = 100;

% Output directory
output_dir = 'eval_output\eval_openloop\';

if ~exist(output_dir, "dir")
    error("Could not find output directory: %s", output_dir);
end

for shape_id=1:30
    
    % Prepare state variables
    yalmip('clear');
    x = sdpvar(2, 1);
    
    % Load dataset
    %%%%%%%%% Perpare dataset
    % Load dataset
    [data_pos, data_vel, shapename, dt] = plot_shape(shape_id, n_demos, n_samples, false, [0 0]);
    % Normalize dataset
    scale_factor = 1/max(abs(data_pos(:)));
    data_pos = data_pos*scale_factor;
    data_vel = data_vel*scale_factor;
    data_start_idx = linspace(n_samples, n_samples*n_demos, n_demos) - (n_samples - 1);
    % Plot settings
    axis_bounds = zeros(2, 2);
    axis_bounds(1, :) = [min(data_pos(1, :))-plot_padding max(data_pos(1, :))+plot_padding];
    axis_bounds(2, :) = [min(data_pos(2, :))-plot_padding max(data_pos(2, :))+plot_padding];

    
    % Load System
    savepath = fullfile("eval_sys_unstable", strcat(string(shape_id), ".mat"));
    if ~exist(savepath, "file")
        error("System not found!")
    else
        [f, f_coeffs] = func_restore_poly_sys(x, savepath);
        fprintf("Found system for shapeid %i.\n", shape_id);
    end
    
    % Create data for OL dynamics
    x1_range = linspace(axis_bounds(1, 1), axis_bounds(1, 2), n_interpolation);
    x2_range = linspace(axis_bounds(2, 1), axis_bounds(2, 2), n_interpolation);
    x1_dot_surf = zeros(n_interpolation);
    x2_dot_surf = zeros(n_interpolation);

    for i=1:n_interpolation
        for j=1:n_interpolation

            x_curr = [x1_range(i); x2_range(j)];
            x_dot_curr = replace(f, x, x_curr);

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
    
    algoname = 'ol';
    
    data_params.n_samples = n_samples;
    data_params.n_demos = n_demos;
    
    output_file = fullfile(output_dir, strcat('plot_data_', algoname, '_', shapename, '.mat'));
    
    save(output_file, 'shapename', 'algoname', 'data_params', ...
        'x1_dot_surf', 'x2_dot_surf', ...
        'x1_range', 'x2_range', 'target_data');
    
end