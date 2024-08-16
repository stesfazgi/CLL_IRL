%% Script to create eigendynamics that are partly fitted to the dynamics with an additive unstable part
addpath("..\functions");
yalmip('clear');
cd("..");


%% Settings
% Path to store output (relative to root dir of repo)
relative_path = "eval_sys_unstable";
% Maximum degree of f
degF = 3;
% Set coefficients with magnitude lower than fTol to 0
fTol = 1e-3;
% Scaling factor for linear unstable part
add_f = 1.5;
% Scaling factor for fitted polynomial
scale_f = 0.1;
% Number of samples per demo to use for fitting
n_samples = 1000;
% Number of demos to use for fitting
n_demos = 7;



full_path = fullfile(pwd, relative_path);

if ~exist(full_path, "dir")
    % Create folder
    mkdir(full_path);
    
    % Create struct that will later be parsed to JSON
    % Store generation parameters
    sys_struct.degF = degF;
    sys_struct.fTol = fTol;
    sys_struct.add_f = 1.5;
    sys_struct.scale_f = 0.1;
    sys_struct.n_samples = 1000;
    sys_struct.n_demos = 7;
    
    % Create system states
    x = sdpvar(2, 1);
    
    % Define replacement variables for state variables. Those will later
    % function as placeholders while calling "sdisplay" to avoid that
    % YALMIP replaces "x(1)" with "internal(-)". When parsing the string
    % with sympy, expressions without indexing are easier to handle.
    z1 = sdpvar(1, 1);
    z2 = sdpvar(1, 1);
    
    % Create linear, unstable, additive part
    add_sys = add_f*x;
    sys_struct.add_sys = sdisplay(add_sys);
    
    % Iterate over shapes
    for sid=1:30
        % Load dataset
        [data_pos, data_vel, shapename, ~] = plot_shape(sid, n_demos, n_samples, false);
        % Normalize dataset
        scale_factor = 1/max(abs(data_pos(:)));
        data_pos = data_pos*scale_factor;
        data_vel = data_vel*scale_factor;
        
        % Generate path for exporting current system
        curr_filename = strcat(string(sid), ".mat");
        curr_syspath = fullfile(full_path, curr_filename);
        
        % Fit system
        [f, f_coeffs] = func_fit_poly_sys(x, degF, data_vel, ...
            data_pos, fTol, add_sys, scale_f, curr_syspath);
        
        padding = 1;
        axis_bounds = zeros(2, 2);
        axis_bounds(1, :) = [min(data_pos(1, :))-padding max(data_pos(1, :))+padding];
        axis_bounds(2, :) = [min(data_pos(2, :))-padding max(data_pos(2, :))+padding];
        
        % Append data to output struct
        placeholder_f = replace(f, x, [z1; z2]);
        f_string = sdisplay(placeholder_f);
        
        % Sanity check 1
        if any(contains(f_string, "internal"))
            error("Found internal in string while processing shape %s.", shapename);
        end
        % Sanity check 2
        if any(contains(f_string, "x"))
            error("Found x in string while processing shape %s.", shapename);
        end
        
        sys_struct.(shapename).f1 = f_string{1};
        sys_struct.(shapename).f2 = f_string{2};
        sys_struct.(shapename).coeffs1 = f_coeffs(1, :);
        sys_struct.(shapename).coeffs2 = f_coeffs(2, :);
        sys_struct.(shapename).boundsx1 = axis_bounds(1, :);
        sys_struct.(shapename).boundsx2 = axis_bounds(2, :);
        
    end
    
    j_path = fullfile(full_path, "systems.json");
    j_string = jsonencode(sys_struct);
    
    fid = fopen(j_path, 'w');
    fprintf(fid, j_string);
    fclose(fid);
    
    
else
    fprintf("Output folder already exists. Not overwriting. \n");
end
