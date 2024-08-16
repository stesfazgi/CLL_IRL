yalmip_path = 'PATH_TO_YALMIP_MASTER';
sedumi_path = 'PATH_TO_SEDUMI_MASTER';

% Install SeDuMi
addpath(sedumi_path)
install_sedumi


% Install YALMIP
addpath(yalmip_path)
addpath(join([yalmip_path, '\extras']))
addpath(join([yalmip_path, '\solvers']))
addpath(join([yalmip_path, '\modules']))
addpath(join([yalmip_path, '\modules\parametric']))
addpath(join([yalmip_path, '\modules\moment']))
addpath(join([yalmip_path, '\modules\global']))
addpath(join([yalmip_path, '\modules\sos']))
addpath(join([yalmip_path, '\operators']))