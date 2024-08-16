function [pos, vel, shapename, t_sample] = plot_shape(shape_idx, n_demos, n_points, enablePlotting, skip)
% This is a matlab function illustrating 30 human handwriting motions
% recorded from Tablet-PC. These datas can be found in the folder
% 'DataSet'. 
%
% Please acknowledge the authors in any academic publications
% that have made use of this library by citing the following paper:
%
%  S. M. Khansari-Zadeh and A. Billard, "Learning Stable Non-Linear Dynamical 
%  Systems with Gaussian Mixture Models", IEEE Transaction on Robotics, 2011.
%
% To get latest upadate of the software please visit
%                          http://lasa.epfl.ch/khansari
%
% Please send your feedbacks or questions to:
%                           mohammad.khansari_at_epfl.ch


%%
names = {'Angle','BendedLine','CShape','DoubleBendedLine','GShape',...
         'heee','JShape','JShape_2','Khamesh','Leaf_1',...
         'Leaf_2','Line','LShape','NShape','PShape',...
         'RShape','Saeghe','Sharpc','Sine','Snake',...
         'Spoon','Sshape','Trapezoid','Worm','WShape','Zshape',...
         'Multi_Models_1','Multi_Models_2','Multi_Models_3','Multi_Models_4'};

[~, max_idx] = size(names);

if ~exist('enablePlotting', 'var')
    enablePlotting = true;
end

if ~exist('skip', 'var')
    skip = [0 0];
end

if (shape_idx > max_idx || shape_idx < 0)
    fprintf('ERROR: Maximum shape id is %i, but you selected %i', max_idx, shape_idx);
    return;
else
    pos = [];
    vel = [];
    
    %% preprocessing
    load(['LASA/DataSet/' names{shape_idx} '.mat'],'demos','dt') %loading the model
    
    shapename = names{shape_idx};

    
    if enablePlotting
        % plotting the result
        figure('name',names{shape_idx},'position',[100   100   600   800]);
        hold on;
    end
    
    
    max_demos = length(demos);
    if (n_demos > max_demos)
        fprintf('ERROR: Only %i demoes availabe (requested %i)', max_demos, n_demos);
        return
    end
    
    selected = round(linspace(1+skip(1), 1000-skip(2), n_points));
    
    for i=1:n_demos
%         pos = [pos, demos{draw_demos(i)}.pos(:, selected)]; % pos = [pos, demos{i}.pos(:, selected)]; %
%         vel = [vel, demos{draw_demos(i)}.vel(:, selected)]; % vel = [vel, demos{i}.vel(:, selected)]; %
        pos = [pos, demos{i}.pos(:, selected)];
        vel = [vel, demos{i}.vel(:, selected)];
    end
    
    if enablePlotting
        plot(pos(1, :), pos(2, :));
        quiver(pos(1,:),pos(2,:),vel(1,:),vel(2,:));
        hold on;

        xlabel('x (mm)','fontsize',15);
        ylabel('y (mm)','fontsize',15);
        title(names{shape_idx},'fontsize',15)
        hold off;
    end
    
    t_sample = dt;
    
end