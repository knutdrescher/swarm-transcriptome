function createHeatmap( data_matrix, varargin)
%This function creates a heatmap plot from spatiotemporal transcriptome
%data or swarm properties. The first input is the matrix containing the
%data to be plotted, whose first two columns need to be time and radial
%position. 
% The second input is the index of the experiment

% Additional properties may be 
% - the indices of matrix columns to plot
% - the experiment index (used to modify save names and crop heatmaps)
% - gene or property names
% - a path where heatmaps will be saved as .png or .svg files
% - the colorbar range


times = data_matrix(:,1);
timepoints = unique(times);
pos = data_matrix(:,2);

%default settings
indices = 3:size(data_matrix,2);
range = [];
param = [];
exp_index = 1;
pathToSave = [];



switch nargin
    case 2
        indices = varargin{1};
    case 3
        indices = varargin{1};
        range = varargin{2};  
    case 4
        indices = varargin{1};
        range = varargin{2};
        param = varargin{3}; 
    case 5
        indices = varargin{1};
        range = varargin{2};
        param = varargin{3}; 
        exp_index = varargin{4}; 
    case 6
        indices = varargin{1};
        range = varargin{2};
        param = varargin{3}; 
        exp_index = varargin{4}; 
        pathToSave = varargin{5};
end

data_matrix = data_matrix(:,indices);

maxValues = zeros(size(data_matrix,2),1);
minValues = zeros(size(data_matrix,2),1);
allValues = zeros(size(data_matrix,2),1);


c = size(allValues,2);
for j=1:size(data_matrix,2)
    vec = data_matrix(:,j);
    allValues(j,c:c+length(vec)-1) = vec;
    maxval = max(vec);
    minval = min(vec);
    if maxval>maxValues(j)
        maxValues(j) = maxval;
    end
    if minval<minValues(j)
        minValues(j) = minval;
    end
end


map = parula(1000);
figures = gobjects(size(data_matrix,2)+1,1);


for f = 1:size(data_matrix,2)
    if f==1 || f>1
        h = figure;
        figures(f) = h;
        if ~isempty(param)
            title(param{f});
        end
        xlabel('Time (h)');
        ylabel('Position (mm)');

        box('on');

        ylim([-0.1 35])
        switch exp_index
            case 1
                xlim([0 6]);
            case 2
                xlim([-0.2 5.35]);
            case 3
                xlim([0.35 5.9]);
        end
        
        set(gca, 'FontSize', 18, 'LineWidth', 1.5)
        
        if f<=size(data_matrix,2)
            minz = prctile(allValues(f,1:end),5);
            maxz = prctile(allValues(f,1:end),95);
            minValues(f) = minz;
            maxValues(f) = maxz;
            set(gca,'colormap', map);
            if minz<maxz
                set(gca, 'CLim', [minz, maxz]);
            end
            if ~isempty(range)
                set(gca, 'CLim', [range(f,1), range(f,2)]);
            end
            c = colorbar;
            ylabel(c,'log_2(L_{RNA})'); 
            set(c, 'FontSize', 24, 'LineWidth', 1)

        end
    end
end
widths = zeros(length(timepoints),1);


for j = 1:length(timepoints)
    disp(j);
    t = timepoints(j);
    ind = times==t;
    positions = pos(ind);
    widths(j) = max(positions);
    
    [sorted, indices] = sort(positions);
    [t,height] = getHeightAndMiddle_t(timepoints, j);
    
    vals = data_matrix(ind, :);
    
    for l = 1:size(data_matrix,2)
        figure(figures(l));
        
        for n = 1:length(positions)
            [p,width] = getHeightAndMiddle(sorted,n);
            value = vals(indices(n),l);

            color = getColor(value, minValues(l), maxValues(l), map);

            rectangle(gca, 'Position', [t, p/1000, height, width/1000], 'Facecolor', color,'EdgeColor', color);
        end
    end
    
    
end


if ~isempty(pathToSave)
    if ~exist(pathToSave, 'dir')
        mkdir(pathToSave);
    end
for i = 1:numel(param)
    h_fig = figures(i);
    figure(h_fig);
    

    saveName = [param{i}, sprintf('_rep%d',exp_index)];
    print(h_fig, '-dpng','-r200', fullfile(pathToSave,  [saveName, '.png']));
    print(h_fig, '-dsvg','-r200', fullfile(pathToSave,  [saveName, '.svg']));
    delete(h_fig);
end
end



    function  [m,w] = getHeightAndMiddle(vec, index)

            if index ==1 && length(vec)>1
                w = vec(index+1);
                m = 0;
            elseif index==1 && length(vec) == 1
                w = 1500;
                m = -750;
            elseif index == length(vec) && length(vec)>1
                w = vec(index)-vec(index-1);
                m = vec(index)-w/2;
            else
                w = (vec(index+1) - vec(index-1))/2;
                m = 1/4*(vec(index-1) + vec(index+1)+2* vec(index))-w/2;
            end

    end

    function  [m,w] = getHeightAndMiddle_t(vec, index)

            if index ==1 && length(vec)>1
                w = vec(index+1)-vec(index);
                m = vec(index)-w/2;
            elseif index==1 && length(vec) == 1
                w = 1500;
                m = -750;
            elseif index == length(vec) && length(vec)>1
                w = vec(index)-vec(index-1);
                m = vec(index)-w/2;
            else
                w = (vec(index+1) - vec(index-1))/2;
                m = 1/4*(vec(index-1) + vec(index+1)+2* vec(index))-w/2;
            end

    end

    function  col = getColor(value, mini, maxi, map)
            if isnan(value)|| ~isreal(value)
                col = [1 1 1];
                return;
            end
            x = (value-mini)/(maxi-mini);
            x = min(x,0.999);
            x = max(x,0);
            x = x*1000;
            x = round(x);
            try
                col = map(x+1,:);
            catch
                disp('');
            end


    end


end

