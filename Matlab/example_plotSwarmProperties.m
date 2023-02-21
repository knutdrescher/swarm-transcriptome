% This example script generates heatmaps for all three replicates for a
% specific list of swarm properties
clear;
close all

% Load data
result = load(fullfile('data', 'swarmProperties.mat'));
result = result.result;

% The variable mat contains time and radial position in the first and
% fourth column, respectively (second and third are xy coordinates),
% followed by the normalized gene expression levels
mat = result.mat;

% restrict plotting to sampled spaces
keep = ~cellfun(@isempty, {result.info.sample});
mat = mat(keep,:);

% Imaging was started approx. 1h before the start of expansion and
% the first sampling process, so we subtract that hour to make it
% consistent
time = mat(:,1)-1;
dist = mat(:,4);

% The variable exp_id contains the replicate id for each sample 
exp_id = result.exp_id(keep);

% propertyNames is the list of property names in the order as they appear in mat
property_names = result.propertyNames;

% Example properties to plot
propertiesToPlot = {'Cell speed', 'Nematic Order'};

% Find indices of genes
inds = zeros(length(propertiesToPlot),1);
for j = 1:length(propertiesToPlot)
        index_gene = find(cellfun(@(x) strcmp(x, propertiesToPlot{j}), property_names));
        if length(index_gene)==1
            inds(j) = index_gene+4;
        end
end

% Determine colorbar range of genes 
range = zeros(length(propertiesToPlot),2);
for j = 1:length(propertiesToPlot)
        index_gene = find(cellfun(@(x) strcmp(x, propertiesToPlot{j}), property_names));
        if length(index_gene)==1
            inds(j) = index_gene+4;
            range(j,:) = [prctile(mat(:,index_gene),1),prctile(mat(:,index_gene),99)];
        end
end

% Loop over each replicate and plot heatmap
for j = 1:3
    exp_rows = exp_id==j;
    % The matrix to be given to the plotting function needs to contain
    % time and radial position in its first two columns
    matToPlot = [time(exp_rows), dist(exp_rows), mat(exp_rows, inds)];

    % additional arguments here could be a path to save the heatmap  
    createHeatmap( matToPlot, 3:size(matToPlot,2), range, propertiesToPlot, j);

    % As a minimal example, the function can also be called without range
    % or gene names
    createHeatmap(matToPlot);
end