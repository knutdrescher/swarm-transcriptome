% This example script generates heatmaps for all three replicates for a
% specific list of genes
clear;
close all

% Load data
result = load(fullfile('data', 'transcriptomeData_normalized.mat'));
result = result.result;

% The variable mat contains time and radial position in the first two
% columns, followed by the normalized gene expression levels
mat = result.mat;
% Imaging was started approx. 1h before the start of expansion and
% the first sampling process, so we subtract that hour to make it
% consistent
time = mat(:,1)-1;
dist = mat(:,2);

% The variable exp_id contains the replicate id for each sample 
exp_id = result.exp_id;

% geneNames is the list of gene names in the order as they appear in mat
gene_names = result.geneNames;

% Example genes to plot
genesToPlot = {'srfAA', 'sigD'};

% Find indices of genes
inds = zeros(length(genesToPlot),1);
for j = 1:length(genesToPlot)
        index_gene = find(cellfun(@(x) strcmp(x, genesToPlot{j}), gene_names));
        if length(index_gene)==1
            inds(j) = index_gene+4;
        end
end

% Determine colorbar range of genes 
range = zeros(length(genesToPlot),2);
for j = 1:length(genesToPlot)
        index_gene = find(cellfun(@(x) strcmp(x, genesToPlot{j}), gene_names));
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
    createHeatmap( matToPlot, 3:size(matToPlot,2), range, genesToPlot, j);

    % As a minimal example, the function can also be called without range
    % or gene names
    createHeatmap(matToPlot);
end