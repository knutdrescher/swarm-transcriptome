function info = calculateSwarmProperties( path )

pxSize = 0.398;
localDensity_rad = round(5/pxSize);
localNumberDensity_rad = round(10/pxSize);
nematicOrder_rad = round(10/pxSize);
swarm_rad = round(30/pxSize);

localDensityFrame = 120;

orientationDiffInSwarm = 15; % degrees

minVSwarm = 10/(100*pxSize); %10µm/s
maxVNonMotile = 6/(100*pxSize); %6µm/s

folders = dir(path);
folders = folders(3:end);
folders = folders([folders.isdir]);

info = struct.empty;

for k = 1:length(folders)
    files = dir(fullfile(path, folders(k).name));
    files = files(3:end);
    files = files([files.isdir]);
    r = getWidth(files);
    info = struct.empty;
    counter = 1;
    info(1).width = r;
    disp(folders(k).name);
    for f=files'
        sub = fullfile(path,folders(k).name, f.name);
        if ~strcmp(f.name,'videos')
            
            try
                pathToTif = fullfile(sub,'Untitled_1','Untitled_1_MMStack_Pos0.ome.tif');
                pathToLabels = fullfile(sub, 'labels_thresh0_3_whole.tif');
                pathToProb = fullfile(sub, 'prob_whole.tif');
                
                info(counter).path = pathToTif;
                [rad, x, y] = getPos(f.name);
                info(counter).x = x;
                info(counter).y = y;
                info(counter).radius = rad;
                disp(f.name)
                tic;
                
                imgInfo = imfinfo(pathToTif);
                labelInfo = imfinfo(pathToLabels);
                frames = size(imgInfo,1);
                frames = min(48,frames);
                
                % cell arrays to store single cell properties
                speeds = cell(frames-1,1);
                density_biomass = cell(frames-1, 1);
                density_nCell = cell(frames-1, 1);
                cellLengths = cell(frames-1,1);
                cellVolume = cell(frames-1,1);
                cellAspectRatio = cell(frames-1,1);
                nematic_order_orientations = cell(frames-1,1);
                nonMotClusterSizes = cell(frames-1,1);
                
                
                % vectors to store global properties
                numberOfCells = zeros(frames-1,1);
                biomassDensity = zeros(frames-1,1);
                biomassDensityFluctuations = zeros(frames-1,1);
                fractionOfNonMotilecells = zeros(frames-1,1);
                swarmingFact = zeros(frames-1,1);
                nNonMotClusters = zeros(frames-1,1);
                fractionOfSwarmingCells = zeros(frames-1,1);
                
                %start of calculations
                for j=1:frames-1
                    try
                    im1 = imread(pathToTif,'Index',j, 'Info', imgInfo);
                    im2 = imread(pathToTif,'Index',j+1, 'Info', imgInfo);
                    labels = imread(pathToLabels, 'Index',j, 'Info', labelInfo);
                    labels = labels';
                    cellnum = max(unique(labels));
                    mask = labels>0;
                    
                    % global fields
                    numberOfCells(j) = max(labels(:));
                    biomassDensity(j) = mean(mask(:));
                    densities = getLocalDensity( localDensityFrame, localDensityFrame, mask );
                    biomassDensityFluctuations(j) = std(densities);
                    
                    [ux,uy] = getFlowField(im1,im2);
                    speed_cells = sqrt(ux.^2+uy.^2);
                    vx = regionprops(labels, ux, 'MeanIntensity');
                    vy = regionprops(labels, uy, 'MeanIntensity');
                    
                    centroids = regionprops(labels, 'Centroid', 'Orientation','BoundingBox','Extent', 'MajorAxisLength', 'MinorAxisLength', 'Area');
                    
                    
                    cents = [centroids.Centroid];
                    cx = cents(1:2:end);
                    cy = cents(2:2:end);
                    
                    %determine potential neighbors based on knn
                    nNeighs = 2*(cellnum/(size(ux,1)*size(ux,2)))*4/3*pi*(max([localDensity_rad,swarm_rad,nematicOrder_rad, localNumberDensity_rad])^2);
                    neighIds = knnsearch([cy;cx]',[cy;cx]', 'K', ceil(nNeighs));
                    
                    % determine covered area in circle around each cell
                    px_neigh = arrayfun(@(x) max(round(x)-localDensity_rad,1):min(round(x)+localDensity_rad,size(im1,2)), cx, 'UniformOutput', false);
                    py_neigh = arrayfun(@(x) max(round(x)-localDensity_rad,1):min(round(x)+localDensity_rad,size(im1,1)), cy, 'UniformOutput', false);
                    comb = cellfun(@(x,y) combvec( x, y), py_neigh, px_neigh, 'UniformOutput', false);
                    inCircle = cellfun(@(x,y,z) x(:,hypot(x(1,:)-y,x(2,:)-z)<localDensity_rad), comb, num2cell(cy), num2cell(cx),'UniformOutput', false);
                    inCircle = cellfun(@(x) sub2ind(size(labels), x(1,:),x(2,:)), inCircle,'UniformOutput', false);
                    areaFractionCovered = cellfun(@(x) sum(labels(x)>0)/numel(x), inCircle);
                    
                    % determine speed per cell
                    props = regionprops(labels,speed_cells, 'MeanIntensity', 'MaxIntensity');
                    speedDist = [props.MeanIntensity];
        
                    % determine nematic order for each cell based on
                    % orientation or based on movement
                    S_orientation = zeros(length(cellnum),1);
                    S_motility = zeros(length(cellnum),1);
                    vx = [vx.MeanIntensity];
                    vy = [vy.MeanIntensity];
                    nNeighbors = zeros(length(cellnum),1);
                    orientations = [centroids.Orientation];
                    swarm_neighs = zeros(size(neighIds));
                    swarm_fact = zeros(length(cellnum),1);
                    for i = 1:cellnum
                        
                        potentialNeighbors = neighIds(i,:);
                        dists_potentialNeighbors = arrayfun(@(x,y) hypot(x-cx(i),y-cy(i)), cx(potentialNeighbors), cy(potentialNeighbors));

                        % determine number of neighbors
                        neigh = dists_potentialNeighbors < localNumberDensity_rad & dists_potentialNeighbors>0;
                        nNeighbors(i) = sum(neigh);
                        
                        % determine nematic order based on cells with
                        % centroid closer than given radius
                        orientation_diff = orientations(potentialNeighbors)-orientations(i);
                        orientation_diff = min([abs(orientation_diff);abs(180+orientation_diff); abs(180-orientation_diff)]);
                        
                        neigh = dists_potentialNeighbors < nematicOrder_rad & dists_potentialNeighbors > 0;
                        
                        S_orientation(i) = mean(1.5*cosd(orientation_diff(neigh)).^2-0.5);
                        
                        neigh = potentialNeighbors(neigh);
                        S_motility(i) = mean(arrayfun(@(x,y) 1.5*(dot([vx(i) vy(i)],[x y])/hypot(vx(i),vy(i))/hypot(x,y)).^2-0.5, vx(neigh), vy(neigh)));
                        
                        speed_potentialNeighbors = speedDist(potentialNeighbors);
                        
                        % determine swarm neighbors - use same distance as
                        % for nematic order but this can be changed
                        if speedDist(i)>minVSwarm
                            swarm_neighs(i,:) = (orientation_diff < orientationDiffInSwarm) & (dists_potentialNeighbors < swarm_rad) & (dists_potentialNeighbors>0) & (speed_potentialNeighbors>minVSwarm); 
                            swarm_fact(i) = sum(swarm_neighs(i,:))/sum((dists_potentialNeighbors < swarm_rad) & (dists_potentialNeighbors>0));
                        else
                            swarm_fact(i) = 0;
                        end
                    
                    end
                    
                    % Non-motile cells and cell clusters - be aware of
                    % faulty segmentation
                     [nonMotClusterSize, numberOfNonMotClusters, ~] = getNonMotClusterSize( [cx', cy'],  10, speedDist < maxVNonMotile);
                    % Swarm cluster detection based on clustering approach
                    angles = arrayfun(@angle, vx+1i*vy);
                    matForSwarmDetection = [cx', cy', angles', [centroids.Orientation]'];
                    nSwarming = estimateSwarms( matForSwarmDetection, nanmean(swarm_fact), 3, speedDist>minVSwarm);
                    

                    % save properties in cell array
                    speeds{j} = speedDist;
                    cellLengths{j} = [centroids.MajorAxisLength];
                    cellAspectRatio{j} = [centroids.MajorAxisLength]./[centroids.MinorAxisLength];
                    cellVolume{j} = [centroids.Area];
                    
                    density_nCell{j} = nNeighbors;
                    density_biomass{j} = areaFractionCovered;
                    
                    nematic_order_orientations{j} = S_orientation;
                    
                    swarmingFact(j) = nanmean(swarm_fact);
                    fractionOfSwarmingCells(j) = nSwarming/cellnum;

                     
                    fractionOfNonMotilecells(j) = sum(speedDist < maxVNonMotile)/numel(speedDist);
                    nNonMotClusters(j) = numberOfNonMotClusters;
                    nonMotClusterSizes{j} = nonMotClusterSize;

                   
                    catch err
                        disp(err.message);
                    end
                end
                
                 toc;
%               
                % basic parameters
                info(counter).speed = speeds;
                info(counter).cellLength = cellLengths;
                info(counter).cellVolume = cellVolume;
                info(counter).density_nCells = density_nCell;
                info(counter).density_biomass = density_biomass;
                info(counter).cellAspectRatio = cellAspectRatio;
                
                %emergent parameters
                info(counter).nematicOrder = nematic_order_orientations;
            
                info(counter).swarmingFact = swarmingFact;
                info(counter).fractionOfNonMotilecells = fractionOfNonMotilecells;
                info(counter).fractionOfSwarmingCells = fractionOfSwarmingCells;
                info(counter).nNonMotClust = nNonMotClusters;
                info(counter).nonMotClusterSize = nonMotClusterSizes;
                
                % global parameters
                info(counter).numberOfCells_global = numberOfCells;
                info(counter).biomassDensity_global = biomassDensity; 
                info(counter).biomassDensityFluctuations_global = biomassDensityFluctuations;

            catch err
                disp(err.message);
            end
            
            counter = counter+1;
        end
    end
    
    props = whos('info');


%     if props.bytes*10^(-9)>2
%         save(fullfile(path, folders(k).name,'results.mat'), 'info', '-v7.3');
%     else
%         save(fullfile(path, folders(k).name,'results.mat'), 'info');
%     end

    
end



    function [pos, x, y] = getPos(name)
        match = regexp(name, ...
            '(?<=^x)(?<pos_x>-?\d+)y(?<pos_y>-?\d+)(?=.*$)','names');
        
        x = str2double(match.pos_x);
        y = str2double(match.pos_y);
        pos = sqrt(x^2 + y^2);
    end

    function [rad,x,y] = getWidth(files)
        rad = 0;
        for fi=files'
            match = regexp(fi.name, ...
                '(?<=^x)(?<pos_x>-?\d+)y(?<pos_y>-?\d+)(?=.*$)','names');
            
            x = str2double(match.pos_x);
            y = str2double(match.pos_y);
            radN = sqrt(x^2 + y^2);
            if radN>rad
                rad = radN;
            end
        end
    end
    





end

