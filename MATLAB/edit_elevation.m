%% Edit Elevation - Shape Terrain to Match Course Profile
% This script modifies terrain elevation along the flight path so that
% a pilot flying at constant AGL will experience the climbs and descents
% defined in the course structure.
%
% How it works:
%   - Climb event: terrain DROPS so pilot must climb to maintain AGL
%   - Descent event: terrain RISES so pilot must descend to maintain AGL
%   - Straight event: terrain stays level (constant slope = 0)
%
% The terrain is smoothly blended at corridor edges to look natural.
%
% Prerequisite: Run generate_groundtrack.m first to create track_geometry.mat
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc; close all;

%% Load Project Configuration
config = load_project_config();
dataFolder = config.rawFolder;
outputFolder = config.editedFolder;

fprintf('Project: %s\n', config.projectName);

% Reference AGL altitude (the "floor" of the corridor)
% Terrain will be shaped so this clearance exists along the track
referenceAGL = 15;  % meters - pilot flies this high above reshaped terrain

% Edge blending width (smooth transition from modified to original terrain)
blendWidth = 100;  % meters - how wide the smooth transition zone is

%% Load track geometry
geometryFile = fullfile(outputFolder, 'track_geometry.mat');
if ~exist(geometryFile, 'file')
    error('Track geometry file not found. Run generate_groundtrack.m first.');
end
fprintf('Loading track geometry from: %s\n', geometryFile);
geomData = load(geometryFile);
geom = geomData.trackGeometry;

%% Load course structure for transition parameters
structureFile = fullfile(outputFolder, 'event_structure.mat');
if ~exist(structureFile, 'file')
    error('Course structure file not found. Run define_course_structure.m first.');
end
structData = load(structureFile);
transition = structData.eventStructure.transition;

% Extract geometry parameters
waypoints = geom.waypoints;
corridorWidth = geom.corridorWidth;
tifBaseNames = geom.tifBaseNames;
xMinTotal = geom.xMinTerrain;
xMaxTotal = geom.xMaxTerrain;
yMinTotal = geom.yMinTerrain;
yMaxTotal = geom.yMaxTerrain;

halfWidth = corridorWidth / 2;
numWaypoints = length(waypoints);

fprintf('\nTrack parameters:\n');
fprintf('  Seed: %d\n', geom.randomSeed);
fprintf('  Waypoints: %d\n', numWaypoints);
fprintf('  Total length: %.0f m\n', geom.totalLength);
fprintf('  Corridor width: %.0f m\n', corridorWidth);
fprintf('  Reference AGL: %.0f m\n', referenceAGL);
fprintf('  Blend width: %.0f m\n', blendWidth);

%% Build dense track path with elevation profile
% The waypoint z-values define the RELATIVE elevation profile:
%   z=0 at start, z increases during climb, z decreases during descent
% We invert this for terrain: terrain = base - z (so climb = terrain drops)

fprintf('\nBuilding dense track path with elevation profile...\n');

trackPoints = [];  % [x, y, targetTerrainZ, distanceAlongTrack, isTransition]
cumulativeDistance = 0;

for i = 1:(numWaypoints - 1)
    wp1 = waypoints(i);
    wp2 = waypoints(i + 1);
    
    % Distance between waypoints
    segmentLength = sqrt((wp2.x - wp1.x)^2 + (wp2.y - wp1.y)^2);
    
    % Number of points (1 per meter for high resolution)
    numPoints = max(2, ceil(segmentLength));
    
    % Interpolate
    for j = 0:(numPoints - 1)
        t = j / (numPoints - 1);
        px = wp1.x + t * (wp2.x - wp1.x);
        py = wp1.y + t * (wp2.y - wp1.y);
        pz = wp1.z + t * (wp2.z - wp1.z);  % Relative flight profile elevation
        
        dist = cumulativeDistance + t * segmentLength;
        trackPoints = [trackPoints; px, py, pz, dist, 0];  % 0 = not transition
    end
    
    cumulativeDistance = cumulativeDistance + segmentLength;
end

fprintf('  Dense track points: %d\n', size(trackPoints, 1));
fprintf('  Flight profile range: %.1f to %.1f m (relative)\n', ...
    min(trackPoints(:,3)), max(trackPoints(:,3)));

%% Store end of flight path info for transition calculation
flightPathEndX = trackPoints(end, 1);
flightPathEndY = trackPoints(end, 2);
flightPathEndZ = trackPoints(end, 3);  % Relative elevation at end of flight path
flightPathEndDist = trackPoints(end, 4);

% Calculate heading at end of track (from last two waypoints)
lastWP = waypoints(end);
secondLastWP = waypoints(end-1);
endHeading = atan2d(lastWP.y - secondLastWP.y, lastWP.x - secondLastWP.x);

%% Load all tiles and get reference elevation at track start
fprintf('\nLoading tiles...\n');
numTiles = length(tifBaseNames);
elevationData = cell(numTiles, 1);
geoInfo = cell(numTiles, 1);

for i = 1:numTiles
    baseName = tifBaseNames{i};
    
    geoInfoFile = fullfile(dataFolder, [baseName '_geotiffinfo.mat']);
    geoInfoData = load(geoInfoFile);
    geoInfoFields = fieldnames(geoInfoData);
    geoInfo{i} = geoInfoData.(geoInfoFields{1});
    
    rasterFile = fullfile(dataFolder, [baseName '_readgeoraster.mat']);
    rasterData = load(rasterFile);
    elevationData{i} = rasterData.A;
    
    fprintf('  Loaded: %s (%d x %d)\n', baseName, size(elevationData{i}, 1), size(elevationData{i}, 2));
end

% Get terrain elevation at start point
startX = waypoints(1).x;
startY = waypoints(1).y;
startTerrainElev = getTerrainElevation(startX, startY, elevationData, geoInfo);
fprintf('\nStart point terrain elevation: %.1f m\n', startTerrainElev);

% Base terrain level for the corridor (terrain at start minus reference AGL headroom)
% Actually, we want the terrain to be UNDER the flight path
% Flight path starts at z=0 (relative), terrain should be referenceAGL below that
baseTerrainLevel = startTerrainElev;
fprintf('Base terrain level at start: %.1f m\n', baseTerrainLevel);

%% Add transition segment to blend back to original terrain at exit
% This prevents cliffs/walls where the modified corridor meets original terrain
fprintf('\n--- Adding Transition Segment ---\n');

% Get original terrain elevation at projected exit point
projectedEndX = flightPathEndX + transition.minLength * cosd(endHeading);
projectedEndY = flightPathEndY + transition.minLength * sind(endHeading);
endTerrainElev = getTerrainElevation(projectedEndX, projectedEndY, elevationData, geoInfo);

% Calculate required z to match original terrain:
% Modified terrain = baseTerrainLevel + flightProfileZ
% We need: baseTerrainLevel + targetZ = endTerrainElev
% So: targetZ = endTerrainElev - baseTerrainLevel
targetEndZ = endTerrainElev - baseTerrainLevel;
requiredZChange = targetEndZ - flightPathEndZ;

fprintf('  Flight path ends at z=%.1f m (relative)\n', flightPathEndZ);
fprintf('  Original terrain at exit: %.1f m\n', endTerrainElev);
fprintf('  Required z change: %.1f m to match original terrain\n', requiredZChange);

% Calculate transition length based on gradient limit
minLengthNeeded = abs(requiredZChange) / (transition.maxGradient / 100);
transitionLength = max(transition.minLength, minLengthNeeded);
actualGradient = abs(requiredZChange) / transitionLength * 100;

fprintf('  Transition length: %.0f m (%.1f%% gradient)\n', transitionLength, actualGradient);

% Generate transition track points
% In the transition, we blend from the modified corridor elevation to original terrain
% We store the transition progress (0 to 1) in column 5
numTransitionPoints = max(2, ceil(transitionLength));
for j = 1:numTransitionPoints
    t = j / numTransitionPoints;  % 0 to 1 (exclusive of 0, inclusive of 1)
    
    px = flightPathEndX + t * transitionLength * cosd(endHeading);
    py = flightPathEndY + t * transitionLength * sind(endHeading);
    % Store the transition blend factor (t) - will blend from modified to original terrain
    % pz stores the flight profile z at the END of the main corridor (for reference)
    pz = flightPathEndZ;  % Reference z at corridor end
    dist = flightPathEndDist + t * transitionLength;
    
    trackPoints = [trackPoints; px, py, pz, dist, t];  % t = transition blend factor
end

fprintf('  Transition points added: %d\n', numTransitionPoints);
fprintf('  Final z: %.1f m (should match original terrain offset: %.1f m)\n', ...
    trackPoints(end, 3), targetEndZ);

%% Get ORIGINAL terrain elevation along entire track (before any modification)
fprintf('\nSampling original terrain elevation along track...\n');
originalTerrainAlongTrack = zeros(size(trackPoints, 1), 1);
for tp = 1:size(trackPoints, 1)
    originalTerrainAlongTrack(tp) = getTerrainElevation(trackPoints(tp,1), trackPoints(tp,2), elevationData, geoInfo);
end

%% Process each tile - reshape terrain along track
fprintf('\nReshaping terrain along track...\n');

for i = 1:numTiles
    baseName = tifBaseNames{i};
    fprintf('\n  Tile %d: %s\n', i, baseName);
    
    % Get tile bounds
    bbox = geoInfo{i}.BoundingBox;
    tileXMin = bbox(1, 1);
    tileXMax = bbox(2, 1);
    tileYMin = bbox(1, 2);
    tileYMax = bbox(2, 2);
    cellSize = geoInfo{i}.PixelScale(1);
    [nRows, nCols] = size(elevationData{i});
    
    % Create coordinate grids for this tile
    xCoords = linspace(tileXMin + cellSize/2, tileXMax - cellSize/2, nCols);
    yCoords = linspace(tileYMax - cellSize/2, tileYMin + cellSize/2, nRows);
    [X, Y] = meshgrid(xCoords, yCoords);
    
    % Original elevation
    originalElevation = elevationData{i};
    modifiedElevation = originalElevation;
    
    % For each pixel, find closest track point and calculate target elevation
    fprintf('    Calculating terrain modifications...\n');
    
    % Build KD-tree-like structure for fast lookup (simplified: just track XY)
    trackXY = trackPoints(:, 1:2);
    
    % Process in chunks for memory efficiency
    for row = 1:nRows
        for col = 1:nCols
            px = X(row, col);
            py = Y(row, col);
            
            % Quick bounds check - skip if far from track
            if px < min(trackXY(:,1)) - halfWidth - blendWidth || ...
               px > max(trackXY(:,1)) + halfWidth + blendWidth || ...
               py < min(trackXY(:,2)) - halfWidth - blendWidth || ...
               py > max(trackXY(:,2)) + halfWidth + blendWidth
                continue;
            end
            
            % Find distance to nearest track point
            distances = sqrt((trackXY(:,1) - px).^2 + (trackXY(:,2) - py).^2);
            [minDist, nearestIdx] = min(distances);
            
            % Only modify pixels within corridor + blend zone
            if minDist > halfWidth + blendWidth
                continue;
            end
            
            % Get the flight profile elevation at this track point
            flightProfileZ = trackPoints(nearestIdx, 3);
            isTransition = trackPoints(nearestIdx, 5);  % 0 = corridor, >0 = transition blend factor
            
            % Calculate target terrain elevation
            if isTransition == 0
                % Inside main corridor: use flight profile
                % - At start (flightProfileZ=0): terrain = baseTerrainLevel
                % - During climb (flightProfileZ>0): terrain RISES (so pilot climbs to follow)
                % - During descent (flightProfileZ<0): terrain DROPS (so pilot descends to follow)
                targetTerrainElev = baseTerrainLevel + flightProfileZ;
            else
                % In transition zone: blend from corridor elevation to original terrain
                % At t=0 (start of transition): use corridor elevation
                % At t=1 (end of transition): use original terrain
                corridorTerrainElev = baseTerrainLevel + flightProfileZ;
                transitionBlend = isTransition;  % 0 to 1
                % Smoothly blend from corridor to original using cosine interpolation
                smoothBlend = 0.5 * (1 - cos(pi * transitionBlend));
                targetTerrainElev = corridorTerrainElev * (1 - smoothBlend) + origElev * smoothBlend;
            end
            
            % Calculate blend factor based on distance from track center
            if minDist <= halfWidth
                % Inside corridor: full modification
                blendFactor = 1.0;
            else
                % In blend zone: smooth transition
                blendDist = minDist - halfWidth;
                blendFactor = 1.0 - (blendDist / blendWidth);
                blendFactor = max(0, min(1, blendFactor));
                % Smooth the blend with cosine interpolation
                blendFactor = 0.5 * (1 + cos(pi * (1 - blendFactor)));
            end
            
            % Apply blended elevation
            origElev = originalElevation(row, col);
            modifiedElevation(row, col) = origElev + blendFactor * (targetTerrainElev - origElev);
        end
    end
    
    % Store modified data
    elevationData{i} = modifiedElevation;
    
    % Statistics
    changed = abs(modifiedElevation - originalElevation) > 0.01;
    numChanged = sum(changed(:));
    if numChanged > 0
        maxChange = max(abs(modifiedElevation(changed) - originalElevation(changed)));
        fprintf('    Pixels modified: %d (%.2f%%)\n', numChanged, 100 * numChanged / numel(changed));
        fprintf('    Max elevation change: %.1f m\n', maxChange);
    else
        fprintf('    No pixels modified in this tile\n');
    end
end

%% Save edited data
fprintf('\n--- Saving Edited Data ---\n');

for i = 1:numTiles
    baseName = tifBaseNames{i};
    
    % Save geotiffinfo as *_geotiffinfo_edited.mat
    geoInfoEdited = geoInfo{i};
    geoInfoFile = fullfile(outputFolder, [baseName '_geotiffinfo_edited.mat']);
    save(geoInfoFile, 'geoInfoEdited');
    fprintf('  Saved: %s_geotiffinfo_edited.mat\n', baseName);
    
    % Save elevation data as *_readgeoraster_edited.mat
    A = elevationData{i};
    R = [];
    rasterFile = fullfile(outputFolder, [baseName '_readgeoraster_edited.mat']);
    save(rasterFile, 'A', 'R');
    fprintf('  Saved: %s_readgeoraster_edited.mat\n', baseName);
end

%% Visualize the result
fprintf('\n--- Generating Visualization ---\n');

% Merge tiles for visualization
cellSize = geoInfo{1}.PixelScale(1);
nColsTotal = round((xMaxTotal - xMinTotal) / cellSize);
nRowsTotal = round((yMaxTotal - yMinTotal) / cellSize);
mergedElevation = NaN(nRowsTotal, nColsTotal);

for i = 1:numTiles
    bbox = geoInfo{i}.BoundingBox;
    tileXMin = bbox(1,1);
    tileYMax = bbox(2,2);
    
    colStart = round((tileXMin - xMinTotal) / cellSize) + 1;
    rowStart = round((yMaxTotal - tileYMax) / cellSize) + 1;
    
    [tileRows, tileCols] = size(elevationData{i});
    rowEnd = min(rowStart + tileRows - 1, nRowsTotal);
    colEnd = min(colStart + tileCols - 1, nColsTotal);
    
    mergedElevation(rowStart:rowEnd, colStart:colEnd) = ...
        elevationData{i}(1:(rowEnd-rowStart+1), 1:(colEnd-colStart+1));
end

% Create coordinate grids
[X, Y] = meshgrid(linspace(xMinTotal, xMaxTotal, nColsTotal), ...
                   linspace(yMaxTotal, yMinTotal, nRowsTotal));

%% Contour plot
figure('Name', 'Terrain Profile - Course Applied', 'Position', [200 200 900 700]);
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Reshaped Terrain - Course Profile Applied (Seed: %d)', geom.randomSeed));
axis equal tight;
hold on;

% Overlay track with event colors
eventColors = struct(...
    'start', [0.2 0.8 0.2], ...
    'straight', [0.2 0.4 1.0], ...
    'turn', [1.0 0.0 1.0], ...
    'climb', [1.0 0.5 0.0], ...
    'descent', [0.0 0.8 0.8], ...         % Cyan
    'transition_up', [0.5 0.5 0.5], ...   % Dark Gray
    'transition_down', [0.5 0.5 0.5]);    % Dark Gray

for i = 2:length(waypoints)
    wp1 = waypoints(i-1);
    wp2 = waypoints(i);
    segType = wp2.type;
    if isfield(eventColors, segType)
        c = eventColors.(segType);
    else
        c = [0.5 0.5 0.5];
    end
    % Make transition segments thicker and dashed for visibility
    if contains(segType, 'transition')
        plot([wp1.x, wp2.x], [wp1.y, wp2.y], '--', 'LineWidth', 5, 'Color', c);
    else
        plot([wp1.x, wp2.x], [wp1.y, wp2.y], '-', 'LineWidth', 3, 'Color', c);
    end
end

plot(waypoints(1).x, waypoints(1).y, 'o', 'MarkerSize', 14, 'MarkerFaceColor', [0 0.8 0], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 2);
plot(waypoints(end).x, waypoints(end).y, 's', 'MarkerSize', 14, 'MarkerFaceColor', [0.8 0 0], ...
    'MarkerEdgeColor', 'k', 'LineWidth', 2);

% Plot transition segment (from flight path end to actual corridor end)
transitionStartIdx = find(trackPoints(:,4) > flightPathEndDist, 1);
if ~isempty(transitionStartIdx)
    transitionX = trackPoints(transitionStartIdx:end, 1);
    transitionY = trackPoints(transitionStartIdx:end, 2);
    plot(transitionX, transitionY, '--', 'LineWidth', 4, 'Color', [0.5 0.5 0.5]);
    % Mark the true corridor end
    plot(trackPoints(end,1), trackPoints(end,2), 'd', 'MarkerSize', 12, ...
        'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerEdgeColor', 'k', 'LineWidth', 2);
end

%% Plot elevation profile along track
figure('Name', 'Track Elevation Profile', 'Position', [250 150 900 500]);

% Get EDITED terrain elevation along track
editedTerrainAlongTrack = zeros(size(trackPoints, 1), 1);
for tp = 1:size(trackPoints, 1)
    editedTerrainAlongTrack(tp) = getTerrainElevation(trackPoints(tp,1), trackPoints(tp,2), elevationData, geoInfo);
end

% Calculate expected flight path elevation
% - In corridor: baseTerrainLevel + flightProfileZ
% - In transition: blend from corridor elevation to original terrain
flightPathElevation = zeros(size(trackPoints, 1), 1);
for tp = 1:size(trackPoints, 1)
    isTransition = trackPoints(tp, 5);
    flightProfileZ = trackPoints(tp, 3);
    
    if isTransition == 0
        % In corridor: use flight profile
        flightPathElevation(tp) = baseTerrainLevel + flightProfileZ;
    else
        % In transition: blend to original terrain
        corridorElev = baseTerrainLevel + flightProfileZ;
        origElev = originalTerrainAlongTrack(tp);
        smoothBlend = 0.5 * (1 - cos(pi * isTransition));
        flightPathElevation(tp) = corridorElev * (1 - smoothBlend) + origElev * smoothBlend;
    end
end

% Plot all three
hold on;
h1 = plot(trackPoints(:,4), originalTerrainAlongTrack, 'b-', 'LineWidth', 2);
h2 = plot(trackPoints(:,4), editedTerrainAlongTrack, 'r-', 'LineWidth', 2);
h3 = plot(trackPoints(:,4), flightPathElevation, 'g--', 'LineWidth', 2);

% Mark transition start
xline(flightPathEndDist, '--', 'Transition', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');

xlabel('Distance along track (m)');
ylabel('Elevation (m)');
title('Elevation Profile Along Track');
legend([h1 h2 h3], 'Original Terrain', 'Edited Terrain', 'Flight Path (relative)', 'Location', 'best');
grid on;
hold off;

% Print verification
fprintf('\n--- Transition Verification ---\n');
fprintf('At end of corridor:\n');
fprintf('  Original terrain: %.2f m\n', originalTerrainAlongTrack(end));
fprintf('  Edited terrain:   %.2f m\n', editedTerrainAlongTrack(end));
fprintf('  Difference:       %.2f m (should be ~0)\n', abs(originalTerrainAlongTrack(end) - editedTerrainAlongTrack(end)));

%% 3D Surface Plot
figure('Name', '3D Terrain View', 'Position', [300 100 900 700]);

% Downsample for faster rendering (every 5th point)
downsample = 5;
Xd = X(1:downsample:end, 1:downsample:end);
Yd = Y(1:downsample:end, 1:downsample:end);
Zd = mergedElevation(1:downsample:end, 1:downsample:end);

% Plot surface
surf(Xd, Yd, Zd, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
zlabel('Elevation (m)');
title(sprintf('3D Terrain with Course Profile (Seed: %d)', geom.randomSeed));
axis tight;
hold on;

% Plot track as 3D line (at terrain surface + small offset for visibility)
trackZ = zeros(length(waypoints), 1);
for i = 1:length(waypoints)
    trackZ(i) = getTerrainElevation(waypoints(i).x, waypoints(i).y, elevationData, geoInfo) + 2;
end
plot3([waypoints.x], [waypoints.y], trackZ, 'r-', 'LineWidth', 3);

% Plot start and end markers
plot3(waypoints(1).x, waypoints(1).y, trackZ(1) + 5, 'go', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
plot3(waypoints(end).x, waypoints(end).y, trackZ(end) + 5, 'rs', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

% Set view angle
view(45, 30);
lighting gouraud;
camlight('headlight');

fprintf('\n--- Edit Complete ---\n');
fprintf('Terrain reshaped to match course profile.\n');
fprintf('  Corridor width: %.0f m\n', corridorWidth);
fprintf('  Flight profile range: %.1f to %.1f m\n', min(trackPoints(:,3)), max(trackPoints(:,3)));
fprintf('\nRun export_geotiff.m to generate TIF files for Trian3D.\n');

%% ========================================================================
%  HELPER FUNCTION: Get terrain elevation at a point
%  ========================================================================
function elev = getTerrainElevation(x, y, elevationData, geoInfo)
    elev = NaN;
    for i = 1:length(geoInfo)
        bbox = geoInfo{i}.BoundingBox;
        if x >= bbox(1,1) && x <= bbox(2,1) && y >= bbox(1,2) && y <= bbox(2,2)
            % Point is in this tile
            cellSize = geoInfo{i}.PixelScale(1);
            [nRows, nCols] = size(elevationData{i});
            
            col = round((x - bbox(1,1)) / cellSize) + 1;
            row = round((bbox(2,2) - y) / cellSize) + 1;
            
            col = max(1, min(nCols, col));
            row = max(1, min(nRows, row));
            
            elev = elevationData{i}(row, col);
            return;
        end
    end
end
