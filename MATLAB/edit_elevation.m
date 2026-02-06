%% Edit Elevation - Apply Waypoint-Based Track Geometry
% This script loads a track geometry from generate_groundtrack.m and applies
% elevation changes along the track path defined by waypoints.
%
% Prerequisite: Run generate_groundtrack.m first to create track_geometry.mat
%
% Author: Tim Jusko
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
dataFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Raw');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');

% Elevation drop for the track (how deep to carve)
baseElevationDrop = -30;  % meters (will be combined with track z-profile)

%% Load track geometry
geometryFile = fullfile(outputFolder, 'track_geometry.mat');
if ~exist(geometryFile, 'file')
    error('Track geometry file not found. Run generate_groundtrack.m first.');
end
fprintf('Loading track geometry from: %s\n', geometryFile);
geomData = load(geometryFile);
geom = geomData.trackGeometry;

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

fprintf('\nTrack parameters (loaded from file):\n');
fprintf('  Seed: %d\n', geom.randomSeed);
fprintf('  Waypoints: %d\n', numWaypoints);
fprintf('  Total length: %.0f m\n', geom.totalLength);
fprintf('  Corridor width: %.0f m\n', corridorWidth);
fprintf('  Base elevation drop: %.0f m\n', baseElevationDrop);

%% Build dense track path from waypoints
% Interpolate between waypoints to get a dense set of points along the track
fprintf('\nBuilding dense track path...\n');

trackPoints = [];  % Will store [x, y, z] for dense points along track

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
        pz = wp1.z + t * (wp2.z - wp1.z);
        trackPoints = [trackPoints; px, py, pz];
    end
end

fprintf('  Dense track points: %d\n', size(trackPoints, 1));

%% Load all tiles
fprintf('\nLoading tiles...\n');
numTiles = length(tifBaseNames);
elevationData = cell(numTiles, 1);
geoInfo = cell(numTiles, 1);

for i = 1:numTiles
    baseName = tifBaseNames{i};
    
    % Load geotiffinfo
    geoInfoFile = fullfile(dataFolder, [baseName '_geotiffinfo.mat']);
    geoInfoData = load(geoInfoFile);
    geoInfoFields = fieldnames(geoInfoData);
    geoInfo{i} = geoInfoData.(geoInfoFields{1});
    
    % Load elevation data
    rasterFile = fullfile(dataFolder, [baseName '_readgeoraster.mat']);
    rasterData = load(rasterFile);
    elevationData{i} = rasterData.A;
    
    fprintf('  Loaded: %s (%d x %d)\n', baseName, size(elevationData{i}, 1), size(elevationData{i}, 2));
end

%% Process each tile
fprintf('\nProcessing tiles...\n');

for i = 1:numTiles
    baseName = tifBaseNames{i};
    fprintf('\n  Tile %d: %s\n', i, baseName);
    
    % Get tile bounds
    bbox = geoInfo{i}.BoundingBox;
    tileXMin = bbox(1, 1);
    tileXMax = bbox(2, 1);
    tileYMin = bbox(1, 2);
    tileYMax = bbox(2, 2);
    
    % Get pixel scale
    cellSize = geoInfo{i}.PixelScale(1);
    
    % Get dimensions
    [nRows, nCols] = size(elevationData{i});
    
    % Create coordinate grids for this tile
    xCoords = linspace(tileXMin + cellSize/2, tileXMax - cellSize/2, nCols);
    yCoords = linspace(tileYMax - cellSize/2, tileYMin + cellSize/2, nRows);
    [X, Y] = meshgrid(xCoords, yCoords);
    
    % Initialize mask and elevation adjustment arrays
    trackMask = false(nRows, nCols);
    elevAdjust = zeros(nRows, nCols);
    
    % For each pixel, check distance to track and calculate elevation adjustment
    fprintf('    Calculating track mask (this may take a moment)...\n');
    
    % Vectorized approach: for each track point, mark nearby pixels
    for tp = 1:size(trackPoints, 1)
        px = trackPoints(tp, 1);
        py = trackPoints(tp, 2);
        pz = trackPoints(tp, 3);
        
        % Skip track points outside this tile (with margin)
        if px < tileXMin - halfWidth || px > tileXMax + halfWidth || ...
           py < tileYMin - halfWidth || py > tileYMax + halfWidth
            continue;
        end
        
        % Find pixels within corridor width of this track point
        dist = sqrt((X - px).^2 + (Y - py).^2);
        nearbyMask = dist <= halfWidth;
        
        % Update mask
        trackMask = trackMask | nearbyMask;
        
        % For pixels in the track, set their elevation adjustment
        newPixels = nearbyMask & (elevAdjust == 0);
        elevAdjust(newPixels) = baseElevationDrop + pz;
    end
    
    % Count affected pixels
    numAffected = sum(trackMask(:));
    fprintf('    Pixels in track: %d (%.2f%%)\n', numAffected, 100 * numAffected / numel(trackMask));
    
    % Apply elevation change
    originalElevation = elevationData{i};
    modifiedElevation = originalElevation;
    modifiedElevation(trackMask) = originalElevation(trackMask) + elevAdjust(trackMask);
    
    % Store modified data
    elevationData{i} = modifiedElevation;
    
    fprintf('    Original elevation range: %.2f to %.2f\n', min(originalElevation(:)), max(originalElevation(:)));
    fprintf('    Modified elevation range: %.2f to %.2f\n', min(modifiedElevation(:)), max(modifiedElevation(:)));
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
figure('Name', 'Elevation Contours - Track Edit', 'Position', [200 200 800 600]);
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Elevation Contour Map - Track Edit (Seed: %d)', geom.randomSeed));
axis equal tight;
hold on;

% Overlay track centerline
trackX = [waypoints.x];
trackY = [waypoints.y];
plot(trackX, trackY, 'r-', 'LineWidth', 2, 'DisplayName', 'Track');
plot(waypoints(1).x, waypoints(1).y, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(waypoints(end).x, waypoints(end).y, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');
legend('Location', 'best');

fprintf('\n--- Edit Complete ---\n');
fprintf('Generated contour plot of edited elevation data.\n');
fprintf('Track carved: %.0f m wide, base drop %.0f m\n', corridorWidth, abs(baseElevationDrop));
fprintf('Run export_geotiff.m to generate TIF files for Trian3D.\n');
