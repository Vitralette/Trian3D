%% Create Canyon Edit - Diagonal Racetrack
% This script creates a canyon by reducing elevation along a diagonal path
% from bottom-left to top-right across both tiles.
%
% Parameters:
%   - Track width: 50m
%   - Elevation reduction: -100m
%
% Author: Generated for Trian3D Sample Project
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
dataFolder = fullfile('..', 'TRIAN3D', 'SampleProject');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');

% Canyon parameters
trackWidth = 50;        % Width of the canyon in meters
elevationDrop = -100;   % Elevation change in meters (negative = dig down)

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output folder: %s\n', outputFolder);
end

% Define the base names of the TIF files
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};

%% Load all tiles
fprintf('Loading tiles...\n');
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

%% Determine combined extent
allXLimits = zeros(numTiles, 2);
allYLimits = zeros(numTiles, 2);

for i = 1:numTiles
    bbox = geoInfo{i}.BoundingBox;
    allXLimits(i, :) = [bbox(1,1), bbox(2,1)];  % [xMin, xMax]
    allYLimits(i, :) = [bbox(1,2), bbox(2,2)];  % [yMin, yMax]
end

xMinTotal = min(allXLimits(:, 1));
xMaxTotal = max(allXLimits(:, 2));
yMinTotal = min(allYLimits(:, 1));
yMaxTotal = max(allYLimits(:, 2));

fprintf('\nCombined extent:\n');
fprintf('  X: %.0f to %.0f (%.0f m)\n', xMinTotal, xMaxTotal, xMaxTotal - xMinTotal);
fprintf('  Y: %.0f to %.0f (%.0f m)\n', yMinTotal, yMaxTotal, yMaxTotal - yMinTotal);

%% Define the diagonal track (bottom-left to top-right)
% Start point: bottom-left corner
trackStartX = xMinTotal;
trackStartY = yMinTotal;

% End point: top-right corner
trackEndX = xMaxTotal;
trackEndY = yMaxTotal;

fprintf('\nCanyon track:\n');
fprintf('  Start: (%.0f, %.0f)\n', trackStartX, trackStartY);
fprintf('  End:   (%.0f, %.0f)\n', trackEndX, trackEndY);
fprintf('  Width: %.0f m\n', trackWidth);
fprintf('  Elevation drop: %.0f m\n', elevationDrop);

% Calculate line parameters for distance calculation
dx = trackEndX - trackStartX;
dy = trackEndY - trackStartY;
lineLength = sqrt(dx^2 + dy^2);

% Line coefficients: a*x + b*y + c = 0
a = dy;
b = -dx;
c = dx * trackStartY - dy * trackStartX;
normFactor = sqrt(a^2 + b^2);

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
    % Note: In GeoTIFF, row 1 is typically at yMax (top)
    xCoords = linspace(tileXMin + cellSize/2, tileXMax - cellSize/2, nCols);
    yCoords = linspace(tileYMax - cellSize/2, tileYMin + cellSize/2, nRows);
    [X, Y] = meshgrid(xCoords, yCoords);
    
    % Calculate distance from each pixel to the track centerline
    distanceToLine = abs(a * X + b * Y + c) / normFactor;
    
    % Project point onto line and check if projection is within segment
    t = ((X - trackStartX) * dx + (Y - trackStartY) * dy) / (lineLength^2);
    
    % Create mask: within track width AND within segment
    halfWidth = trackWidth / 2;
    withinWidth = distanceToLine <= halfWidth;
    withinSegment = (t >= -0.01) & (t <= 1.01);
    
    canyonMask = withinWidth & withinSegment;
    
    % Count affected pixels
    numAffected = sum(canyonMask(:));
    fprintf('    Pixels in canyon: %d (%.2f%%)\n', numAffected, 100 * numAffected / numel(canyonMask));
    
    % Apply elevation change
    originalElevation = elevationData{i};
    modifiedElevation = originalElevation;
    modifiedElevation(canyonMask) = originalElevation(canyonMask) + elevationDrop;
    
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

% Plot 1: 3D Surface
figure('Name', 'Canyon Edit Result', 'Position', [100 100 1400 600]);

subplot(1, 2, 1);
surf(X, Y, mergedElevation, 'EdgeColor', 'none');
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
zlabel('Elevation (m)');
title('3D View - Canyon from Bottom-Left to Top-Right');
axis tight;
view(45, 30);
lighting gouraud;
camlight('headlight');

% Draw the track line on top
hold on;
zLine = max(mergedElevation(:)) + 50;
plot3([trackStartX trackEndX], [trackStartY trackEndY], [zLine zLine], ...
    'r-', 'LineWidth', 3);
legend('Terrain', 'Canyon Centerline', 'Location', 'best');

% Plot 2: Top-down view
subplot(1, 2, 2);
imagesc([xMinTotal xMaxTotal], [yMaxTotal yMinTotal], mergedElevation);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title('Top View with Canyon Track');
axis xy;
axis equal tight;
hold on;

% Draw track boundaries
halfWidth = trackWidth / 2;
theta = atan2(dy, dx);
perpX = halfWidth * sin(theta);
perpY = halfWidth * cos(theta);

% Left edge
plot([trackStartX - perpX, trackEndX - perpX], ...
     [trackStartY + perpY, trackEndY + perpY], 'r-', 'LineWidth', 2);
% Right edge  
plot([trackStartX + perpX, trackEndX + perpX], ...
     [trackStartY - perpY, trackEndY - perpY], 'r-', 'LineWidth', 2);
% Centerline
plot([trackStartX trackEndX], [trackStartY trackEndY], 'r--', 'LineWidth', 1);

fprintf('\n--- Edit Complete ---\n');
fprintf('Canyon created: %.0f m wide, %.0f m deep\n', trackWidth, abs(elevationDrop));
fprintf('Run export_geotiff.m to generate TIF files for Trian3D.\n');
