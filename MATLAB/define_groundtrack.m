%% Define Ground Track Geometry
% This script defines the geometry of a racetrack/canyon path and saves
% it to a .mat file for use by other scripts (edit_elevation.m, generate_trees_kml.m)
%
% The geometry is defined as a straight diagonal line from bottom-left 
% to top-right of the terrain extent.
%
% Author: Tim Jusko
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
dataFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Raw');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output folder: %s\n', outputFolder);
end

%% Track parameters (modify these as needed)
trackWidth = 50;        % Width of the track in meters
elevationDrop = -100;   % Elevation change in meters (negative = dig down)

%% Define the base names of the TIF files
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};

%% Load tile info to determine terrain extent
fprintf('Loading tile information...\n');
numTiles = length(tifBaseNames);
geoInfo = cell(numTiles, 1);
elevationData = cell(numTiles, 1);

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
    
    fprintf('  Loaded: %s\n', baseName);
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

fprintf('\nTerrain extent:\n');
fprintf('  X: %.0f to %.0f (%.0f m)\n', xMinTotal, xMaxTotal, xMaxTotal - xMinTotal);
fprintf('  Y: %.0f to %.0f (%.0f m)\n', yMinTotal, yMaxTotal, yMaxTotal - yMinTotal);

%% Define the track path (bottom-left to top-right diagonal)
% Start point: bottom-left corner
trackStartX = xMinTotal;
trackStartY = yMinTotal;

% End point: top-right corner
trackEndX = xMaxTotal;
trackEndY = yMaxTotal;

% Calculate derived geometry values
dx = trackEndX - trackStartX;
dy = trackEndY - trackStartY;
trackLength = sqrt(dx^2 + dy^2);

% Normalized direction vector
dirX = dx / trackLength;
dirY = dy / trackLength;

% Perpendicular vector (pointing left of travel direction)
perpX = -dirY;
perpY = dirX;

halfWidth = trackWidth / 2;

% Line coefficients for distance calculation: a*x + b*y + c = 0
lineA = dy;
lineB = -dx;
lineC = dx * trackStartY - dy * trackStartX;
lineNormFactor = sqrt(lineA^2 + lineB^2);

fprintf('\nTrack geometry:\n');
fprintf('  Start: (%.0f, %.0f)\n', trackStartX, trackStartY);
fprintf('  End:   (%.0f, %.0f)\n', trackEndX, trackEndY);
fprintf('  Length: %.0f m\n', trackLength);
fprintf('  Width: %.0f m\n', trackWidth);
fprintf('  Elevation drop: %.0f m\n', elevationDrop);

%% Save track geometry to .mat file
trackGeometry.trackWidth = trackWidth;
trackGeometry.elevationDrop = elevationDrop;
trackGeometry.trackStartX = trackStartX;
trackGeometry.trackStartY = trackStartY;
trackGeometry.trackEndX = trackEndX;
trackGeometry.trackEndY = trackEndY;
trackGeometry.xMinTotal = xMinTotal;
trackGeometry.xMaxTotal = xMaxTotal;
trackGeometry.yMinTotal = yMinTotal;
trackGeometry.yMaxTotal = yMaxTotal;
trackGeometry.dx = dx;
trackGeometry.dy = dy;
trackGeometry.trackLength = trackLength;
trackGeometry.dirX = dirX;
trackGeometry.dirY = dirY;
trackGeometry.perpX = perpX;
trackGeometry.perpY = perpY;
trackGeometry.halfWidth = halfWidth;
trackGeometry.lineA = lineA;
trackGeometry.lineB = lineB;
trackGeometry.lineC = lineC;
trackGeometry.lineNormFactor = lineNormFactor;
trackGeometry.tifBaseNames = tifBaseNames;

geometryFile = fullfile(outputFolder, 'track_geometry.mat');
save(geometryFile, 'trackGeometry');
fprintf('\n--- Geometry Saved ---\n');
fprintf('Saved: %s\n', geometryFile);

%% Visualize: Contour plot with track overlay
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
figure('Name', 'Elevation Contours - Ground Track Definition', 'Position', [200 200 800 600]);
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title('Elevation Contour Map - Ground Track Definition');
axis equal tight;
hold on;

% Draw track boundaries
leftEdgeX = [trackStartX + halfWidth*perpX, trackEndX + halfWidth*perpX];
leftEdgeY = [trackStartY + halfWidth*perpY, trackEndY + halfWidth*perpY];
rightEdgeX = [trackStartX - halfWidth*perpX, trackEndX - halfWidth*perpX];
rightEdgeY = [trackStartY - halfWidth*perpY, trackEndY - halfWidth*perpY];

plot(leftEdgeX, leftEdgeY, 'r-', 'LineWidth', 2, 'DisplayName', 'Track Boundary');
plot(rightEdgeX, rightEdgeY, 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');

% Draw centerline
plot([trackStartX, trackEndX], [trackStartY, trackEndY], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Centerline');

% Mark start and end points
plot(trackStartX, trackStartY, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(trackEndX, trackEndY, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

legend('Location', 'best');

fprintf('\nNext step: Run edit_elevation.m to apply elevation changes.\n');
