%% Plot Elevation Surface from GeoTIFF Files
% This script loads digital elevation model (DEM) data from pre-saved
% MATLAB objects (geotiffinfo and readgeoraster outputs) and visualizes
% them as a 3D surface plot.
%
% Note: Since the Mapping Toolbox is not available, we use pre-generated
% .mat files containing the geotiffinfo and readgeoraster outputs.
%
% Author: Tim Jusko
% Date: 2026-02-06

clear; clc; close all;

%% Load Project Configuration
config = load_project_config();
dataFolder = config.rawFolder;

fprintf('Project: %s\n', config.projectName);

% Define the base names of the TIF files (without extension)
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};

fprintf('Loading data for %d tiles:\n', length(tifBaseNames));
for i = 1:length(tifBaseNames)
    fprintf('  %d. %s\n', i, tifBaseNames{i});
end

%% Load pre-saved geotiffinfo and readgeoraster objects
elevationData = cell(length(tifBaseNames), 1);
geoInfo = cell(length(tifBaseNames), 1);

for i = 1:length(tifBaseNames)
    baseName = tifBaseNames{i};
    
    % Load geotiffinfo object
    geoInfoFile = fullfile(dataFolder, [baseName '_geotiffinfo.mat']);
    geoInfoData = load(geoInfoFile);
    % Get the variable name (first field in the struct)
    geoInfoFields = fieldnames(geoInfoData);
    geoInfo{i} = geoInfoData.(geoInfoFields{1});
    
    % Load readgeoraster object (contains elevation data as 'A')
    rasterFile = fullfile(dataFolder, [baseName '_readgeoraster.mat']);
    rasterData = load(rasterFile);
    elevationData{i} = rasterData.A;
    
    fprintf('\nFile: %s.tif\n', baseName);
    fprintf('  Size: %d x %d pixels\n', size(elevationData{i}, 1), size(elevationData{i}, 2));
    fprintf('  BoundingBox: [%.0f, %.0f] to [%.0f, %.0f]\n', ...
        geoInfo{i}.BoundingBox(1,1), geoInfo{i}.BoundingBox(1,2), ...
        geoInfo{i}.BoundingBox(2,1), geoInfo{i}.BoundingBox(2,2));
    fprintf('  Elevation range: %.2f to %.2f\n', ...
        min(elevationData{i}(:)), max(elevationData{i}(:)));
end

%% Merge elevation data into a single surface
% Determine the spatial extent of all tiles using BoundingBox from geotiffinfo
% BoundingBox format: [xMin yMin; xMax yMax]
numTiles = length(tifBaseNames);
allXLimits = zeros(numTiles, 2);
allYLimits = zeros(numTiles, 2);

for i = 1:numTiles
    bbox = geoInfo{i}.BoundingBox;
    allXLimits(i, :) = [bbox(1,1), bbox(2,1)];  % [xMin, xMax]
    allYLimits(i, :) = [bbox(1,2), bbox(2,2)];  % [yMin, yMax]
end

% Calculate combined extent
xMin = min(allXLimits(:, 1));
xMax = max(allXLimits(:, 2));
yMin = min(allYLimits(:, 1));
yMax = max(allYLimits(:, 2));

fprintf('\n--- Combined Extent ---\n');
fprintf('X: %.2f to %.2f (%.2f m)\n', xMin, xMax, xMax - xMin);
fprintf('Y: %.2f to %.2f (%.2f m)\n', yMin, yMax, yMax - yMin);

% Get the resolution from the first file (PixelScale)
cellSize = geoInfo{1}.PixelScale(1);  % PixelScale is [scaleX, scaleY, scaleZ]
fprintf('Cell size: %.2f m\n', cellSize);

% Create the merged grid
nCols = round((xMax - xMin) / cellSize);
nRows = round((yMax - yMin) / cellSize);

mergedElevation = NaN(nRows, nCols);

% Place each tile in the merged grid
for i = 1:numTiles
    % Get bounding box for this tile
    bbox = geoInfo{i}.BoundingBox;
    tileXMin = bbox(1,1);
    tileYMax = bbox(2,2);
    
    % Calculate the position of this tile in the merged grid
    colStart = round((tileXMin - xMin) / cellSize) + 1;
    rowStart = round((yMax - tileYMax) / cellSize) + 1;
    
    [tileRows, tileCols] = size(elevationData{i});
    
    rowEnd = rowStart + tileRows - 1;
    colEnd = colStart + tileCols - 1;
    
    % Ensure indices are within bounds
    rowEnd = min(rowEnd, nRows);
    colEnd = min(colEnd, nCols);
    
    mergedElevation(rowStart:rowEnd, colStart:colEnd) = ...
        elevationData{i}(1:(rowEnd-rowStart+1), 1:(colEnd-colStart+1));
end

%% Create coordinate grids for plotting
[X, Y] = meshgrid(linspace(xMin, xMax, nCols), linspace(yMax, yMin, nRows));

%% Contour plot
figure('Name', 'Elevation Contours', 'Position', [200 200 800 600]);
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title('Elevation Contour Map');
axis equal tight;

fprintf('\n--- Visualization Complete ---\n');
fprintf('Generated contour plot of elevation data.\n');
fprintf('\nNext step: Run edit_canyon.m to apply terrain modifications.\n');