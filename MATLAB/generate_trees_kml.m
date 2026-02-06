%% Generate KML with Trees Along Canyon Boundaries
% This script generates a KML file with trees placed along the left and right
% boundaries of the canyon racetrack created by edit_canyon.m
%
% Trees are placed at equidistant intervals (default: 20m)
%
% Author: Generated for Trian3D Sample Project
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
editedFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Modified');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Tree placement parameters
treeSpacing = 20;       % Distance between trees in meters
trackWidth = 50;        % Must match edit_canyon.m
treeName = 'trees';     % Object name recognized by Trian3D (from KML bible)

% Define the base names of the TIF files
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};

%% Load edited tile info to get extent
fprintf('Loading tile information...\n');
numTiles = length(tifBaseNames);
geoInfo = cell(numTiles, 1);
elevationData = cell(numTiles, 1);

for i = 1:numTiles
    baseName = tifBaseNames{i};
    
    % Load from edited folder
    geoInfoFile = fullfile(editedFolder, [baseName '_geotiffinfo_edited.mat']);
    if ~exist(geoInfoFile, 'file')
        % Fall back to original if edited doesn't exist
        geoInfoFile = fullfile('..', 'TRIAN3D', 'SampleProject', [baseName '_geotiffinfo.mat']);
    end
    geoInfoData = load(geoInfoFile);
    geoInfoFields = fieldnames(geoInfoData);
    geoInfo{i} = geoInfoData.(geoInfoFields{1});
    
    % Load elevation data
    rasterFile = fullfile(editedFolder, [baseName '_readgeoraster_edited.mat']);
    if ~exist(rasterFile, 'file')
        rasterFile = fullfile('..', 'TRIAN3D', 'SampleProject', [baseName '_readgeoraster.mat']);
    end
    rasterData = load(rasterFile);
    elevationData{i} = rasterData.A;
end

%% Determine combined extent (same as edit_canyon.m)
allXLimits = zeros(numTiles, 2);
allYLimits = zeros(numTiles, 2);

for i = 1:numTiles
    bbox = geoInfo{i}.BoundingBox;
    allXLimits(i, :) = [bbox(1,1), bbox(2,1)];
    allYLimits(i, :) = [bbox(1,2), bbox(2,2)];
end

xMinTotal = min(allXLimits(:, 1));
xMaxTotal = max(allXLimits(:, 2));
yMinTotal = min(allYLimits(:, 1));
yMaxTotal = max(allYLimits(:, 2));

fprintf('Combined extent (UTM):\n');
fprintf('  X: %.0f to %.0f\n', xMinTotal, xMaxTotal);
fprintf('  Y: %.0f to %.0f\n', yMinTotal, yMaxTotal);

%% Define the diagonal track (same as edit_canyon.m)
trackStartX = xMinTotal;
trackStartY = yMinTotal;
trackEndX = xMaxTotal;
trackEndY = yMaxTotal;

% Track direction vector
dx = trackEndX - trackStartX;
dy = trackEndY - trackStartY;
trackLength = sqrt(dx^2 + dy^2);

% Normalized direction
dirX = dx / trackLength;
dirY = dy / trackLength;

% Perpendicular vector (pointing left of travel direction)
perpX = -dirY;
perpY = dirX;

halfWidth = trackWidth / 2;

fprintf('\nTrack parameters:\n');
fprintf('  Length: %.0f m\n', trackLength);
fprintf('  Width: %.0f m\n', trackWidth);

%% Calculate tree positions along both boundaries
numTrees = floor(trackLength / treeSpacing) + 1;
fprintf('\nGenerating %d trees per side (%d total)\n', numTrees, numTrees * 2);

% Arrays to store tree positions (UTM coordinates)
leftTreesUTM = zeros(numTrees, 2);
rightTreesUTM = zeros(numTrees, 2);

for i = 1:numTrees
    % Distance along the track
    dist = (i - 1) * treeSpacing;
    
    % Position along centerline
    centerX = trackStartX + dist * dirX;
    centerY = trackStartY + dist * dirY;
    
    % Left boundary (offset perpendicular)
    leftTreesUTM(i, 1) = centerX + halfWidth * perpX;
    leftTreesUTM(i, 2) = centerY + halfWidth * perpY;
    
    % Right boundary (offset perpendicular other direction)
    rightTreesUTM(i, 1) = centerX - halfWidth * perpX;
    rightTreesUTM(i, 2) = centerY - halfWidth * perpY;
end

%% Convert UTM to WGS84 (Latitude/Longitude)
% The data is in EPSG:25832 (UTM Zone 32N)
% KML requires WGS84 coordinates (longitude, latitude)

fprintf('\nConverting UTM to WGS84...\n');

% UTM Zone 32N parameters
utmZone = 32;
hemisphere = 'N';

% Convert all tree positions
leftTreesWGS84 = zeros(numTrees, 2);  % [lon, lat]
rightTreesWGS84 = zeros(numTrees, 2);

for i = 1:numTrees
    [leftTreesWGS84(i, 2), leftTreesWGS84(i, 1)] = utm2wgs84(leftTreesUTM(i, 1), leftTreesUTM(i, 2), utmZone, hemisphere);
    [rightTreesWGS84(i, 2), rightTreesWGS84(i, 1)] = utm2wgs84(rightTreesUTM(i, 1), rightTreesUTM(i, 2), utmZone, hemisphere);
end

fprintf('  Sample tree position (left, first): %.6f, %.6f (lon, lat)\n', leftTreesWGS84(1, 1), leftTreesWGS84(1, 2));
fprintf('  Sample tree position (right, last): %.6f, %.6f (lon, lat)\n', rightTreesWGS84(end, 1), rightTreesWGS84(end, 2));

%% Generate KML file
kmlFile = fullfile(outputFolder, 'canyon_trees.kml');
fprintf('\nWriting KML file: %s\n', kmlFile);

fid = fopen(kmlFile, 'w');

% KML header
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid, '  <Document>\n');

% Write left boundary trees
fprintf(fid, '    <!-- Left boundary trees -->\n');
for i = 1:numTrees
    fprintf(fid, '    <Placemark>\n');
    fprintf(fid, '      <name>%s</name>\n', treeName);
    fprintf(fid, '      <Point>\n');
    fprintf(fid, '        <coordinates>\n');
    fprintf(fid, '          %.14f,%.14f,0\n', leftTreesWGS84(i, 1), leftTreesWGS84(i, 2));
    fprintf(fid, '        </coordinates>\n');
    fprintf(fid, '      </Point>\n');
    fprintf(fid, '    </Placemark>\n');
end

% Write right boundary trees
fprintf(fid, '    <!-- Right boundary trees -->\n');
for i = 1:numTrees
    fprintf(fid, '    <Placemark>\n');
    fprintf(fid, '      <name>%s</name>\n', treeName);
    fprintf(fid, '      <Point>\n');
    fprintf(fid, '        <coordinates>\n');
    fprintf(fid, '          %.14f,%.14f,0\n', rightTreesWGS84(i, 1), rightTreesWGS84(i, 2));
    fprintf(fid, '        </coordinates>\n');
    fprintf(fid, '      </Point>\n');
    fprintf(fid, '    </Placemark>\n');
end

% KML footer
fprintf(fid, '  </Document>\n');
fprintf(fid, '</kml>\n');

fclose(fid);
fprintf('KML file written with %d trees.\n', numTrees * 2);

%% Visualize: Contour plot with tree markers
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

% Create figure
figure('Name', 'Canyon with Trees', 'Position', [100 100 1200 900]);

% Plot 1: Contour with trees (top-down)
subplot(1, 2, 1);
contourf(linspace(xMinTotal, xMaxTotal, nColsTotal), ...
         linspace(yMinTotal, yMaxTotal, nRowsTotal), ...
         flipud(mergedElevation), 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
hold on;

% Plot tree positions (UTM coordinates)
plot(leftTreesUTM(:, 1), leftTreesUTM(:, 2), 'g^', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
plot(rightTreesUTM(:, 1), rightTreesUTM(:, 2), 'g^', 'MarkerSize', 6, 'MarkerFaceColor', 'g');

% Draw track boundaries
plot([trackStartX + halfWidth*perpX, trackEndX + halfWidth*perpX], ...
     [trackStartY + halfWidth*perpY, trackEndY + halfWidth*perpY], 'r-', 'LineWidth', 2);
plot([trackStartX - halfWidth*perpX, trackEndX - halfWidth*perpX], ...
     [trackStartY - halfWidth*perpY, trackEndY - halfWidth*perpY], 'r-', 'LineWidth', 2);
plot([trackStartX, trackEndX], [trackStartY, trackEndY], 'r--', 'LineWidth', 1);

xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Contour Map with Trees (spacing: %dm)', treeSpacing));
axis equal tight;
legend('Terrain', 'Left Trees', 'Right Trees', 'Track Boundary', 'Location', 'best');

% Plot 2: 2D elevation map with trees
subplot(1, 2, 2);
imagesc([xMinTotal xMaxTotal], [yMaxTotal yMinTotal], mergedElevation);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
hold on;

% Plot tree positions
plot(leftTreesUTM(:, 1), leftTreesUTM(:, 2), 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(rightTreesUTM(:, 1), rightTreesUTM(:, 2), 'g^', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

% Draw track boundaries
plot([trackStartX + halfWidth*perpX, trackEndX + halfWidth*perpX], ...
     [trackStartY + halfWidth*perpY, trackEndY + halfWidth*perpY], 'r-', 'LineWidth', 2);
plot([trackStartX - halfWidth*perpX, trackEndX - halfWidth*perpX], ...
     [trackStartY - halfWidth*perpY, trackEndY - halfWidth*perpY], 'r-', 'LineWidth', 2);

xlabel('Easting (m)');
ylabel('Northing (m)');
title('Elevation Map with Tree Positions');
axis xy;
axis equal tight;

fprintf('\n--- Complete ---\n');
fprintf('KML file: %s\n', kmlFile);
fprintf('Total trees: %d (%d per side)\n', numTrees * 2, numTrees);
fprintf('Tree spacing: %d m\n', treeSpacing);

%% ========================================================================
%  FUNCTION: utm2wgs84
%  Converts UTM coordinates to WGS84 (latitude, longitude)
%  ========================================================================
function [lat, lon] = utm2wgs84(easting, northing, zone, hemisphere)
    % WGS84 ellipsoid parameters
    a = 6378137.0;              % Semi-major axis
    f = 1 / 298.257223563;      % Flattening
    e2 = 2*f - f^2;             % First eccentricity squared
    e_prime2 = e2 / (1 - e2);   % Second eccentricity squared
    
    % UTM parameters
    k0 = 0.9996;                % Scale factor
    
    % Remove false easting
    x = easting - 500000;
    
    % Remove false northing for southern hemisphere
    if upper(hemisphere) == 'S'
        y = northing - 10000000;
    else
        y = northing;
    end
    
    % Central meridian
    lon0 = (zone - 1) * 6 - 180 + 3;
    
    % Footpoint latitude
    M = y / k0;
    mu = M / (a * (1 - e2/4 - 3*e2^2/64 - 5*e2^3/256));
    
    e1 = (1 - sqrt(1 - e2)) / (1 + sqrt(1 - e2));
    
    phi1 = mu + (3*e1/2 - 27*e1^3/32) * sin(2*mu) ...
              + (21*e1^2/16 - 55*e1^4/32) * sin(4*mu) ...
              + (151*e1^3/96) * sin(6*mu) ...
              + (1097*e1^4/512) * sin(8*mu);
    
    % Calculate latitude and longitude
    N1 = a / sqrt(1 - e2 * sin(phi1)^2);
    T1 = tan(phi1)^2;
    C1 = e_prime2 * cos(phi1)^2;
    R1 = a * (1 - e2) / (1 - e2 * sin(phi1)^2)^1.5;
    D = x / (N1 * k0);
    
    lat = phi1 - (N1 * tan(phi1) / R1) * (D^2/2 ...
          - (5 + 3*T1 + 10*C1 - 4*C1^2 - 9*e_prime2) * D^4/24 ...
          + (61 + 90*T1 + 298*C1 + 45*T1^2 - 252*e_prime2 - 3*C1^2) * D^6/720);
    
    lon = lon0 + (D - (1 + 2*T1 + C1) * D^3/6 ...
          + (5 - 2*C1 + 28*T1 - 3*C1^2 + 8*e_prime2 + 24*T1^2) * D^5/120) / cos(phi1);
    
    % Convert to degrees
    lat = rad2deg(lat);
    lon = rad2deg(lon);
end
