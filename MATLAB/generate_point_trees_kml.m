%% Generate Point Trees KML
% This script generates a KML file with trees placed as points along the 
% left and right boundaries of the track.
%
% Output: point_trees.kml
%
% Author: Tim Jusko
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
editedFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Export');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Placement parameters
treeSpacing = 20;       % Distance between trees in meters
objectName = 'trees';   % Object name recognized by Trian3D (from KML bible)
placementMethod = 'point';  % KML geometry type

%% Load track geometry
geometryFile = fullfile(editedFolder, 'track_geometry.mat');
if ~exist(geometryFile, 'file')
    error('Track geometry file not found. Run generate_groundtrack.m first.');
end
fprintf('Loading track geometry from: %s\n', geometryFile);
geomData = load(geometryFile);
geom = geomData.trackGeometry;

% Extract geometry parameters
waypoints = geom.waypoints;
corridorWidth = geom.corridorWidth;
halfWidth = corridorWidth / 2;
totalLength = geom.totalLength;
xMinTotal = geom.xMinTerrain;
xMaxTotal = geom.xMaxTerrain;
yMinTotal = geom.yMinTerrain;
yMaxTotal = geom.yMaxTerrain;

fprintf('\nTrack parameters (loaded from file):\n');
fprintf('  Seed: %d\n', geom.randomSeed);
fprintf('  Waypoints: %d\n', length(waypoints));
fprintf('  Total length: %.0f m\n', totalLength);
fprintf('  Corridor width: %.0f m\n', corridorWidth);

%% Build dense track path and calculate tree positions
fprintf('\nCalculating tree positions along track boundaries...\n');

% Arrays to store tree positions (UTM coordinates)
leftTreesUTM = [];
rightTreesUTM = [];

% Walk along track and place trees at regular intervals
distanceTraveled = 0;
nextTreeDistance = 0;
numWaypoints = length(waypoints);

for i = 1:(numWaypoints - 1)
    wp1 = waypoints(i);
    wp2 = waypoints(i + 1);
    
    % Segment vector
    segDx = wp2.x - wp1.x;
    segDy = wp2.y - wp1.y;
    segLength = sqrt(segDx^2 + segDy^2);
    
    if segLength < 0.1
        continue;  % Skip zero-length segments
    end
    
    % Normalized direction
    dirX = segDx / segLength;
    dirY = segDy / segLength;
    
    % Perpendicular vector (left of travel direction)
    perpX = -dirY;
    perpY = dirX;
    
    % Place trees along this segment
    segStart = 0;
    while distanceTraveled + segStart <= distanceTraveled + segLength
        % Check if we need to place a tree
        if distanceTraveled + segStart >= nextTreeDistance
            % Calculate position along segment
            t = segStart / segLength;
            centerX = wp1.x + t * segDx;
            centerY = wp1.y + t * segDy;
            
            % Left boundary tree
            leftX = centerX + halfWidth * perpX;
            leftY = centerY + halfWidth * perpY;
            leftTreesUTM = [leftTreesUTM; leftX, leftY];
            
            % Right boundary tree
            rightX = centerX - halfWidth * perpX;
            rightY = centerY - halfWidth * perpY;
            rightTreesUTM = [rightTreesUTM; rightX, rightY];
            
            nextTreeDistance = nextTreeDistance + treeSpacing;
        end
        segStart = segStart + 1;  % Move 1 meter at a time
    end
    
    distanceTraveled = distanceTraveled + segLength;
end

numTrees = size(leftTreesUTM, 1);
fprintf('  Generated %d trees per side (%d total)\n', numTrees, numTrees * 2);

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
kmlFilename = sprintf('%s_%s.kml', placementMethod, objectName);
kmlFile = fullfile(outputFolder, kmlFilename);
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
    fprintf(fid, '      <name>%s</name>\n', objectName);
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
    fprintf(fid, '      <name>%s</name>\n', objectName);
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

% Load elevation data for visualization
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};
numTiles = length(tifBaseNames);
geoInfo = cell(numTiles, 1);
elevationData = cell(numTiles, 1);

for i = 1:numTiles
    baseName = tifBaseNames{i};
    
    % Load from edited folder
    geoInfoFile = fullfile(editedFolder, [baseName '_geotiffinfo_edited.mat']);
    geoInfoData = load(geoInfoFile);
    geoInfoFields = fieldnames(geoInfoData);
    geoInfo{i} = geoInfoData.(geoInfoFields{1});
    
    rasterFile = fullfile(editedFolder, [baseName '_readgeoraster_edited.mat']);
    rasterData = load(rasterFile);
    elevationData{i} = rasterData.A;
end

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

%% Contour plot
figure('Name', sprintf('Elevation Contours - Trees (Seed: %d)', geom.randomSeed), 'Position', [200 200 800 600]);
[X, Y] = meshgrid(linspace(xMinTotal, xMaxTotal, nColsTotal), ...
                   linspace(yMaxTotal, yMinTotal, nRowsTotal));
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Elevation Contour Map - Tree Placement (Seed: %d)', geom.randomSeed));
axis equal tight;
hold on;

% Plot track centerline
trackX = [waypoints.x];
trackY = [waypoints.y];
plot(trackX, trackY, 'r-', 'LineWidth', 2, 'DisplayName', 'Track');
plot(waypoints(1).x, waypoints(1).y, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(waypoints(end).x, waypoints(end).y, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

%% Read back KML and plot tree positions for verification
fprintf('\nVerifying KML coordinates by reading back from file...\n');

% Read the KML file
kmlContent = fileread(kmlFile);

% Extract all coordinate entries using regex
coordPattern = '<coordinates>\s*([-\d.]+),([-\d.]+),';
[tokens, ~] = regexp(kmlContent, coordPattern, 'tokens', 'match');

% Convert WGS84 coordinates back to UTM and plot
numTreesFromKML = length(tokens);
fprintf('  Found %d tree coordinates in KML file\n', numTreesFromKML);

treeEastings = zeros(numTreesFromKML, 1);
treeNorthings = zeros(numTreesFromKML, 1);

for i = 1:numTreesFromKML
    lon = str2double(tokens{i}{1});
    lat = str2double(tokens{i}{2});
    [treeEastings(i), treeNorthings(i)] = wgs842utm(lat, lon, utmZone);
end

% Plot tree markers from KML
plot(treeEastings, treeNorthings, 'g^', 'MarkerSize', 6, 'MarkerFaceColor', 'g', ...
    'DisplayName', sprintf('Trees from KML (%d)', numTreesFromKML));
legend('Location', 'best');

fprintf('  Sample tree from KML (first): E=%.2f, N=%.2f\n', treeEastings(1), treeNorthings(1));
fprintf('  Sample tree from KML (last):  E=%.2f, N=%.2f\n', treeEastings(end), treeNorthings(end));

fprintf('\n--- Complete ---\n');
fprintf('Generated contour plot of elevation data with tree positions.\n');
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
    
    % Central meridian (in radians)
    lon0_deg = (zone - 1) * 6 - 180 + 3;
    lon0 = deg2rad(lon0_deg);
    
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

%% ========================================================================
%  FUNCTION: wgs842utm
%  Converts WGS84 (latitude, longitude) to UTM coordinates
%  ========================================================================
function [easting, northing] = wgs842utm(lat, lon, zone)
    % WGS84 ellipsoid parameters
    a = 6378137.0;              % Semi-major axis
    f = 1 / 298.257223563;      % Flattening
    e2 = 2*f - f^2;             % First eccentricity squared
    e_prime2 = e2 / (1 - e2);   % Second eccentricity squared
    
    % UTM parameters
    k0 = 0.9996;                % Scale factor
    
    % Convert to radians
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);
    
    % Central meridian
    lon0_deg = (zone - 1) * 6 - 180 + 3;
    lon0 = deg2rad(lon0_deg);
    
    % Calculate intermediate values
    N = a / sqrt(1 - e2 * sin(lat_rad)^2);
    T = tan(lat_rad)^2;
    C = e_prime2 * cos(lat_rad)^2;
    A = (lon_rad - lon0) * cos(lat_rad);
    
    % Meridional arc
    M = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * lat_rad ...
           - (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*lat_rad) ...
           + (15*e2^2/256 + 45*e2^3/1024) * sin(4*lat_rad) ...
           - (35*e2^3/3072) * sin(6*lat_rad));
    
    % Calculate easting and northing
    easting = k0 * N * (A + (1 - T + C) * A^3/6 ...
              + (5 - 18*T + T^2 + 72*C - 58*e_prime2) * A^5/120) + 500000;
    
    northing = k0 * (M + N * tan(lat_rad) * (A^2/2 ...
               + (5 - T + 9*C + 4*C^2) * A^4/24 ...
               + (61 - 58*T + T^2 + 600*C - 330*e_prime2) * A^6/720));
    
    % Add false northing for southern hemisphere
    if lat < 0
        northing = northing + 10000000;
    end
end
