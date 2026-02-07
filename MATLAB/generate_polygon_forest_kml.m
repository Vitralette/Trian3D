%% Generate Polygon Forest KML
% This script generates a KML file with a forest polygon "tube" surrounding
% the race track. The forest extends from the track corridor edge outward.
%
% Interactive mode: View the placement, then choose to REROLL or SAVE.
%
% Output: polygon_forest.kml
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc; close all;

%% Load Project Configuration
config = load_project_config();
editedFolder = config.editedFolder;
outputFolder = config.exportFolder;

fprintf('Project: %s\n', config.projectName);

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Forest tube parameters
forestInnerOffset = 0;          % Distance from track edge to forest inner boundary (meters)
                                 % 0 = forest starts right at corridor edge
forestWidth = 150;               % Width of the forest tube (meters)
pointSpacing = 10;               % Spacing between polygon vertices along track (meters)
objectName = 'forest';           % Object name recognized by Trian3D
placementMethod = 'polygon';     % KML geometry type

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
fprintf('  Track Seed: %d\n', geom.randomSeed);
fprintf('  Waypoints: %d\n', length(waypoints));
fprintf('  Total length: %.0f m\n', totalLength);
fprintf('  Corridor width: %.0f m\n', corridorWidth);

fprintf('\nForest tube parameters:\n');
fprintf('  Inner offset from track edge: %.0f m\n', forestInnerOffset);
fprintf('  Forest width: %.0f m\n', forestWidth);
fprintf('  Point spacing: %.0f m\n', pointSpacing);

%% Build dense track path for polygon generation
fprintf('\nBuilding dense track path...\n');

numWaypoints = length(waypoints);

% Calculate cumulative distances
cumulativeDistance = zeros(numWaypoints, 1);
for i = 2:numWaypoints
    segDx = waypoints(i).x - waypoints(i-1).x;
    segDy = waypoints(i).y - waypoints(i-1).y;
    segLength = sqrt(segDx^2 + segDy^2);
    cumulativeDistance(i) = cumulativeDistance(i-1) + segLength;
end

% Generate evenly spaced points along track
numPoints = ceil(totalLength / pointSpacing) + 1;
trackDistances = linspace(0, totalLength, numPoints);

% Arrays for track centerline and direction vectors
centerX = zeros(numPoints, 1);
centerY = zeros(numPoints, 1);
dirX = zeros(numPoints, 1);  % Direction along track (unit vector)
dirY = zeros(numPoints, 1);

for p = 1:numPoints
    targetDist = trackDistances(p);
    
    % Find which segment contains this distance
    segIdx = find(cumulativeDistance <= targetDist, 1, 'last');
    if segIdx >= numWaypoints
        segIdx = numWaypoints - 1;
    end
    
    % Distance into this segment
    distIntoSeg = targetDist - cumulativeDistance(segIdx);
    
    % Segment parameters
    wp1 = waypoints(segIdx);
    wp2 = waypoints(segIdx + 1);
    segDx = wp2.x - wp1.x;
    segDy = wp2.y - wp1.y;
    segLength = sqrt(segDx^2 + segDy^2);
    
    if segLength < 0.1
        segLength = 0.1;
    end
    
    % Normalized direction along track
    dirX(p) = segDx / segLength;
    dirY(p) = segDy / segLength;
    
    % Center point
    t = distIntoSeg / segLength;
    t = max(0, min(1, t));  % Clamp to [0, 1]
    centerX(p) = wp1.x + t * segDx;
    centerY(p) = wp1.y + t * segDy;
end

% Smooth the direction vectors to avoid abrupt changes at segment boundaries
smoothingWindow = max(3, round(50 / pointSpacing));  % ~50m smoothing window
if smoothingWindow > 1
    % Moving average smoothing (preserving endpoints)
    dirXSmooth = dirX;
    dirYSmooth = dirY;
    halfWin = floor(smoothingWindow / 2);
    
    for p = (halfWin + 1):(numPoints - halfWin)
        dirXSmooth(p) = mean(dirX(p-halfWin:p+halfWin));
        dirYSmooth(p) = mean(dirY(p-halfWin:p+halfWin));
    end
    
    % Re-normalize smoothed directions
    for p = 1:numPoints
        len = sqrt(dirXSmooth(p)^2 + dirYSmooth(p)^2);
        if len > 0.01
            dirXSmooth(p) = dirXSmooth(p) / len;
            dirYSmooth(p) = dirYSmooth(p) / len;
        end
    end
    
    dirX = dirXSmooth;
    dirY = dirYSmooth;
end

% Calculate perpendicular directions (left of travel)
perpX = -dirY;
perpY = dirX;

fprintf('  Generated %d points along track\n', numPoints);

%% Calculate forest polygon boundaries
% Inner boundary = track edge + offset
% Outer boundary = inner boundary + forest width

innerOffset = halfWidth + forestInnerOffset;
outerOffset = innerOffset + forestWidth;

% Left side boundaries (positive perpendicular)
leftInnerX = centerX + innerOffset * perpX;
leftInnerY = centerY + innerOffset * perpY;
leftOuterX = centerX + outerOffset * perpX;
leftOuterY = centerY + outerOffset * perpY;

% Right side boundaries (negative perpendicular)
rightInnerX = centerX - innerOffset * perpX;
rightInnerY = centerY - innerOffset * perpY;
rightOuterX = centerX - outerOffset * perpX;
rightOuterY = centerY - outerOffset * perpY;

%% Build forest polygon boundaries as simple closed loops
% A "tube" polygon goes: outer edge forward, then inner edge backward, then closes
fprintf('  Building polygon boundaries...\n');

% LEFT SIDE POLYGON: outer edge forward -> inner edge backward -> close
% This creates a closed loop: start at outer[1], go to outer[end], 
% then to inner[end], back to inner[1], close to outer[1]
leftPolyX = [leftOuterX; flipud(leftInnerX); leftOuterX(1)];
leftPolyY = [leftOuterY; flipud(leftInnerY); leftOuterY(1)];

% RIGHT SIDE POLYGON: inner edge forward -> outer edge backward -> close
rightPolyX = [rightInnerX; flipud(rightOuterX); rightInnerX(1)];
rightPolyY = [rightInnerY; flipud(rightOuterY); rightInnerY(1)];

fprintf('  Left polygon: %d vertices\n', length(leftPolyX));
fprintf('  Right polygon: %d vertices\n', length(rightPolyX));

% Collect polygons for KML export (simple - no polyshape processing)
allPolygons = {};
allPolygons{1} = struct('x', leftPolyX, 'y', leftPolyY, 'side', 'left');
allPolygons{2} = struct('x', rightPolyX, 'y', rightPolyY, 'side', 'right');

fprintf('  Total polygon regions to export: %d\n', length(allPolygons));

%% Convert UTM to WGS84
fprintf('\nConverting UTM to WGS84...\n');

utmZone = 32;
hemisphere = 'N';

% Convert all polygon regions
for p = 1:length(allPolygons)
    numVerts = length(allPolygons{p}.x);
    allPolygons{p}.lon = zeros(numVerts, 1);
    allPolygons{p}.lat = zeros(numVerts, 1);
    for i = 1:numVerts
        [allPolygons{p}.lat(i), allPolygons{p}.lon(i)] = ...
            utm2wgs84(allPolygons{p}.x(i), allPolygons{p}.y(i), utmZone, hemisphere);
    end
end

fprintf('  Conversion complete\n');

%% Generate KML file
kmlFilename = sprintf('%s_%s.kml', placementMethod, objectName);
kmlFile = fullfile(outputFolder, kmlFilename);
fprintf('\nWriting KML file: %s\n', kmlFile);

fid = fopen(kmlFile, 'w');

% KML header
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid, '  <Document>\n');

% Write all polygon regions
for p = 1:length(allPolygons)
    fprintf(fid, '    <Placemark>\n');
    fprintf(fid, '      <name>%s</name>\n', objectName);
    fprintf(fid, '      <Polygon>\n');
    fprintf(fid, '        <outerBoundaryIs>\n');
    fprintf(fid, '          <LinearRing>\n');
    fprintf(fid, '            <coordinates>\n');
    for i = 1:length(allPolygons{p}.lon)
        fprintf(fid, '              %.14f,%.14f,0\n', allPolygons{p}.lon(i), allPolygons{p}.lat(i));
    end
    fprintf(fid, '            </coordinates>\n');
    fprintf(fid, '          </LinearRing>\n');
    fprintf(fid, '        </outerBoundaryIs>\n');
    fprintf(fid, '      </Polygon>\n');
    fprintf(fid, '    </Placemark>\n');
end

% KML footer
fprintf(fid, '  </Document>\n');
fprintf(fid, '</kml>\n');

fclose(fid);
fprintf('KML file written with %d forest polygon(s).\n', length(allPolygons));

%% Read back KML and save reconstructed data for visualization
fprintf('\nReading back KML to verify and save visualization data...\n');

% Read the KML file we just wrote
kmlContent = fileread(kmlFile);

% Extract each Placemark's polygon coordinates
% More specific pattern: match coordinates within LinearRing
polyPattern = '<coordinates>\s*([\d\s.,\-\neE]+?)\s*</coordinates>';
[polyMatches, ~] = regexp(kmlContent, polyPattern, 'tokens', 'match');

numPolysFromKML = length(polyMatches);
fprintf('  Found %d polygon(s) in KML file\n', numPolysFromKML);

% Convert WGS84 coordinates back to UTM (this is what Trian3D will interpret)
kmlPolygons = cell(numPolysFromKML, 1);
for p = 1:numPolysFromKML
    coordStr = strtrim(polyMatches{p}{1});
    
    % Split by newlines or whitespace to get individual coordinate triplets
    coordLines = regexp(coordStr, '\s+', 'split');
    coordLines = coordLines(~cellfun('isempty', coordLines));
    
    numVerts = length(coordLines);
    eastings = zeros(numVerts, 1);
    northings = zeros(numVerts, 1);
    
    fprintf('  Polygon %d: %d vertices\n', p, numVerts);
    
    for j = 1:numVerts
        parts = strsplit(coordLines{j}, ',');
        if length(parts) >= 2
            lon = str2double(parts{1});
            lat = str2double(parts{2});
            [eastings(j), northings(j)] = wgs842utm(lat, lon, utmZone);
        end
    end
    
    kmlPolygons{p}.x = eastings;
    kmlPolygons{p}.y = northings;
end

% Save reconstructed-from-KML data as .mat for visualization
matFilename = 'polygon_forest_data.mat';
matFile = fullfile(editedFolder, matFilename);

% Create structure with polygon data reconstructed FROM KML
forestData.polygons = kmlPolygons;  % Each has .x, .y (UTM) - reconstructed from KML
forestData.numPolygons = numPolysFromKML;
forestData.forestInnerOffset = forestInnerOffset;
forestData.forestWidth = forestWidth;
forestData.trackSeed = geom.randomSeed;

save(matFile, 'forestData');
fprintf('Reconstructed polygon data saved to: %s\n', matFile);

%% Visualize: Contour plot with forest polygons
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

%% Contour plot with forest polygons
figure('Name', sprintf('Forest Tube (Track Seed: %d)', geom.randomSeed), 'Position', [200 200 900 700]);
[X, Y] = meshgrid(linspace(xMinTotal, xMaxTotal, nColsTotal), ...
                   linspace(yMaxTotal, yMinTotal, nRowsTotal));
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Elevation Contour Map - Forest Tube\nTrack Seed: %d | Inner Offset: %.0fm | Forest Width: %.0fm', ...
    geom.randomSeed, forestInnerOffset, forestWidth));
axis equal tight;
hold on;

% Plot track centerline
trackX = [waypoints.x];
trackY = [waypoints.y];
plot(trackX, trackY, 'r-', 'LineWidth', 2, 'DisplayName', 'Track');
plot(waypoints(1).x, waypoints(1).y, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(waypoints(end).x, waypoints(end).y, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

% Plot corridor boundaries
plot(leftInnerX, leftInnerY, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(rightInnerX, rightInnerY, 'r--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Plot forest polygons RECONSTRUCTED FROM KML (exactly what Trian3D will see)
forestColor = [0.2 0.6 0.2];  % Dark green
for p = 1:numPolysFromKML
    if p == 1
        patch(kmlPolygons{p}.x, kmlPolygons{p}.y, forestColor, ...
            'FaceAlpha', 0.4, 'EdgeColor', forestColor, 'LineWidth', 1, ...
            'DisplayName', sprintf('Forest from KML (%d polys)', numPolysFromKML));
    else
        patch(kmlPolygons{p}.x, kmlPolygons{p}.y, forestColor, ...
            'FaceAlpha', 0.4, 'EdgeColor', forestColor, 'LineWidth', 1, ...
            'HandleVisibility', 'off');
    end
end

legend('Location', 'best');

%% Summary
fprintf('\n--- Complete ---\n');
fprintf('Generated contour plot with forest tube.\n');
fprintf('KML file: %s\n', kmlFile);
fprintf('Forest inner offset: %.0f m from track edge\n', forestInnerOffset);
fprintf('Forest width: %.0f m\n', forestWidth);
fprintf('Total polygon regions: %d\n', length(allPolygons));

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
    a = 6378137.0;
    f = 1 / 298.257223563;
    e2 = 2*f - f^2;
    e_prime2 = e2 / (1 - e2);
    
    k0 = 0.9996;
    
    lat_rad = deg2rad(lat);
    lon_rad = deg2rad(lon);
    
    lon0_deg = (zone - 1) * 6 - 180 + 3;
    lon0 = deg2rad(lon0_deg);
    
    N = a / sqrt(1 - e2 * sin(lat_rad)^2);
    T = tan(lat_rad)^2;
    C = e_prime2 * cos(lat_rad)^2;
    A = (lon_rad - lon0) * cos(lat_rad);
    
    M = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * lat_rad ...
           - (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*lat_rad) ...
           + (15*e2^2/256 + 45*e2^3/1024) * sin(4*lat_rad) ...
           - (35*e2^3/3072) * sin(6*lat_rad));
    
    easting = k0 * N * (A + (1 - T + C) * A^3/6 ...
              + (5 - 18*T + T^2 + 72*C - 58*e_prime2) * A^5/120) + 500000;
    
    northing = k0 * (M + N * tan(lat_rad) * (A^2/2 ...
               + (5 - T + 9*C + 4*C^2) * A^4/24 ...
               + (61 - 58*T + T^2 + 600*C - 330*e_prime2) * A^6/720));
    
    if lat < 0
        northing = northing + 10000000;
    end
end