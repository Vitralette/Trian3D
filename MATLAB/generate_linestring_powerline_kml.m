%% Generate LineString Powerline KML
% This script generates a KML file with powerlines placed orthogonally
% across the track as obstacles for pilots (fly over or under).
%
% Interactive mode: View the placement, then choose to REROLL or SAVE.
%
% Output: linestring_powerline.kml
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc; close all;

%% Configuration
editedFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Export');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Powerline placement parameters
numPowerlines = 5;              % Number of powerline obstacles
powerlineSeed = 1;              % Starting seed (will increment on reroll)
powerlineExtension = 20;       % How far powerlines extend beyond track (meters per side)
minDistanceBetween = 200;       % Minimum distance between powerlines (meters)
minDistanceFromEnds = 100;      % Minimum distance from track start/end (meters)
objectName = 'powerline';       % Object name recognized by Trian3D
placementMethod = 'linestring'; % KML geometry type

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

fprintf('\nPowerline parameters:\n');
fprintf('  Number of powerlines: %d\n', numPowerlines);
fprintf('  Powerline seed: %d\n', powerlineSeed);
fprintf('  Extension beyond track: %.0f m per side\n', powerlineExtension);
fprintf('  Min distance between: %.0f m\n', minDistanceBetween);

%% ========================================================================
%  INTERACTIVE LOOP - Generate, view, reroll or save
%  ========================================================================
userDone = false;

while ~userDone
    
% Close any existing figures from previous iteration
close all;

%% Build cumulative distance along track
numWaypoints = length(waypoints);
cumulativeDistance = zeros(numWaypoints, 1);

for i = 2:numWaypoints
    segDx = waypoints(i).x - waypoints(i-1).x;
    segDy = waypoints(i).y - waypoints(i-1).y;
    segLength = sqrt(segDx^2 + segDy^2);
    cumulativeDistance(i) = cumulativeDistance(i-1) + segLength;
end

%% Generate random powerline positions along track
fprintf('\nGenerating random powerline positions...\n');

rng(powerlineSeed);  % Set random seed for reproducibility

% Valid range for powerline placement
validStartDist = minDistanceFromEnds;
validEndDist = totalLength - minDistanceFromEnds;
validRange = validEndDist - validStartDist;

if validRange < (numPowerlines - 1) * minDistanceBetween
    warning('Track too short or too many powerlines requested. Reducing constraints.');
    minDistanceBetween = max(50, validRange / (numPowerlines + 1));
end

% Generate random positions with minimum spacing constraint
powerlineDistances = zeros(numPowerlines, 1);
maxAttempts = 1000;

for i = 1:numPowerlines
    attempts = 0;
    validPosition = false;
    
    while ~validPosition && attempts < maxAttempts
        % Generate random distance along track
        candidateDist = validStartDist + rand() * validRange;
        
        % Check distance from all existing powerlines
        if i == 1
            validPosition = true;
        else
            minDistToExisting = min(abs(powerlineDistances(1:i-1) - candidateDist));
            if minDistToExisting >= minDistanceBetween
                validPosition = true;
            end
        end
        attempts = attempts + 1;
    end
    
    if ~validPosition
        warning('Could not find valid position for powerline %d after %d attempts', i, maxAttempts);
    end
    
    powerlineDistances(i) = candidateDist;
end

% Sort powerlines by distance along track
powerlineDistances = sort(powerlineDistances);

fprintf('  Powerline positions along track (meters):\n');
for i = 1:numPowerlines
    fprintf('    Powerline %d: %.1f m (%.1f%%)\n', i, powerlineDistances(i), ...
        100 * powerlineDistances(i) / totalLength);
end

%% Calculate powerline endpoints (UTM coordinates)
fprintf('\nCalculating powerline geometry...\n');

% Structure to store powerline data
powerlines = struct('distance', {}, 'centerX', {}, 'centerY', {}, ...
                    'leftX', {}, 'leftY', {}, 'rightX', {}, 'rightY', {}, ...
                    'dirX', {}, 'dirY', {});

for p = 1:numPowerlines
    targetDist = powerlineDistances(p);
    
    % Find which segment contains this distance
    segIdx = find(cumulativeDistance <= targetDist, 1, 'last');
    if segIdx >= numWaypoints
        segIdx = numWaypoints - 1;
    end
    
    % Distance along this segment
    distIntoSeg = targetDist - cumulativeDistance(segIdx);
    
    % Segment parameters
    wp1 = waypoints(segIdx);
    wp2 = waypoints(segIdx + 1);
    segDx = wp2.x - wp1.x;
    segDy = wp2.y - wp1.y;
    segLength = sqrt(segDx^2 + segDy^2);
    
    if segLength < 0.1
        segLength = 0.1;  % Prevent division by zero
    end
    
    % Normalized direction along track
    dirX = segDx / segLength;
    dirY = segDy / segLength;
    
    % Perpendicular direction (left of travel = powerline direction)
    perpX = -dirY;
    perpY = dirX;
    
    % Center point on track
    t = distIntoSeg / segLength;
    centerX = wp1.x + t * segDx;
    centerY = wp1.y + t * segDy;
    
    % Powerline endpoints (extending beyond track corridor)
    extensionTotal = halfWidth + powerlineExtension;
    leftX = centerX + extensionTotal * perpX;
    leftY = centerY + extensionTotal * perpY;
    rightX = centerX - extensionTotal * perpX;
    rightY = centerY - extensionTotal * perpY;
    
    % Store powerline data
    powerlines(p).distance = targetDist;
    powerlines(p).centerX = centerX;
    powerlines(p).centerY = centerY;
    powerlines(p).leftX = leftX;
    powerlines(p).leftY = leftY;
    powerlines(p).rightX = rightX;
    powerlines(p).rightY = rightY;
    powerlines(p).dirX = dirX;
    powerlines(p).dirY = dirY;
    powerlines(p).perpX = perpX;
    powerlines(p).perpY = perpY;
end

fprintf('  Calculated %d powerline positions\n', numPowerlines);

%% Convert UTM to WGS84 (Latitude/Longitude)
fprintf('\nConverting UTM to WGS84...\n');

% UTM Zone 32N parameters
utmZone = 32;
hemisphere = 'N';

% Convert powerline endpoints
for p = 1:numPowerlines
    [powerlines(p).leftLat, powerlines(p).leftLon] = ...
        utm2wgs84(powerlines(p).leftX, powerlines(p).leftY, utmZone, hemisphere);
    [powerlines(p).rightLat, powerlines(p).rightLon] = ...
        utm2wgs84(powerlines(p).rightX, powerlines(p).rightY, utmZone, hemisphere);
    [powerlines(p).centerLat, powerlines(p).centerLon] = ...
        utm2wgs84(powerlines(p).centerX, powerlines(p).centerY, utmZone, hemisphere);
end

fprintf('  Sample powerline (first):\n');
fprintf('    Left:  %.6f, %.6f (lon, lat)\n', powerlines(1).leftLon, powerlines(1).leftLat);
fprintf('    Right: %.6f, %.6f (lon, lat)\n', powerlines(1).rightLon, powerlines(1).rightLat);

%% Generate KML file
kmlFilename = sprintf('%s_%s.kml', placementMethod, objectName);
kmlFile = fullfile(outputFolder, kmlFilename);
fprintf('\nWriting KML file: %s\n', kmlFile);

fid = fopen(kmlFile, 'w');

% KML header
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<kml xmlns="http://www.opengis.net/kml/2.2">\n');
fprintf(fid, '  <Document>\n');

% Write powerlines as LineStrings
for p = 1:numPowerlines
    fprintf(fid, '    <Placemark>\n');
    fprintf(fid, '      <name>%s</name>\n', objectName);
    fprintf(fid, '      <LineString>\n');
    fprintf(fid, '        <coordinates>\n');
    % Write left endpoint, then right endpoint (simple 2-point line)
    fprintf(fid, '          %.14f,%.14f,0\n', powerlines(p).leftLon, powerlines(p).leftLat);
    fprintf(fid, '          %.14f,%.14f,0\n', powerlines(p).rightLon, powerlines(p).rightLat);
    fprintf(fid, '        </coordinates>\n');
    fprintf(fid, '      </LineString>\n');
    fprintf(fid, '    </Placemark>\n');
end

% KML footer
fprintf(fid, '  </Document>\n');
fprintf(fid, '</kml>\n');

fclose(fid);
fprintf('KML file written with %d powerlines.\n', numPowerlines);

%% Visualize: Contour plot with powerline markers
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

%% Contour plot with powerlines
figure('Name', sprintf('Powerline Placement (Seed: %d)', powerlineSeed), 'Position', [200 200 900 700]);
[X, Y] = meshgrid(linspace(xMinTotal, xMaxTotal, nColsTotal), ...
                   linspace(yMaxTotal, yMinTotal, nRowsTotal));
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Elevation Contour Map - Powerline Placement\nTrack Seed: %d | Powerline Seed: %d | N=%d', ...
    geom.randomSeed, powerlineSeed, numPowerlines));
axis equal tight;
hold on;

% Plot track centerline
trackX = [waypoints.x];
trackY = [waypoints.y];
plot(trackX, trackY, 'r-', 'LineWidth', 2, 'DisplayName', 'Track');
plot(waypoints(1).x, waypoints(1).y, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(waypoints(end).x, waypoints(end).y, 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

% Plot corridor boundaries (faint)
leftBoundX = [];
leftBoundY = [];
rightBoundX = [];
rightBoundY = [];
for i = 1:(numWaypoints - 1)
    wp1 = waypoints(i);
    wp2 = waypoints(i + 1);
    segDx = wp2.x - wp1.x;
    segDy = wp2.y - wp1.y;
    segLen = sqrt(segDx^2 + segDy^2);
    if segLen > 0.1
        perpX = -segDy / segLen;
        perpY = segDx / segLen;
        leftBoundX = [leftBoundX, wp1.x + halfWidth * perpX];
        leftBoundY = [leftBoundY, wp1.y + halfWidth * perpY];
        rightBoundX = [rightBoundX, wp1.x - halfWidth * perpX];
        rightBoundY = [rightBoundY, wp1.y - halfWidth * perpY];
    end
end
plot(leftBoundX, leftBoundY, 'r--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(rightBoundX, rightBoundY, 'r--', 'LineWidth', 0.5, 'HandleVisibility', 'off');

%% Read back KML and plot powerline waypoints from file
fprintf('\nVerifying KML coordinates by reading back from file...\n');

% Read the KML file
kmlContent = fileread(kmlFile);

% Count powerline placemarks
numPowerlinesInKML = length(regexp(kmlContent, '<name>powerline</name>'));
fprintf('  Found %d powerline placemarks in KML file\n', numPowerlinesInKML);

% Extract LineString coordinates for each powerline
linePattern = '<LineString>\s*<coordinates>\s*(.*?)\s*</coordinates>';
[lineMatches, ~] = regexp(kmlContent, linePattern, 'tokens', 'match');

fprintf('  Found %d LineString elements\n', length(lineMatches));

% Plot each powerline's waypoints from KML
kmlMarkerColors = lines(length(lineMatches));
for i = 1:length(lineMatches)
    coordStr = strtrim(lineMatches{i}{1});
    
    % Parse coordinate string (lon,lat,alt separated by whitespace/newlines)
    coordLines = strsplit(coordStr);
    coordLines = coordLines(~cellfun('isempty', coordLines));  % Remove empty entries
    
    numPoints = length(coordLines);
    lons = zeros(numPoints, 1);
    lats = zeros(numPoints, 1);
    eastings = zeros(numPoints, 1);
    northings = zeros(numPoints, 1);
    
    for j = 1:numPoints
        parts = strsplit(coordLines{j}, ',');
        lons(j) = str2double(parts{1});
        lats(j) = str2double(parts{2});
        % Convert to UTM for plotting
        [eastings(j), northings(j)] = wgs842utm(lats(j), lons(j), utmZone);
    end
    
    % Plot thin line connecting waypoints
    plot(eastings, northings, '-', 'Color', kmlMarkerColors(i,:), 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    
    % Plot markers at each waypoint
    plot(eastings, northings, 'o', 'Color', kmlMarkerColors(i,:), ...
        'MarkerSize', 8, 'MarkerFaceColor', kmlMarkerColors(i,:), ...
        'DisplayName', sprintf('PL%d KML (%d pts)', i, numPoints));
    
    fprintf('  Powerline %d: %d waypoints\n', i, numPoints);
    for j = 1:numPoints
        fprintf('    Point %d: lon=%.6f, lat=%.6f -> E=%.1f, N=%.1f\n', ...
            j, lons(j), lats(j), eastings(j), northings(j));
    end
end

legend('Location', 'best');

%% User Decision: Reroll or Save
fprintf('\n=== Powerlines Generated (Seed: %d) ===\n', powerlineSeed);
fprintf('  [R] Reroll - generate new random placement\n');
fprintf('  [S] Save   - keep this placement and continue\n');
userChoice = input('Your choice: ', 's');

if strcmpi(userChoice, 'S')
    % Save and exit loop
    fprintf('Saving powerlines with seed %d...\n', powerlineSeed);
    userDone = true;
else
    % Reroll: increment seed, loop again
    powerlineSeed = powerlineSeed + 1;
    fprintf('\nRerolling with new seed %d...\n\n', powerlineSeed);
end

end  % End of while ~userDone loop

%% Summary
fprintf('\n--- Complete ---\n');
fprintf('Generated contour plot with powerline positions.\n');
fprintf('KML file: %s\n', kmlFile);
fprintf('Total powerlines: %d\n', numPowerlines);
fprintf('Powerline seed: %d\n', powerlineSeed);
fprintf('Powerline extension: %.0f m beyond corridor\n', powerlineExtension);

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
