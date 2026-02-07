%% Generate Ground Track from Course Structure
% This script takes a course structure (from define_course_structure.m) and
% generates a random EVENT SEQUENCE from the event pool, then builds the track.
%
% Interactive mode: View the track, then choose to REROLL or SAVE.
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc; close all;

%% Configuration
dataFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Raw');
editedFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');

%% Initial seed
randomSeed = 1;  % Starting seed (will increment on reroll)

%% Load Course Structure (once, outside loop)
structureFile = fullfile(editedFolder, 'event_structure.mat');
if ~exist(structureFile, 'file')
    error('Course structure file not found. Run define_course_structure.m first.');
end
load(structureFile, 'eventStructure');

eventPool = eventStructure.eventPool;
params = eventStructure.params;
corridorWidth = eventStructure.corridorWidth;
transition = eventStructure.transition;
forceStartStraight = eventStructure.forceStartStraight;
forceEndStraight = eventStructure.forceEndStraight;

%% Load terrain data (once, outside loop)
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};
numTiles = length(tifBaseNames);
geoInfo = cell(numTiles, 1);
elevationData = cell(numTiles, 1);

for i = 1:numTiles
    baseName = tifBaseNames{i};
    geoInfoFile = fullfile(dataFolder, [baseName '_geotiffinfo.mat']);
    geoInfoData = load(geoInfoFile);
    geoInfoFields = fieldnames(geoInfoData);
    geoInfo{i} = geoInfoData.(geoInfoFields{1});
    
    rasterFile = fullfile(dataFolder, [baseName '_readgeoraster.mat']);
    rasterData = load(rasterFile);
    elevationData{i} = rasterData.A;
end

% Calculate terrain bounds
allXLimits = zeros(numTiles, 2);
allYLimits = zeros(numTiles, 2);
for i = 1:numTiles
    bbox = geoInfo{i}.BoundingBox;
    allXLimits(i, :) = [bbox(1,1), bbox(2,1)];
    allYLimits(i, :) = [bbox(1,2), bbox(2,2)];
end
xMinTerrain = min(allXLimits(:, 1));
xMaxTerrain = max(allXLimits(:, 2));
yMinTerrain = min(allYLimits(:, 1));
yMaxTerrain = max(allYLimits(:, 2));

terrainWidth = xMaxTerrain - xMinTerrain;
terrainHeight = yMaxTerrain - yMinTerrain;

fprintf('Terrain extent: %.0f x %.0f m\n\n', terrainWidth, terrainHeight);

%% ========================================================================
%  INTERACTIVE LOOP - Generate, view, reroll or save
%  ========================================================================
userDone = false;
figHandle = [];

while ~userDone
    clc;
    rng(randomSeed);
    fprintf('=== Ground Track Generator ===\n');
    fprintf('Random seed: %d\n', randomSeed);
    fprintf('=========================================\n\n');

    %% Generate Random Event Sequence from Pool
    fprintf('Generating random event sequence from pool...\n');
    
    % Build list of all events from pool
    allEvents = {};
    for i = 1:size(eventPool, 1)
        eventType = eventPool{i, 1};
        eventCount = eventPool{i, 2};
        for j = 1:eventCount
            allEvents{end+1} = eventType;
        end
    end
    
    % Shuffle the events
    shuffledIdx = randperm(length(allEvents));
    eventSequence = allEvents(shuffledIdx);
    
    % Handle forced start/end straights
    if forceStartStraight
        % Find a straight and move it to the front
        straightIdx = find(strcmp(eventSequence, 'straight'), 1);
        if ~isempty(straightIdx) && straightIdx ~= 1
            eventSequence([1, straightIdx]) = eventSequence([straightIdx, 1]);
        end
    end

    if forceEndStraight
        % Find a straight (not the first one) and move it to the end
        straightIndices = find(strcmp(eventSequence, 'straight'));
        if length(straightIndices) > 1
            lastStraightIdx = straightIndices(end);
            if lastStraightIdx ~= length(eventSequence)
                eventSequence([lastStraightIdx, end]) = eventSequence([end, lastStraightIdx]);
            end
        elseif length(straightIndices) == 1 && ~forceStartStraight
            % Only one straight, move it to end
            straightIdx = straightIndices(1);
            if straightIdx ~= length(eventSequence)
                eventSequence([straightIdx, end]) = eventSequence([end, straightIdx]);
            end
        end
    end

    fprintf('Generated sequence: %s\n', strjoin(eventSequence, ' -> '));
    fprintf('Total events: %d\n\n', length(eventSequence));

    %% Generate Track by Processing Events
    % SMART START: Position near one end, heading towards the other end
    % This maximizes usable track length on a rectangular terrain

    margin = 100;  % meters from terrain edge (reduced for tight terrain)

    % For a wide terrain (width > height), start on left/right edge
    % For a tall terrain (height > width), start on top/bottom edge
    if terrainWidth > terrainHeight
        % Wide terrain: start on left side, head right (east)
        startX = xMinTerrain + margin + rand() * 100;
        startY = yMinTerrain + terrainHeight/3 + rand() * terrainHeight/3;  % Middle third
        heading = (rand() - 0.5) * 30;  % Generally east, +/- 15 degrees
    else
        % Tall terrain: start on bottom, head north
        startX = xMinTerrain + terrainWidth/3 + rand() * terrainWidth/3;
        startY = yMinTerrain + margin + rand() * 100;
        heading = 90 + (rand() - 0.5) * 30;  % Generally north
    end

    % Normalize heading to 0-360
    heading = mod(heading, 360);

    % Initialize track as waypoints
    waypoints = [];
    waypoints(1).x = startX;
waypoints(1).y = startY;
waypoints(1).z = 0;  % Relative elevation (will be carved into terrain)
waypoints(1).type = 'start';

currentX = startX;
currentY = startY;
currentZ = 0;
currentHeading = heading;

% Track previous turn direction for 'opposite' handling
lastTurnDir = 0;  % 0 = no previous turn, 1 = right, -1 = left

fprintf('Generating track...\n');
fprintf('Start: (%.0f, %.0f), Heading: %.0f deg\n\n', currentX, currentY, currentHeading);

generatedEvents = struct('type', {}, 'actualParams', {});

% Define safe boundaries (with margin)
safeXMin = xMinTerrain + margin;
safeXMax = xMaxTerrain - margin;
safeYMin = yMinTerrain + margin;
safeYMax = yMaxTerrain - margin;

for i = 1:length(eventSequence)
    eventType = eventSequence{i};
    fprintf('Event %d: %s\n', i, upper(eventType));
    
    switch eventType
        case 'straight'
            % Fixed length
            len = params.straight.length;
            
            % Calculate end point
            endX = currentX + len * cosd(currentHeading);
            endY = currentY + len * sind(currentHeading);
            
            % Add waypoint
            wp.x = endX;
            wp.y = endY;
            wp.z = currentZ;
            wp.type = 'straight';
            waypoints(end+1) = wp;
            
            % Update state
            currentX = endX;
            currentY = endY;
            
            % Store actual parameters
            generatedEvents(i).type = 'straight';
            generatedEvents(i).actualParams.length = len;
            
            fprintf('  Length: %.0f m\n', len);
            
        case {'turn', 'turn_l', 'turn_r'}
            % Get parameters based on turn type
            if strcmp(eventType, 'turn_l')
                radius = params.turn_l.radius;
                angle = params.turn_l.angle;
                turnDir = -1;
                dirStr = 'left';
            elseif strcmp(eventType, 'turn_r')
                radius = params.turn_r.radius;
                angle = params.turn_r.angle;
                turnDir = 1;
                dirStr = 'right';
            else
                % Generic turn - alternate direction
                radius = params.turn.radius;
                angle = params.turn.angle;
                
                if lastTurnDir == 0
                    % First turn: random direction
                    if rand() > 0.5
                        turnDir = 1;
                        dirStr = 'right';
                    else
                        turnDir = -1;
                        dirStr = 'left';
                    end
                else
                    % Alternate from last turn
                    turnDir = -lastTurnDir;
                    if turnDir == 1
                        dirStr = 'right';
                    else
                        dirStr = 'left';
                    end
                end
            end
            
            % Remember this turn direction
            lastTurnDir = turnDir;
            
            % ============================================================
            % TURN GEOMETRY (correct derivation from first principles)
            % ============================================================
            % Aircraft at position P, heading H (degrees, 0=East, 90=North)
            % Turn through angle A with radius R
            %
            % RIGHT TURN (turnDir = +1):
            %   - Turn center is perpendicular RIGHT of heading
            %   - Center direction = H - 90°
            %   - Aircraft rotates CLOCKWISE around center
            %   - Exit heading = H - A (heading decreases)
            %
            % LEFT TURN (turnDir = -1):
            %   - Turn center is perpendicular LEFT of heading
            %   - Center direction = H + 90°
            %   - Aircraft rotates COUNTER-CLOCKWISE around center
            %   - Exit heading = H + A (heading increases)
            %
            % Unified: center direction = H - turnDir * 90°
            %          exit heading = H - turnDir * A
            % ============================================================
            
            % Find turn center
            centerDir = currentHeading - turnDir * 90;
            turnCenterX = currentX + radius * cosd(centerDir);
            turnCenterY = currentY + radius * sind(centerDir);
            
            % Starting angle (from center's perspective, pointing to aircraft)
            startAngle = currentHeading + turnDir * 90;  % Opposite of centerDir
            
            % Generate arc waypoints
            numArcPoints = max(10, round(angle / 5));
            
            for j = 1:numArcPoints
                progress = j / numArcPoints;
                
                % Angular position along arc
                % Right turn: angle decreases (clockwise)
                % Left turn: angle increases (counter-clockwise)
                arcAngle = startAngle - turnDir * angle * progress;
                
                wpX = turnCenterX + radius * cosd(arcAngle);
                wpY = turnCenterY + radius * sind(arcAngle);
                
                wp.x = wpX;
                wp.y = wpY;
                wp.z = currentZ;
                wp.type = 'turn';
                waypoints(end+1) = wp;
            end
            
            % Update position
            currentX = waypoints(end).x;
            currentY = waypoints(end).y;
            
            % Exit heading (analytically correct)
            currentHeading = currentHeading - turnDir * angle;
            
            % Store actual parameters
            generatedEvents(i).type = 'turn';
            generatedEvents(i).actualParams.radius = radius;
            generatedEvents(i).actualParams.angle = angle;
            generatedEvents(i).actualParams.direction = dirStr;
            
            fprintf('  Radius: %.0f m, Angle: %.0f deg, Direction: %s\n', radius, angle, dirStr);
            
        case 'climb'
            % Fixed parameters
            len = params.climb.length;
            gradient = params.climb.gradient;
            
            % Calculate elevation change
            elevChange = len * (gradient / 100);
            
            % Calculate end point
            endX = currentX + len * cosd(currentHeading);
            endY = currentY + len * sind(currentHeading);
            endZ = currentZ + elevChange;
            
            % Add waypoint
            wp.x = endX;
            wp.y = endY;
            wp.z = endZ;
            wp.type = 'climb';
            waypoints(end+1) = wp;
            
            % Update state
            currentX = endX;
            currentY = endY;
            currentZ = endZ;
            
            % Store actual parameters
            generatedEvents(i).type = 'climb';
            generatedEvents(i).actualParams.length = len;
            generatedEvents(i).actualParams.gradient = gradient;
            generatedEvents(i).actualParams.elevChange = elevChange;
            
            fprintf('  Length: %.0f m, Gradient: %.1f%%, Elev change: +%.1f m\n', len, gradient, elevChange);
            
        case 'descent'
            % Fixed parameters
            len = params.descent.length;
            gradient = params.descent.gradient;
            
            % Calculate elevation change (negative for descent)
            elevChange = -len * (gradient / 100);
            
            % Calculate end point
            endX = currentX + len * cosd(currentHeading);
            endY = currentY + len * sind(currentHeading);
            endZ = currentZ + elevChange;
            
            % Add waypoint
            wp.x = endX;
            wp.y = endY;
            wp.z = endZ;
            wp.type = 'descent';
            waypoints(end+1) = wp;
            
            % Update state
            currentX = endX;
            currentY = endY;
            currentZ = endZ;
            
            % Store actual parameters
            generatedEvents(i).type = 'descent';
            generatedEvents(i).actualParams.length = len;
            generatedEvents(i).actualParams.gradient = gradient;
            generatedEvents(i).actualParams.elevChange = elevChange;
            
            fprintf('  Length: %.0f m, Gradient: %.1f%%, Elev change: %.1f m\n', len, gradient, elevChange);
            
        case {'climb_turn', 'climb_turn_l', 'climb_turn_r'}
            % Get parameters based on turn type
            if strcmp(eventType, 'climb_turn_l')
                radius = params.climb_turn_l.radius;
                angle = params.climb_turn_l.angle;
                gradient = params.climb_turn_l.gradient;
                turnDir = -1;
                dirStr = 'left';
            elseif strcmp(eventType, 'climb_turn_r')
                radius = params.climb_turn_r.radius;
                angle = params.climb_turn_r.angle;
                gradient = params.climb_turn_r.gradient;
                turnDir = 1;
                dirStr = 'right';
            else
                % Generic climb turn - alternate direction
                radius = params.climb_turn.radius;
                angle = params.climb_turn.angle;
                gradient = params.climb_turn.gradient;
                
                if lastTurnDir == 0
                    if rand() > 0.5
                        turnDir = 1;
                        dirStr = 'right';
                    else
                        turnDir = -1;
                        dirStr = 'left';
                    end
                else
                    turnDir = -lastTurnDir;
                    if turnDir == 1
                        dirStr = 'right';
                    else
                        dirStr = 'left';
                    end
                end
            end
            
            lastTurnDir = turnDir;
            
            % Calculate arc length and elevation change
            arcLength = radius * deg2rad(angle);
            elevChange = arcLength * (gradient / 100);
            
            % Turn geometry (same as regular turn)
            centerDir = currentHeading - turnDir * 90;
            turnCenterX = currentX + radius * cosd(centerDir);
            turnCenterY = currentY + radius * sind(centerDir);
            startAngle = currentHeading + turnDir * 90;
            
            % Generate arc waypoints with increasing elevation
            numArcPoints = max(10, round(angle / 5));
            
            for j = 1:numArcPoints
                progress = j / numArcPoints;
                arcAngle = startAngle - turnDir * angle * progress;
                
                wpX = turnCenterX + radius * cosd(arcAngle);
                wpY = turnCenterY + radius * sind(arcAngle);
                wpZ = currentZ + elevChange * progress;  % Climb during turn
                
                wp.x = wpX;
                wp.y = wpY;
                wp.z = wpZ;
                wp.type = 'climb_turn';
                waypoints(end+1) = wp;
            end
            
            % Update position
            currentX = waypoints(end).x;
            currentY = waypoints(end).y;
            currentZ = waypoints(end).z;
            currentHeading = currentHeading - turnDir * angle;
            
            % Store actual parameters
            generatedEvents(i).type = 'climb_turn';
            generatedEvents(i).actualParams.radius = radius;
            generatedEvents(i).actualParams.angle = angle;
            generatedEvents(i).actualParams.gradient = gradient;
            generatedEvents(i).actualParams.direction = dirStr;
            generatedEvents(i).actualParams.elevChange = elevChange;
            
            fprintf('  Radius: %.0f m, Angle: %.0f deg, Direction: %s, Gradient: %.1f%%, Elev: +%.1f m\n', ...
                radius, angle, dirStr, gradient, elevChange);
            
        case {'descent_turn', 'descent_turn_l', 'descent_turn_r'}
            % Get parameters based on turn type
            if strcmp(eventType, 'descent_turn_l')
                radius = params.descent_turn_l.radius;
                angle = params.descent_turn_l.angle;
                gradient = params.descent_turn_l.gradient;
                turnDir = -1;
                dirStr = 'left';
            elseif strcmp(eventType, 'descent_turn_r')
                radius = params.descent_turn_r.radius;
                angle = params.descent_turn_r.angle;
                gradient = params.descent_turn_r.gradient;
                turnDir = 1;
                dirStr = 'right';
            else
                % Generic descent turn - alternate direction
                radius = params.descent_turn.radius;
                angle = params.descent_turn.angle;
                gradient = params.descent_turn.gradient;
                
                if lastTurnDir == 0
                    if rand() > 0.5
                        turnDir = 1;
                        dirStr = 'right';
                    else
                        turnDir = -1;
                        dirStr = 'left';
                    end
                else
                    turnDir = -lastTurnDir;
                    if turnDir == 1
                        dirStr = 'right';
                    else
                        dirStr = 'left';
                    end
                end
            end
            
            lastTurnDir = turnDir;
            
            % Calculate arc length and elevation change (negative for descent)
            arcLength = radius * deg2rad(angle);
            elevChange = -arcLength * (gradient / 100);
            
            % Turn geometry (same as regular turn)
            centerDir = currentHeading - turnDir * 90;
            turnCenterX = currentX + radius * cosd(centerDir);
            turnCenterY = currentY + radius * sind(centerDir);
            startAngle = currentHeading + turnDir * 90;
            
            % Generate arc waypoints with decreasing elevation
            numArcPoints = max(10, round(angle / 5));
            
            for j = 1:numArcPoints
                progress = j / numArcPoints;
                arcAngle = startAngle - turnDir * angle * progress;
                
                wpX = turnCenterX + radius * cosd(arcAngle);
                wpY = turnCenterY + radius * sind(arcAngle);
                wpZ = currentZ + elevChange * progress;  % Descent during turn
                
                wp.x = wpX;
                wp.y = wpY;
                wp.z = wpZ;
                wp.type = 'descent_turn';
                waypoints(end+1) = wp;
            end
            
            % Update position
            currentX = waypoints(end).x;
            currentY = waypoints(end).y;
            currentZ = waypoints(end).z;
            currentHeading = currentHeading - turnDir * angle;
            
            % Store actual parameters
            generatedEvents(i).type = 'descent_turn';
            generatedEvents(i).actualParams.radius = radius;
            generatedEvents(i).actualParams.angle = angle;
            generatedEvents(i).actualParams.gradient = gradient;
            generatedEvents(i).actualParams.direction = dirStr;
            generatedEvents(i).actualParams.elevChange = elevChange;
            
            fprintf('  Radius: %.0f m, Angle: %.0f deg, Direction: %s, Gradient: %.1f%%, Elev: %.1f m\n', ...
                radius, angle, dirStr, gradient, elevChange);
    end
    fprintf('\n');
end

%% Calculate total track length
totalLength = 0;
for i = 2:length(waypoints)
    dx = waypoints(i).x - waypoints(i-1).x;
    dy = waypoints(i).y - waypoints(i-1).y;
    totalLength = totalLength + sqrt(dx^2 + dy^2);
end

fprintf('=== Track Summary ===\n');
fprintf('Total waypoints: %d\n', length(waypoints));
fprintf('Total length: %.0f m\n', totalLength);
fprintf('Start: (%.0f, %.0f)\n', waypoints(1).x, waypoints(1).y);
fprintf('End:   (%.0f, %.0f)\n', waypoints(end).x, waypoints(end).y);
fprintf('Elevation range: %.1f to %.1f m (relative)\n', ...
    min([waypoints.z]), max([waypoints.z]));

%% Check track bounds (informational only - track will be drawn regardless)
allX = [waypoints.x];
allY = [waypoints.y];
xExtent = [min(allX), max(allX)];
yExtent = [min(allY), max(allY)];

fprintf('Track X extent: %.0f to %.0f (terrain: %.0f to %.0f)\n', xExtent(1), xExtent(2), xMinTerrain, xMaxTerrain);
fprintf('Track Y extent: %.0f to %.0f (terrain: %.0f to %.0f)\n', yExtent(1), yExtent(2), yMinTerrain, yMaxTerrain);

if xExtent(1) < xMinTerrain || xExtent(2) > xMaxTerrain || yExtent(1) < yMinTerrain || yExtent(2) > yMaxTerrain
    fprintf('NOTE: Track extends outside terrain - more tiles needed!\n');
else
    fprintf('Track fits within terrain bounds.\n');
end

%% Save track geometry
trackGeometry.waypoints = waypoints;
trackGeometry.generatedEvents = generatedEvents;
trackGeometry.eventSequence = eventSequence;
trackGeometry.corridorWidth = corridorWidth;
trackGeometry.randomSeed = randomSeed;
trackGeometry.totalLength = totalLength;
trackGeometry.xMinTerrain = xMinTerrain;
trackGeometry.xMaxTerrain = xMaxTerrain;
trackGeometry.yMinTerrain = yMinTerrain;
trackGeometry.yMaxTerrain = yMaxTerrain;
trackGeometry.tifBaseNames = tifBaseNames;

geometryFile = fullfile(editedFolder, 'track_geometry.mat');
save(geometryFile, 'trackGeometry');
fprintf('\nSaved: %s\n', geometryFile);

%% Visualize: Contour plot with track overlay
fprintf('\n--- Generating Visualization ---\n');

% Merge tiles for visualization
cellSize = geoInfo{1}.PixelScale(1);
nColsTotal = round((xMaxTerrain - xMinTerrain) / cellSize);
nRowsTotal = round((yMaxTerrain - yMinTerrain) / cellSize);
mergedElevation = NaN(nRowsTotal, nColsTotal);

for i = 1:numTiles
    bbox = geoInfo{i}.BoundingBox;
    tileXMin = bbox(1,1);
    tileYMax = bbox(2,2);
    
    colStart = round((tileXMin - xMinTerrain) / cellSize) + 1;
    rowStart = round((yMaxTerrain - tileYMax) / cellSize) + 1;
    
    [tileRows, tileCols] = size(elevationData{i});
    rowEnd = min(rowStart + tileRows - 1, nRowsTotal);
    colEnd = min(colStart + tileCols - 1, nColsTotal);
    
    mergedElevation(rowStart:rowEnd, colStart:colEnd) = ...
        elevationData{i}(1:(rowEnd-rowStart+1), 1:(colEnd-colStart+1));
end

% Create coordinate grids
[X, Y] = meshgrid(linspace(xMinTerrain, xMaxTerrain, nColsTotal), ...
                   linspace(yMaxTerrain, yMinTerrain, nRowsTotal));

%% Contour plot
if ~isempty(figHandle) && isvalid(figHandle)
    close(figHandle);
end
figHandle = figure('Name', sprintf('Generated Track (Seed: %d)', randomSeed), 'Position', [200 200 900 700]);
contourf(X, Y, mergedElevation, 20);
colormap(parula);
cb = colorbar; cb.Label.String = 'Elevation (m)';
xlabel('Easting (m)');
ylabel('Northing (m)');
title(sprintf('Generated Ground Track (Seed: %d)', randomSeed));
axis equal tight;
hold on;

% Define colors for each event type
eventColors = struct(...
    'start', [0.2 0.8 0.2], ...           % Green
    'straight', [0.2 0.4 1.0], ...        % Blue
    'turn', [1.0 0.0 1.0], ...            % Magenta
    'climb', [1.0 0.5 0.0], ...           % Orange
    'descent', [0.0 0.8 0.8], ...         % Cyan
    'transition_up', [0.6 0.6 0.6], ...   % Gray (blend back up)
    'transition_down', [0.6 0.6 0.6]);    % Gray (blend back down)

% Plot track segments with color per event type
legendHandles = [];
legendLabels = {};
plottedTypes = {};

for i = 2:length(waypoints)
    wp1 = waypoints(i-1);
    wp2 = waypoints(i);
    
    % Use the type of the destination waypoint for segment color
    segType = wp2.type;
    if isfield(eventColors, segType)
        c = eventColors.(segType);
    else
        c = [0.5 0.5 0.5];  % Gray for unknown
    end
    
    % Plot segment
    h = plot([wp1.x, wp2.x], [wp1.y, wp2.y], '-', 'LineWidth', 3, 'Color', c);
    
    % Add to legend only once per type
    if ~ismember(segType, plottedTypes)
        plottedTypes{end+1} = segType;
        legendHandles(end+1) = h;
        legendLabels{end+1} = upper(segType);
    end
end

    % Plot start and end markers
    plot(waypoints(1).x, waypoints(1).y, 'o', 'MarkerSize', 14, 'MarkerFaceColor', [0 0.8 0], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2, 'DisplayName', 'Start');
    plot(waypoints(end).x, waypoints(end).y, 's', 'MarkerSize', 14, 'MarkerFaceColor', [0.8 0 0], ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2);

    %% User Decision: Reroll or Save
    fprintf('\n=== Track Generated (Seed: %d) ===\n', randomSeed);
    fprintf('  [R] Reroll - generate new random track\n');
    fprintf('  [S] Save   - keep this track and continue\n');
    userChoice = input('Your choice: ', 's');
    
    if strcmpi(userChoice, 'S')
        % Save and exit
        fprintf('Saving track with seed %d...\n', randomSeed);
        userDone = true;
    else
        % Reroll: increment seed, loop again (figure closed at start of next iteration)
        randomSeed = randomSeed + 1;
        fprintf('\nRerolling with new seed %d...\n\n', randomSeed);
    end
    
end  % End of while ~userDone loop

fprintf('\n=== Complete ===\n');
fprintf('Next step: Run edit_elevation.m to carve the track into terrain.\n');

%% ========================================================================
%  FUNCTION: checkSelfIntersection
%  Checks if track has any self-intersections
%  ========================================================================
function hasIntersection = checkSelfIntersection(waypoints, minSeparation)
    hasIntersection = false;
    n = length(waypoints);
    
    if n < 4
        return;
    end
    
    % Check each segment against all non-adjacent segments
    for i = 1:(n-1)
        p1 = [waypoints(i).x, waypoints(i).y];
        p2 = [waypoints(i+1).x, waypoints(i+1).y];
        
        for j = (i+2):(n-1)
            % Skip adjacent segments
            if j == i+1 || j == i-1
                continue;
            end
            
            p3 = [waypoints(j).x, waypoints(j).y];
            p4 = [waypoints(j+1).x, waypoints(j+1).y];
            
            % Check if segments intersect or come too close
            if segmentsIntersect(p1, p2, p3, p4) || ...
               minDistBetweenSegments(p1, p2, p3, p4) < minSeparation
                hasIntersection = true;
                return;
            end
        end
    end
end

%% ========================================================================
%  FUNCTION: segmentsIntersect
%  Checks if two line segments intersect
%  ========================================================================
function intersects = segmentsIntersect(p1, p2, p3, p4)
    % Using cross product method
    d1 = crossDir(p3, p4, p1);
    d2 = crossDir(p3, p4, p2);
    d3 = crossDir(p1, p2, p3);
    d4 = crossDir(p1, p2, p4);
    
    if ((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) && ...
       ((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))
        intersects = true;
    else
        intersects = false;
    end
    
    function d = crossDir(pi, pj, pk)
        d = (pk(1) - pi(1)) * (pj(2) - pi(2)) - (pj(1) - pi(1)) * (pk(2) - pi(2));
    end
end

%% ========================================================================
%  FUNCTION: minDistBetweenSegments
%  Returns minimum distance between two line segments
%  ========================================================================
function dist = minDistBetweenSegments(p1, p2, p3, p4)
    % Simplified: check distance from each endpoint to the other segment
    d1 = pointToSegmentDist(p1, p3, p4);
    d2 = pointToSegmentDist(p2, p3, p4);
    d3 = pointToSegmentDist(p3, p1, p2);
    d4 = pointToSegmentDist(p4, p1, p2);
    dist = min([d1, d2, d3, d4]);
end

function dist = pointToSegmentDist(p, a, b)
    ab = b - a;
    ap = p - a;
    t = max(0, min(1, dot(ap, ab) / dot(ab, ab)));
    projection = a + t * ab;
    dist = norm(p - projection);
end

%% ========================================================================
%  FUNCTION: getTerrainElevAtPoint
%  Gets terrain elevation at a specific X,Y coordinate from tile data
%  ========================================================================
function elev = getTerrainElevAtPoint(x, y, elevationData, geoInfo)
    elev = NaN;
    numTiles = length(elevationData);
    
    for i = 1:numTiles
        bbox = geoInfo{i}.BoundingBox;
        xMin = bbox(1,1); xMax = bbox(2,1);
        yMin = bbox(1,2); yMax = bbox(2,2);
        
        if x >= xMin && x <= xMax && y >= yMin && y <= yMax
            % Point is in this tile
            [rows, cols] = size(elevationData{i});
            cellSizeX = (xMax - xMin) / cols;
            cellSizeY = (yMax - yMin) / rows;
            
            col = round((x - xMin) / cellSizeX) + 1;
            row = round((yMax - y) / cellSizeY) + 1;  % Y is inverted in raster
            
            col = max(1, min(cols, col));
            row = max(1, min(rows, row));
            
            elev = elevationData{i}(row, col);
            return;
        end
    end
    
    % If not found, return NaN (or could extrapolate)
    warning('Point (%.1f, %.1f) not found in any tile', x, y);
end