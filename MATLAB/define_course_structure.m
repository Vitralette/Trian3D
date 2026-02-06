%% Define Course Structure
% This script defines the STRUCTURE of a low-level flight course as a
% sequence of events. The actual parameters are randomized within bounds
% by generate_groundtrack.m
%
% Event Types:
%   'straight'  - Level straight segment
%   'turn'      - Horizontal turn (left or right)
%   'climb'     - Ascending straight segment
%   'descent'   - Descending straight segment
%
% Each event has:
%   - type: event type string
%   - params: struct with min/max bounds for randomization
%
% Author: Tim Jusko
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% MTE Constraints (modify based on your helicopter performance data)
mteConstraints.minTurnRadius = 100;      % meters (minimum safe turn radius)
mteConstraints.maxTurnRadius = 300;      % meters
mteConstraints.minClimbGradient = 3;     % percent
mteConstraints.maxClimbGradient = 12;    % percent
mteConstraints.minDescentGradient = 3;   % percent
mteConstraints.maxDescentGradient = 15;  % percent
mteConstraints.minSegmentLength = 150;   % meters
mteConstraints.maxSegmentLength = 600;   % meters
mteConstraints.corridorWidth = 50;       % meters (track width)

%% Define Event Sequence
% This is the FIXED structure - same events in same order every run
% Only the parameters within bounds will be randomized
%
% NOTE: Segment lengths are tuned for ~2km x 1km terrain (2 tiles)
% Adjust if using larger/smaller terrain
%
% DESIGN: S-pattern with alternating turns to prevent self-intersection
% Turn directions alternate: R -> L -> R (or L -> R -> L based on first turn)

events = struct('type', {}, 'params', {});

% Event 1: Initial straight to establish course
events(1).type = 'straight';
events(1).params.lengthMin = 150;
events(1).params.lengthMax = 200;

% Event 2: First turn (random direction, this sets the pattern)
events(2).type = 'turn';
events(2).params.radiusMin = 80;
events(2).params.radiusMax = 120;
events(2).params.angleMin = 30;      % degrees - SMALLER angles for tighter terrain
events(2).params.angleMax = 50;
events(2).params.directionRandom = true;
events(2).params.direction = 'right';

% Event 3: Climb segment
events(3).type = 'climb';
events(3).params.lengthMin = 150;
events(3).params.lengthMax = 200;
events(3).params.gradientMin = 5;    % percent
events(3).params.gradientMax = 10;

% Event 4: Straight at altitude
events(4).type = 'straight';
events(4).params.lengthMin = 100;
events(4).params.lengthMax = 150;

% Event 5: Second turn (OPPOSITE direction to prevent crossing)
events(5).type = 'turn';
events(5).params.radiusMin = 80;
events(5).params.radiusMax = 120;
events(5).params.angleMin = 30;
events(5).params.angleMax = 50;
events(5).params.directionRandom = false;  % FORCED opposite
events(5).params.direction = 'opposite';   % Special: will be set opposite to previous turn

% Event 6: Descent segment
events(6).type = 'descent';
events(6).params.lengthMin = 100;
events(6).params.lengthMax = 150;
events(6).params.gradientMin = 6;
events(6).params.gradientMax = 12;

% Event 7: Low-level straight
events(7).type = 'straight';
events(7).params.lengthMin = 100;
events(7).params.lengthMax = 150;

%% Display Course Structure
fprintf('=== Course Structure ===\n');
fprintf('Number of events: %d\n', length(events));
fprintf('Corridor width: %.0f m\n\n', mteConstraints.corridorWidth);

for i = 1:length(events)
    e = events(i);
    fprintf('Event %d: %s\n', i, upper(e.type));
    switch e.type
        case 'straight'
            fprintf('  Length: %.0f - %.0f m\n', e.params.lengthMin, e.params.lengthMax);
        case 'turn'
            fprintf('  Radius: %.0f - %.0f m\n', e.params.radiusMin, e.params.radiusMax);
            fprintf('  Angle: %.0f - %.0f deg\n', e.params.angleMin, e.params.angleMax);
            if e.params.directionRandom
                fprintf('  Direction: RANDOM\n');
            else
                fprintf('  Direction: %s\n', e.params.direction);
            end
        case {'climb', 'descent'}
            fprintf('  Length: %.0f - %.0f m\n', e.params.lengthMin, e.params.lengthMax);
            fprintf('  Gradient: %.0f - %.0f %%\n', e.params.gradientMin, e.params.gradientMax);
    end
    fprintf('\n');
end

%% Save Course Structure
courseStructure.events = events;
courseStructure.mteConstraints = mteConstraints;
courseStructure.corridorWidth = mteConstraints.corridorWidth;

structureFile = fullfile(outputFolder, 'course_structure.mat');
save(structureFile, 'courseStructure');

fprintf('=== Structure Saved ===\n');
fprintf('File: %s\n', structureFile);
fprintf('\nNext step: Run generate_groundtrack.m to create a random track instance.\n');
fprintf('Change the seed in generate_groundtrack.m to "reroll" the course.\n');
