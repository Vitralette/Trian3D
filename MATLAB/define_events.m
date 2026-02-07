%% Define Course Structure
% This script defines the PARAMETERS for each event type (fixed values)
% and a POOL of events. The sequence is randomized by generate_groundtrack.m
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc; close all;

%% Configuration
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ========================================================================
%  EVENT POOL - Define which events can appear and how many
%  ========================================================================
%  The sequence will be randomly generated from this pool.
%  Format: {'event_type', count}
%  
%  Available events:
%    'straight'  - Level straight segment
%    'turn'      - Horizontal turn (alternates L/R automatically)
%    'turn_l'    - Force left turn
%    'turn_r'    - Force right turn  
%    'climb'     - Ascending straight segment
%    'descent'   - Descending straight segment
%  ========================================================================

eventPool = {
    'straight', 3      % 3 straight segments
    'turn',     2      % 2 turns (will alternate direction)
    'climb',    1      % 1 climb
    'descent',  1      % 1 descent
};

% Always start and end with straight? (recommended for clean entry/exit)
forceStartStraight = true;
forceEndStraight = true;

%  ========================================================================
%  End of pool definition
%  ========================================================================

%% Fixed Parameters for Each Event Type (no min/max - just fixed values)

params.straight.length = 175;           % meters

params.turn.radius = 100;               % meters
params.turn.angle = 40;                 % degrees

params.turn_l.radius = 100;
params.turn_l.angle = 40;
params.turn_l.direction = 'left';

params.turn_r.radius = 100;
params.turn_r.angle = 40;
params.turn_r.direction = 'right';

params.climb.length = 150;              % meters
params.climb.gradient = 8;              % percent

params.descent.length = 150;            % meters
params.descent.gradient = 8;            % percent

%% Corridor Width
corridorWidth = 50;  % meters

%% Transition Segment (blends back to original terrain at end of course)
transition.maxGradient = 1;     % percent (maximum slope for transition)
transition.minLength = 100;     % meters (minimum transition length)

%% Display Configuration
fprintf('=== Course Structure Definition ===\n\n');

fprintf('Event Pool:\n');
totalEvents = 0;
for i = 1:size(eventPool, 1)
    fprintf('  %s: %d\n', eventPool{i,1}, eventPool{i,2});
    totalEvents = totalEvents + eventPool{i,2};
end
fprintf('  Total events: %d\n\n', totalEvents);

fprintf('Fixed Parameters:\n');
fprintf('  Straight: %.0f m\n', params.straight.length);
fprintf('  Turn: %.0f m radius, %.0f deg\n', params.turn.radius, params.turn.angle);
fprintf('  Climb: %.0f m, %.0f%% gradient\n', params.climb.length, params.climb.gradient);
fprintf('  Descent: %.0f m, %.0f%% gradient\n', params.descent.length, params.descent.gradient);
fprintf('  Corridor width: %.0f m\n', corridorWidth);
fprintf('  Transition: max %.0f%% gradient, min %.0f m length\n\n', transition.maxGradient, transition.minLength);

fprintf('Sequence rules:\n');
fprintf('  Force start with straight: %s\n', mat2str(forceStartStraight));
fprintf('  Force end with straight: %s\n', mat2str(forceEndStraight));

%% Save Course Structure
eventStructure.eventPool = eventPool;
eventStructure.params = params;
eventStructure.corridorWidth = corridorWidth;
eventStructure.transition = transition;
eventStructure.forceStartStraight = forceStartStraight;
eventStructure.forceEndStraight = forceEndStraight;

structureFile = fullfile(outputFolder, 'event_structure.mat');
save(structureFile, 'eventStructure');

fprintf('\n=== Structure Saved ===\n');
fprintf('File: %s\n', structureFile);
fprintf('\nNext step: Run generate_groundtrack.m to create a random sequence.\n');
fprintf('Change the seed in generate_groundtrack.m to "reroll" the sequence.\n');
