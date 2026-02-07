%% Define Course Structure
% This script defines the PARAMETERS for each event type (fixed values)
% and a POOL of events. The sequence is randomized by generate_groundtrack.m
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc; close all;

%% Load Project Configuration
config = load_project_config();
outputFolder = config.editedFolder;

fprintf('Project: %s\n', config.projectName);

%% ========================================================================
%  EVENT POOL - Define which events can appear and how many
%  ========================================================================
%  The sequence will be randomly generated from this pool.
%  Format: {'event_type', count}
%  
%  Available events:
%    'straight'     - Level straight segment
%    'turn'         - Horizontal turn (alternates L/R automatically)
%    'turn_l'       - Force left turn
%    'turn_r'       - Force right turn  
%    'climb'        - Ascending straight segment
%    'descent'      - Descending straight segment
%    'climb_turn'   - Climbing turn (alternates L/R automatically)
%    'climb_turn_l' - Climbing left turn
%    'climb_turn_r' - Climbing right turn
%    'descent_turn' - Descending turn (alternates L/R automatically)
%    'descent_turn_l' - Descending left turn
%    'descent_turn_r' - Descending right turn
%  ========================================================================

eventPool = {
    'straight', 3      % 3 straight segments
    'turn',     4      % 4 turns (will alternate direction)
    'climb',    1      % 1 climb
    'descent',  1      % 1 descent
    'climb_turn', 2    % 2 climbing turns
    'descent_turn', 2  % 2 descending turns
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

% Climbing turn parameters (turn + climb combined)
params.climb_turn.radius = 100;         % meters
params.climb_turn.angle = 40;           % degrees
params.climb_turn.gradient = 6;         % percent (climb during turn)

params.climb_turn_l.radius = 100;
params.climb_turn_l.angle = 40;
params.climb_turn_l.gradient = 6;
params.climb_turn_l.direction = 'left';

params.climb_turn_r.radius = 100;
params.climb_turn_r.angle = 40;
params.climb_turn_r.gradient = 6;
params.climb_turn_r.direction = 'right';

% Descending turn parameters (turn + descent combined)
params.descent_turn.radius = 100;       % meters
params.descent_turn.angle = 40;         % degrees
params.descent_turn.gradient = 6;       % percent (descent during turn)

params.descent_turn_l.radius = 100;
params.descent_turn_l.angle = 40;
params.descent_turn_l.gradient = 6;
params.descent_turn_l.direction = 'left';

params.descent_turn_r.radius = 100;
params.descent_turn_r.angle = 40;
params.descent_turn_r.gradient = 6;
params.descent_turn_r.direction = 'right';

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
fprintf('  Climb Turn: %.0f m radius, %.0f deg, %.0f%% gradient\n', params.climb_turn.radius, params.climb_turn.angle, params.climb_turn.gradient);
fprintf('  Descent Turn: %.0f m radius, %.0f deg, %.0f%% gradient\n', params.descent_turn.radius, params.descent_turn.angle, params.descent_turn.gradient);
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
