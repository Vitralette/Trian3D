%% Define Project Configuration
% This script sets the active project name and saves it to a configuration
% file. All other scripts will read this configuration to determine paths.
%
% Usage: Run this script first to set which project you're working on.
%        Then run any other script (generate_groundtrack, edit_elevation, etc.)
%
% Author: Tim Jusko
% Date: 2026-02-07

clear; clc;

%% Configuration
% Set the project name here - this determines all relative paths
projectName = 'SampleProject';  % Change to 'ARCADE' or other project name

%% Build paths
baseFolder = fullfile('..', 'PROJECTS');
projectFolder = fullfile(baseFolder, projectName);
rawFolder = fullfile(projectFolder, 'Raw');
editedFolder = fullfile(projectFolder, 'Edited');
exportFolder = fullfile(projectFolder, 'Export');

%% Verify project folder exists
if ~exist(projectFolder, 'dir')
    fprintf('WARNING: Project folder does not exist: %s\n', projectFolder);
    fprintf('Creating project folder structure...\n');
    mkdir(projectFolder);
    mkdir(rawFolder);
    mkdir(editedFolder);
    mkdir(exportFolder);
    fprintf('  Created: %s\n', rawFolder);
    fprintf('  Created: %s\n', editedFolder);
    fprintf('  Created: %s\n', exportFolder);
else
    % Create subfolders if they don't exist
    if ~exist(rawFolder, 'dir'), mkdir(rawFolder); end
    if ~exist(editedFolder, 'dir'), mkdir(editedFolder); end
    if ~exist(exportFolder, 'dir'), mkdir(exportFolder); end
end

%% Save configuration
projectConfig.projectName = projectName;
projectConfig.baseFolder = baseFolder;
projectConfig.projectFolder = projectFolder;
projectConfig.rawFolder = rawFolder;
projectConfig.editedFolder = editedFolder;
projectConfig.exportFolder = exportFolder;
projectConfig.timestamp = datetime('now');

configFile = 'project_config.mat';
save(configFile, 'projectConfig');

%% Display confirmation
fprintf('\n========================================\n');
fprintf('  PROJECT CONFIGURATION SAVED\n');
fprintf('========================================\n');
fprintf('  Project Name:   %s\n', projectName);
fprintf('  Project Folder: %s\n', projectFolder);
fprintf('  Raw Folder:     %s\n', rawFolder);
fprintf('  Edited Folder:  %s\n', editedFolder);
fprintf('  Export Folder:  %s\n', exportFolder);
fprintf('  Config File:    %s\n', configFile);
fprintf('========================================\n');
fprintf('\nAll other scripts will now use these paths.\n');
fprintf('Run define_project.m again to switch projects.\n');
