function config = load_project_config()
%LOAD_PROJECT_CONFIG Load the current project configuration
%   config = load_project_config() returns the project configuration
%   structure containing all paths. If no config exists, it displays
%   an error message instructing the user to run define_project.m first.
%
%   Fields returned:
%       config.projectName   - Name of the active project
%       config.baseFolder    - Base TRIAN3D folder
%       config.projectFolder - Project root folder
%       config.rawFolder     - Raw data folder
%       config.editedFolder  - Edited data folder
%       config.exportFolder  - Export folder (KML files)
%       config.timestamp     - When config was last saved
%
%   Author: Tim Jusko
%   Date: 2026-02-07

    configFile = 'project_config.mat';
    
    if ~exist(configFile, 'file')
        error(['Project configuration not found.\n' ...
               'Please run define_project.m first to set the active project.']);
    end
    
    data = load(configFile);
    config = data.projectConfig;
    
    % Verify project folder still exists
    if ~exist(config.projectFolder, 'dir')
        warning('Project folder no longer exists: %s', config.projectFolder);
    end
end
