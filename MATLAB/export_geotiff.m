%% Export Elevation Data to GeoTIFF
% This script loads edited elevation data (*_edit.mat) and exports
% to GeoTIFF format that can be loaded into Trian3D.
%
% Workflow:
%   1. Run trian3dhack.m first to load, visualize, edit, and save *_edit.mat
%   2. Run this script to convert *_edit.mat to *_matlab.tif
%
% Author: Generated for Trian3D Sample Project
% Date: 2026-02-06

clear; clc; close all;

%% Configuration
dataFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Edited');
outputFolder = fullfile('..', 'TRIAN3D', 'SampleProject', 'Modified');

% Create output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
    fprintf('Created output folder: %s\n', outputFolder);
end

% Define the base names of the TIF files
tifBaseNames = {
    'dgm1_32_495_5802_1_nw_2023'
    'dgm1_32_496_5802_1_nw_2023'
};

%% Process each tile
for i = 1:length(tifBaseNames)
    baseName = tifBaseNames{i};
    fprintf('\n=== Processing tile %d: %s ===\n', i, baseName);
    
    %% Step 1: Load the edited data (*_geotiffinfo_edited.mat and *_readgeoraster_edited.mat)
    geoInfoFile = fullfile(dataFolder, [baseName '_geotiffinfo_edited.mat']);
    rasterFile = fullfile(dataFolder, [baseName '_readgeoraster_edited.mat']);
    
    if ~exist(geoInfoFile, 'file') || ~exist(rasterFile, 'file')
        warning('Edit files not found for: %s\nRun trian3dhack.m first!', baseName);
        continue;
    end
    
    % Load geotiffinfo
    geoInfoData = load(geoInfoFile);
    geoInfo = geoInfoData.geoInfoEdited;
    
    % Load elevation data
    rasterData = load(rasterFile);
    elevationData = rasterData.A;
    
    fprintf('  Loaded edited data: %d x %d pixels\n', size(elevationData, 1), size(elevationData, 2));
    fprintf('  Elevation range: %.2f to %.2f\n', min(elevationData(:)), max(elevationData(:)));
    
    %% Step 2: Export as GeoTIFF
    outputFile = fullfile(outputFolder, [baseName '_matlab.tif']);
    
    % Write GeoTIFF using low-level Tiff class
    writeGeoTiff(outputFile, elevationData, geoInfo);
    
    fprintf('  Exported to: %s\n', outputFile);
end

fprintf('\n=== Export Complete ===\n');
fprintf('Modified tiles saved to: %s\n', outputFolder);

%% ========================================================================
%  FUNCTION: writeGeoTiff
%  Writes elevation data to a GeoTIFF file with proper geospatial tags
%  ========================================================================
function writeGeoTiff(filename, elevationData, geoInfo)
    % Get data dimensions
    [rows, cols] = size(elevationData);
    
    % Determine the data type and corresponding TIFF sample format
    dataClass = class(elevationData);
    switch dataClass
        case 'single'
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.IEEEFP;
        case 'double'
            elevationData = single(elevationData);  % Convert to single for file size
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.IEEEFP;
        case 'int16'
            bitsPerSample = 16;
            sampleFormat = Tiff.SampleFormat.Int;
        case 'uint16'
            bitsPerSample = 16;
            sampleFormat = Tiff.SampleFormat.UInt;
        case 'int32'
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.Int;
        case 'uint32'
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.UInt;
        otherwise
            elevationData = single(elevationData);
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.IEEEFP;
    end
    
    % Create the TIFF file
    t = Tiff(filename, 'w');
    
    % Set basic TIFF tags
    tagstruct.ImageLength = rows;
    tagstruct.ImageWidth = cols;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = bitsPerSample;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.SampleFormat = sampleFormat;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Compression = Tiff.Compression.None;
    
    % Set row-per-strip (for efficient reading)
    tagstruct.RowsPerStrip = 1;
    
    % Set the standard TIFF tags
    t.setTag(tagstruct);
    
    % Set GeoTIFF tags
    % ModelPixelScaleTag (33550): [ScaleX, ScaleY, ScaleZ]
    if isfield(geoInfo, 'PixelScale') && ~isempty(geoInfo.PixelScale)
        t.setTag('ModelPixelScaleTag', geoInfo.PixelScale);
    end
    
    % ModelTiepointTag (33922): [I, J, K, X, Y, Z]
    % Tiepoint connects raster coordinates (I,J,K) to world coordinates (X,Y,Z)
    if isfield(geoInfo, 'TiePoints') && ~isempty(geoInfo.TiePoints)
        tiepoints = geoInfo.TiePoints;
        % TiePoints structure has fields: ImagePoints and WorldPoints
        if isstruct(tiepoints)
            % Format: [I, J, K, X, Y, Z]
            tiepointTag = [tiepoints.ImagePoints.Row, tiepoints.ImagePoints.Col, 0, ...
                          tiepoints.WorldPoints.X, tiepoints.WorldPoints.Y, 0];
        else
            % If it's already an array, use as-is
            tiepointTag = tiepoints(:)';
        end
        t.setTag('ModelTiepointTag', tiepointTag);
    else
        % Calculate tiepoint from BoundingBox
        % Upper-left corner of the image (pixel 0,0) maps to (xMin, yMax)
        bbox = geoInfo.BoundingBox;
        xMin = bbox(1, 1);
        yMax = bbox(2, 2);
        tiepointTag = [0, 0, 0, xMin, yMax, 0];
        t.setTag('ModelTiepointTag', tiepointTag);
    end
    
    % GeoKeyDirectoryTag (34735): Contains the GeoTIFF keys
    if isfield(geoInfo, 'GeoTIFFTags') && isfield(geoInfo.GeoTIFFTags, 'GeoKeyDirectoryTag')
        geoKeyDirStruct = geoInfo.GeoTIFFTags.GeoKeyDirectoryTag;
        % Build GeoKeyDirectory array from struct
        geoKeyDir = buildGeoKeyDirectoryFromStruct(geoKeyDirStruct);
        t.setTag('GeoKeyDirectoryTag', double(geoKeyDir));
    elseif isfield(geoInfo, 'GeoTIFFCodes')
        % Build GeoKeyDirectory from GeoTIFFCodes
        geoKeyDir = buildGeoKeyDirectory(geoInfo);
        t.setTag('GeoKeyDirectoryTag', double(geoKeyDir));
    end
    
    % Note: GeoDoubleParamsTag (34736) and GeoAsciiParamsTag (34737) are not 
    % supported by MATLAB's Tiff class. The essential spatial information is
    % already stored in ModelPixelScaleTag, ModelTiepointTag, and GeoKeyDirectoryTag.
    
    % Write the elevation data
    t.write(elevationData);
    
    % Close the file
    t.close();
    
    fprintf('    Written GeoTIFF: %d x %d, %s, %d bits\n', cols, rows, dataClass, bitsPerSample);
end

%% ========================================================================
%  FUNCTION: buildGeoKeyDirectoryFromStruct
%  Builds a GeoKeyDirectory array from GeoTIFFTags.GeoKeyDirectoryTag struct
%  ========================================================================
function geoKeyDir = buildGeoKeyDirectoryFromStruct(gkdStruct)
    % GeoKeyDirectory format: 2D array with 4 columns
    % Row 1: [KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys]
    % Row 2+: [KeyID, TIFFTagLocation, Count, Value_Offset]
    
    % GeoTIFF Key IDs
    keyMapping = struct(...
        'GTModelTypeGeoKey', 1024, ...
        'GTRasterTypeGeoKey', 1025, ...
        'GTCitationGeoKey', 1026, ...
        'GeographicTypeGeoKey', 2048, ...
        'GeogCitationGeoKey', 2049, ...
        'GeogGeodeticDatumGeoKey', 2050, ...
        'GeogAngularUnitsGeoKey', 2054, ...
        'GeogEllipsoidGeoKey', 2056, ...
        'GeogSemiMajorAxisGeoKey', 2057, ...
        'GeogSemiMinorAxisGeoKey', 2058, ...
        'ProjectedCSTypeGeoKey', 3072, ...
        'ProjLinearUnitsGeoKey', 3076);
    
    % Collect keys
    keys = [];
    
    % Get field names from the struct
    fields = fieldnames(gkdStruct);
    
    for i = 1:length(fields)
        fieldName = fields{i};
        value = gkdStruct.(fieldName);
        
        % Skip string values (they go to GeoAsciiParamsTag)
        if ischar(value) || isstring(value)
            continue;
        end
        
        % Get the key ID
        if isfield(keyMapping, fieldName)
            keyId = keyMapping.(fieldName);
        else
            continue;  % Unknown key, skip
        end
        
        % Add key row: [KeyID, TIFFTagLocation, Count, Value_Offset]
        % TIFFTagLocation = 0 means value is stored directly
        keys = [keys; keyId, 0, 1, value];
    end
    
    numKeys = size(keys, 1);
    
    % Build 2D array: header row + key rows
    % Header: [KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys]
    geoKeyDir = [1, 1, 0, numKeys; keys];
end

%% ========================================================================
%  FUNCTION: buildGeoKeyDirectory
%  Builds a GeoKeyDirectory from GeoTIFFCodes structure
%  ========================================================================
function geoKeyDir = buildGeoKeyDirectory(geoInfo)
    % GeoKeyDirectory format: 2D array with 4 columns
    % Row 1: [KeyDirectoryVersion, KeyRevision, MinorRevision, NumberOfKeys]
    % Row 2+: [KeyID, TIFFTagLocation, Count, Value_Offset]
    
    codes = geoInfo.GeoTIFFCodes;
    
    % Collect keys
    keys = [];
    
    % GTModelTypeGeoKey (1024): ModelType
    if isfield(codes, 'Model') && ~isempty(codes.Model)
        keys = [keys; 1024, 0, 1, codes.Model];
    end
    
    % GTRasterTypeGeoKey (1025): RasterType (1=PixelIsArea, 2=PixelIsPoint)
    keys = [keys; 1025, 0, 1, 1];  % Assume PixelIsArea
    
    % ProjectedCSTypeGeoKey (3072): PCS code
    if isfield(codes, 'PCS') && ~isempty(codes.PCS)
        keys = [keys; 3072, 0, 1, codes.PCS];
    end
    
    numKeys = size(keys, 1);
    
    % Build 2D array: header row + key rows
    geoKeyDir = [1, 1, 0, numKeys; keys];
end
