# Trian3D MATLAB Track Generator

A MATLAB toolchain for generating randomized race tracks with terrain modification and obstacle placement for Trian3D.

## Workflow

### 1. Import Terrain Data
**Script:** `import_geotiff.m`

- Select tiles in Trian3D you want to place the track on
- Export the GeoTIFF files from Trian3D
- Run `import_geotiff.m` to convert `.tif` files to MATLAB `.mat` format
- Outputs: `*_geotiffinfo.mat` and `*_readgeoraster.mat` in `Raw/` folder

### 2. Define Course Events
**Script:** `define_events.m`

Configure the event pool and parameters:
- **Event types:** `straight`, `turn`, `climb`, `descent`, `climb_turn`, `descent_turn`
- **Forced direction variants:** `turn_l`, `turn_r`, `climb_turn_l`, `climb_turn_r`, `descent_turn_l`, `descent_turn_r`
- Set corridor width, turn radius, gradients, etc.
- Outputs: `event_structure.mat`

### 3. Generate Ground Track
**Script:** `generate_groundtrack.m`

- Randomly generates track from the event pool
- **Interactive:** View the track, then `[R]` Reroll or `[S]` Save
- Change seed to get different track layouts
- Outputs: `track_geometry.mat`

### 4. Edit Terrain Elevation
**Script:** `edit_elevation.m`

- Carves the track corridor into the terrain
- Climbs raise terrain, descents lower terrain
- Smooth blending at corridor edges
- Transition segment blends back to original terrain at track end
- Outputs: `*_edited.mat` files

### 5. Generate Obstacles (Optional)

#### Trees along track boundaries
**Script:** `generate_point_trees_kml.m`
- Places trees at regular intervals along left/right track boundaries
- Outputs: `point_trees.kml`

#### Powerline obstacles
**Script:** `generate_linestring_powerline_kml.m`
- Places powerlines orthogonally across the track
- **Interactive:** `[R]` Reroll or `[S]` Save placement
- Outputs: `linestring_powerline.kml`

### 6. Export Modified Terrain
**Script:** `export_geotiff.m`

- Converts edited `.mat` files back to GeoTIFF format
- Outputs: `*_edited.tif` files

### 7. Import into Trian3D
- Load the edited `.tif` files into Trian3D
- Import the generated `.kml` files for obstacles

## Folder Structure

```
TRIAN3D/SampleProject/
├── Raw/           # Original terrain data (.mat)
├── Edited/        # Modified terrain + track geometry (.mat)
└── Export/        # KML files for obstacles
```

## Reference Files

- `triantounrealtest.kml` - Example KML showing supported object types (trees, powerline, stream, etc.)
