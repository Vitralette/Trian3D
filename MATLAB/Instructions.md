
# Trian3D in Matlab 

What we want to do 

    1. We have a sample project in Trian3D

    2. We want to take the Trian3D files of the sample project 
        - tif
            - these files represent digital elevation data 
            - various tif files generate a complete elevation surface 

    3. load in the tif files with matlab
         - they are located at TRIAN3D/SampleProject      
         - we also have the corresponding matlab object obtained by using geotiffinfo command in matlab, but since I dont have the toolbox we dont want to use the command and only make use of the already generated obejctes i obtained   
         - we also have the corresponding matlab object obtained by using readgeoraster command in matlab, but since I dont have the toolbox we dont want to use the command and only make use of the already generated obejctes i obtained 
         - decide which data which can use from which object or the original file to achive the goal  
         
    4. We need a matlab script that reproduces the eleation data in a surface plot or something similar 

DONE 


Next Steps

1. From the Trian Project we additional "objects" as .vec , those a are lines drawn onto the elevation data 
2. in the next step we want to display those lines ontop of the elevation data 


DONE 



lets build a test modification to the tiles 



Next i want to implement the editing process ,

lets create a new script for an example edit, i want to take the two tiles , then define a groundtrack from bottom left to top  right so diagonal line across the two tiles, then we define a width of 50m which we apply to the line, so we have a "racetrack" going from bottom left to top right , then we want to edit the eleveantion data in that racetrack to create a canyon , so reduce the elevation form the orginal data by -100 m 

the flow should be run trian3dhack -> apply the canyon changes (produces _edited files) -> then also plot the edited file to visaulize the changes -> then we can run the export_geotiff to produce the tiff we can load in trian3d