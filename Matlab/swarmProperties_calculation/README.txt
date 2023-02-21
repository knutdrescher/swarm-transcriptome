The files in this folder were used to calculate all swarm properties obtained from images.
The main function to execute is calculateSwarmProperties, which will call all other functions.

The input of calculateSwarmProperties is the path to all imaged obtained during the experiment.
They are organized in subfolders, one folder for each time point, whose folder name contains
the date and time that this imaging round was started. Folders contain images and videos saved
as tif files with the file names indicating the x-y coordinates on the swarming plate. Segmentation
results obtained using StarDist are saved in the same folders as label image tif stacks.

For practical reasons, property names in the code differ slightly from the names given in the Figures. 
The document PropertyNames.xlsx describes which property name in the code relates to which property
in figures, and also provides conversion factors multiplied to the raw data in order to match the displayed
units.
