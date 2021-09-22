# DICE
Discontinuity Intensity Calculator and Estimator

Hello everybody,

you can use the DICE app in two way:

1) download the installation file and install the app. It will appear onto the app panel of the MATLAB software.
2) download the DICE_application folder, open the DICE.mlapp (with MATLAB) and run it.

The DICE.mlapp is the MATLAB-app file that can be modified. If you are able, feel free to download the DICE app and modify it according to your preferences/needs.
To use the app you will need for the Computer Vision Toolbox of MATLAB

The User Manual is coming soon.

You can find also an Example dataset composed by:

1) a colored 3D point cloud of a small Lastoni di Formin rocky outcrops (presented in Menegoni et al., in review) with 3cm point spacing and saved in the txt format.
2) a DXF files that contains all the 3D fracture trace mapped onto a higher resolution point cloud;
3) an excel file that contains all the geometrical information of the circular fracture planes fitted to the fracture traces using the dxf_planefit.m algorithm (you can find it in the DICE_app folder).
4) an excel file containing the attitude (dip and dip direction) and set belonging of each fitted fracture planes.
5) an excel files containing the mean attitude (dip and dip direction) of the fracture sets.


If you want to use your dataset, you must prepare these files and follow these steps:

1) prepare a subsampled point cloud in txt format that contains X, Y and Z coordinates and R, G and B color value for each point of the cloud. To create it, please consider using CloudCompare (https://www.cloudcompare.org/) importing your point cloud, subsampling it and then saving it. No header must be present in the text file.
2) prepare a DXF file containing all the traces of the fractures mapped onto the point cloud.
3) use the dxf_planefit.m algorithm to fit the circular fracture traces. To do it, launch the algorithm, select the DXF file and the folder where the circular discs will be saved. The excel file that contains all the geometric information of the fitted plane will be saved in the same path of the DXF files.
4) use a stereoplot software (e.g. Stereonet, https://www.rickallmendinger.net/stereonet; Dips, https://www.rocscience.com/software/dips; DIPANALYST, http://www.dipanalyst.com/) to recognize the main fracture sets and create an excel file containing the attitude (as dip and dip direction) and the set belonging for each fracture planes (the order must be the same of the previous excel file).
5) create another excel file containing the mean attitude ( as Dip and Dip Direction) and the number of set for each fracture set.
Please, consider having a look to the example dataset files and then produce yours similarly. N.B. the header of the excel files must be the same of the those used in the example dataset.
Then, run the DICE app, use the buttons to load the previous prepared files and selected and use a sampling strategy (SEE BELOW):
--------------3D SCANLINES---------------------
Sample the scanline/s onto a 3D digital outcrop model using the same software used to sample the fractures. The scanline/s must be polyline/s composed by two points. Save each scanline as a single DXF file in a new empty folder.
Launch the app, load the point cloud and fractures files. Then in the scanline panel load the folder where the scanline/s is/are stored. Then push calculate the discontinuity/fracture information and export them.
The results of the scanline analysis are saved in the scanline folder.

-------------3D CIRCULAR SCAN WINDOW-----------------
Launch the app, load the point cloud and fractures files. Then in the circular scan window panel you can perform the manual or automatic analysis.
The manual analysis is performed by defining the Circular Window (CW) radius, selecting the center of the CW by a point picking mechanism ('Select the center of CW' button) and then calculating the fracture statistic. The statistic will be reported in the MATLAB command window.
The automatic analysis is performed by defining the Circular Window (CW) radius and then calculating the P21 for each point of the cloud ('P21 point cloud' button). The results will be stored in a TXT format point cloud in a new folder placed in the same path of the input point cloud.

-----------3D SCAN SPHERE----------------------
Launch the app, load the point cloud and fractures files. Then, in the Scan Volume panel define the radius of the scan sphere and perform the P32 analysis ('P32 analysis' button).

-----------------------------------------------

Cheers,
Niccolò
