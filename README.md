# bem3_aw_b2
This is the three-dimensional acoustic wave scattering analysis program for arbitrary objects irradiated by arbitrary beams. 
This is the extension of "bem3_aw_b1" using iterative solution. 
This program can analyze multiple scattering between objects with less memory than "bem3_aw_b1".
Intel Math Kernel Library, Gnu Scientific Library and libpng are required. 
Gmsh is used to create mesh data for objects.
The acoustic wave analysis program "multi_aw" is used for analyze incident field. 


## Usage of example code  
1. type 'make' command to compile.  
   The executable aw_d3b2_create_matrix, aw_d3b2_bv_solver, example1.out, example2.out, example3.out, example4.out are created. 
   The executable aw_d3b2_create_matrix is the solver of boundary integral equations, it outputs coefficient matrices and its inverse matrices. 
   The executable aw_d3b2_bv_solver is the sovler for boundary value, it analyzes multiple scattering using iterative solution when defined multiple objects. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "bem3_aw_b2". 
   The example2.out is the execubable of source code example2.c, it shows a example of sound pressure intensity analysis. 
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of sound pressure as an image.  
   
2. type './aw_d3b2_create_matrix' with arguments of medium datafile name, mesh datafile name, output datafile name [rotation and translation settings (optional)].  
   For example, './aw_d3b2_create_matrix medium1.txt sphere_m1_05z.msh sph_m1_05z_md1.obj'.
   The medium1.txt is the sample of medium datafile, one medium is defined in it. The domain numbers are assigned to the medium from 1 in order. 
   The sphere_m1_05z.msh is the sample of mesh datafile, it is a compressed sphere in the z-direction.
   It was created by Gmsh geometry file sphere_m1.geo in the mesh_sample folder and the stretch program of a object in the mesh_sample/stretch. 
   The sphere_m1_05z_image.png is the visualization result of the sphere_m1_05z.msh (using Gmsh). 
   The mfb.txt is the sample of incident field datafile, a focused beam is defined in it. Please refer to the "multi_aw" for detail. 
   The aw_d3b2_create_matrix outputs datafile for the object, coefficient matrices and its inverse matrices (binary file with the .cmat file extention) with the specified datafile name.
   The aw_d3b2_bv_solver has optional arguments for rotation and translation of the object. 
   When the vector defining rotation axis is (rx, ry, rz), the rotation angle is theta, the translation vector is (tx, ty, tz), 
   the arguments are './aw_d3b2_bv_solver medium1.txt sphere_m1_05z.msh sph_m1_05z_md1.obj rx ry rz theta tx ty tz' (these arguments are real number). 
   Rodrigues' rotation formula is used.  
   
3. type './aw_d3b2_bv_solver' with arguments of model setting datafile name, output datafile name.  
   For example, './aw_d3b2_bv_solver model_settings.txt ex.dat'. 
   The model_settings.txt is the sample of model setting datafile, three-objects are defined in it. 
   The objects are set by using the output datafile name and an additional translation vector. 
   The incident field can be changed by replacing the datafile, except the density, wave speed and frequency settings.
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this case, the file "ex.particles" is output, and the visualization result is "particle.png" (using gnuplot script gscript_particles.plt).  
   
4. type './example1.out' with an argument of datafile name output by aw_d3b2_bv_solver.  
   For example, './example1.out ex.dat'. 
   This executable calculates sound pressure, radiation force and torque.

5. type './example2.out' with an arguemnt of datafile name output by aw_d3b2_bv_solver.  
   For example, './example2.out ex.dat'. 
   This executable calculates sound pressure intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of intensity distributions, created by Gnuplot script gscript_example2.plt 
   (using ImageMagick to convert eps to png).
   
6. type './example3.out' with an argument of datafile name output by aw_d3b2_bv_solver.  
   For example, './example3.out ex.dat'. 
   This executable calculates instantaneous value of the sound pressure, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component, and number of time steps (ex. xz_p_014.png). 
   The color bar is output as color_bar.png in the same folder. 
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane). 
   The xz_p.gif, yz_p.gif and xy_p.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  
   
7. type './example4.out' with an argument of datafile name output by aw_d3b2_bv_solver.  
   For example, './example4.out ex.dat'. 
   This executable calculates scattered field (sound puressure) intensity distributions in far-field and outputs them to text files. 
   The I_example4.png is the visualization result of the intensity distributions, created by gnuplot script gscript_example4.plt. 
   The I_example4_3d.png is the visualization result of the intensity distributions on a spherical surface, created by gnuplot script gscript_example4_3d.plt.

Please see d3b2_src/bem3_aw_b2.h for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS.

![mesh image](sphere_m1_05z_image.png "mesh image of the object (sphere_m1_05z_image.png)") 
![objects](particles.png "nodes for surface integral (particles.png)")  
![intensity](I_example2.png "intensity distributions (I_example2.png)")  
![xz_p.gif](xz_p.gif "instantaneous value of the p on y=0 plane (xz_p.gif)")![xy_p.gif](xy_p.gif "instantaneous value of the p on z=0 plane (xy_p.gif)")  
![far-field](I_example4_3d.png "far-field intensity distribution (I_example4_3d.png)")  

## Analysis sample 2 (in the folder analysis_sample2)  

This is the analysis result of plane wave scattering by the five objects symetrically arranged in z=0 plane. 
The object data is the same as the above example. 

![mesh image 2](analysis_sample2/sphere_m1_05z_image.png "mesh image of the object (analysis_sample2/sphere_m1_05z_image.png)")  
![objects 2](analysis_sample2/particles.png "nodes for surface integral (analysis_sample2/particles.png)")  
![intensity 2](analysis_sample2/I_example2.png "intensity distributions (analysis_sample2/I_example2.png)")  
![xz_p.gif 2](analysis_sample2/xz_p.gif "instantaneous value of the p on y=0 plane (analysis_sample2/xz_p.gif)")![xy_p.gif 2](analysis_sample2/xy_p.gif "instantaneous value of the p on z=0 plane (analysis_sample2/xy_p.gif)")  


## Verification  

The verification results are in the folder verification. 
This is the analysis results of scattering by three spheres arranged in z-direction.
The result of using "aw_msp_ivf" is in the folder aw_msp_ivf_result.


## About mesh file  

This code can use quadrangular (bi-linear) and triangular (linear triangular) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh data are in the folder mesh_sample. 
The file with extension .geo is the Gmsh geometry file. 
The file with extension .msh is the mesh file created by using Gmsh geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 -format msh2 xxxx.geo' in command line (xxxx.geo is a geometry file). 
The domain number (Physical Surface) 99 is assigned to the open region in Gmsh geometry file, 
because Gmsh can't use the number 0 (assigned to open region in the code). 
Please refer to the manual of Gmsh for detail of geometry file.


## References  

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)  
2. GNU Scientific Library [GSL](https://www.gnu.org/software/gsl/)  
3. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
4. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
5. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
6. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
7. The sound pressure analysis program [multi_aw](https://github.com/akohta/multi_aw)  
8. The acoustic wave scattering analysis program [aw_msp_ivf](https://github.com/akohta/aw_msp_ivf/)  
