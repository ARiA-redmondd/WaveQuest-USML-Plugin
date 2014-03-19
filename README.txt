WaveQuest-USML-Plugin
=====================


This plugin provides a surface reflection loss model and a convenient way to interface data in Cartesian coordinates with USML. In this folder are the following subdirectories: 

boundary_grid_quad

data_grid_quad

quadtree

reflect_loss_hybrid

Interfacing USML and Cartesian data is accomplished through the use of quadtrees and a data_grid-like class called data_grid_quad. The class files needed to implement a quadtree are included in the quadtree directory. These files are quadtree.cc and quadtree.h. Two other files, wedge_cartesian.txt and quad_test.cc, are also included in the quadtree folder to USML's build of the quadtree class. Both quadtree.cc and quadtree.h should be placed in the usml/types directory and wedge_cartesian.txt and quad_test.cc should be placed in the usml/types/test directory. Note that whenever any files are added or removed from your USML build, a new makefile must be generated and USML must be rebuilt with the new makefile. 

The class data_grid_quad is much like data_grid, except that it computes nearest neighbors by querying the quadtree that stores the Cartesian data. These neighbors are then used for interpolation. The data_grid_quad subdirectory contains the class file for data_grid_quad (data_grid_quad.h) and a file for testing USML's build of this class (datagrid_quad_test.cc). The .h file must be copied into the usml/types directory and the types.h file must be modified to include the line "#include \<usml/types/data_grid_quad.h\>" after the line "#include <usml/types/data_grid_svp.h>". The .cc file must be copied into usml/types/test. 

The boundary_grid_quad directory contains the code for the boundary_grid_quad class, which implements the data_grid_quad as a boundary_model. The boundary_grid_quad class is implemented in boundary_grid_quad.h and its build in USML is tested in boundary_test.cc. The .h file must be put in the usml/ocean directory and the ocean.h file must be modified to include the line "#include <usml/ocean/boundary_grid_quad.h>" after the line "#include <usml/ocean/boundary_grid_fast.h>". The .cc file must be copied into usml/ocean/test. It will replace an existing boundary_test.cc file. This is okay since it's the same as the previous one, but with the addition of unit tests for boundary_grid_quad. 

The reflect_loss_hybrid folder contains code to implement and test an implementation of a Modified Eckhart / APL-UW surface loss model. The files reflect_loss_hybrid.h and reflect_loss_hybrid.cc are class files that implement this model. They must be copied into usml/ocean and the line "#include <usml/ocean/reflect_loss_hybrid.h>" must be added to ocean.h after the line "#include <usml/ocean/reflect_loss_rayleigh.h>". The file refect_loss_test.cc must be put in usml/ocean/test. It will replace an existing boundary_test.cc file. This is okay since it's the same as the previous one, but with the addition of unit tests for reflect_loss_hybrid. 
