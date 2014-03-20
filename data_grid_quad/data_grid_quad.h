/**
 * @file data_grid_quad.h
 */
#ifndef USML_TYPES_DATA_GRID_QUAD_H
#define USML_TYPES_DATA_GRID_QUAD_H

#include <string.h>
#include <usml/types/wvector.h>
#include <usml/types/seq_vector.h>
#include <usml/types/quadtree.h>

using namespace usml::ublas;

namespace usml {
namespace types {

/// @ingroup data_grid_quad
/// @{


/**
 * Supports interpolation.
 *
 * @param  DATA_TYPE    Type of data to be interpolated. Must support +,-,*,/
 *                      with itself and double precision scalars.
 * @param  NUM_DIMS     Number of dimensions in this grid.  Specifying this
 *                      at compile time allows for some loop un-wrapping.
 */
template<class DATA_TYPE, unsigned NUM_DIMS> class data_grid_quad
{


protected:

    /** Axis associated with each dimension of the data grid. */
  //    seq_vector* _axis[NUM_DIMS];

    /**
     * Multi-dimensional data stored as a linear array in column major order.
     * This format is used to support an N-dimensional data set
     * with any number of dimensions.
     * This memory is created in the constructor and deleted in the destructor.
     */
  //    DATA_TYPE *_data;

    /** Used during interpolation to hold the axis offsets. */

    quadtree* _tree;

    double _height_at_origin;

    int _num_boundary_points;
    
    double* _boundary_points;


public:

    /**
     * Multi-dimensional interpolation with the derivative calculation.
     * So many calculations are shared between the determination of an
     * interpolate value and its derivative, that it is computationally
     * efficient to compute them both at the same time.
     *
     * Limit interpolation to axis domain if _edge_limit turned on for that
     * dimension.  Allow extrapolation if _edge_limit turned off.
     *
     * @param   location    Location at which field value is desired. Must
     *                      have the same rank as the data grid or higher.
     *                      WARNING: The contents of the location boost::numeric::ublas::vector
     *                      may be modified if edge_limit() is true for
     *                      any dimension.
     * @param   derivative  If this is not null, the first derivative
     *                      of the field at this point will also be computed.
     * @return              Value of the field at this point.
     */
    DATA_TYPE interpolate(double* location, DATA_TYPE* derivative = NULL)
    {
      double interpolated_x = 0;
      double interpolated_y = 0;
      double interpolated_z;
      double xmin, xmax;
      double ymin, ymax;

      std::vector<double> coords = _tree->inverse_transform(*(location),*(location+1),*(location+2),_height_at_origin);

      xmin = round(10*coords.at(0))*0.1;
      xmax = xmin;
      ymin = round(coords.at(1));
      ymax = ymin;

      _tree->traverse(xmin,ymin,xmax,ymax);
      std::vector<double> neighbors = _tree->Search_Results;

      if (neighbors.size() == 3){
	interpolated_x = neighbors.at(0);
	interpolated_y = neighbors.at(1);
	interpolated_z = neighbors.at(2);
      }

      else {
	if (coords.at(0) <= 3999){
	  xmin = round(10*coords.at(0))*0.1-0.5;
	  xmax = round(10*coords.at(0))*0.1+1;
	  ymin = round(coords.at(1));
	  ymax = ymin;
	}
	
	else {
	  xmin = round(10*coords.at(0))*0.1-1.5;
	  xmax = round(10*coords.at(0))*0.1;
	  ymin = round(coords.at(1));
	  ymax = ymin;
	}
	
	_tree->traverse(xmin,ymin,xmax,ymax);
	neighbors = _tree->Search_Results;

	if (neighbors.size() > 3){
	  double sumx = 0;
	  double sumy = 0;
	  double sumz = 0;
	  for (int i = 0; i < neighbors.size(); i += 3){
	    sumx += neighbors.at(i);
	    sumy += neighbors.at(i+1);
	    sumz += neighbors.at(i+2);
	  }
	  double averagex = 3*sumx/neighbors.size();
	  double averagey = 3*sumy/neighbors.size();
	  double averagez = 3*sumz/neighbors.size();
	  interpolated_x = averagex;
	  interpolated_y = averagey;
	  interpolated_z = averagez;
	}
	
	else {
	  double x = *(_boundary_points);
	  double y = *(_boundary_points+1);
	  double z = *(_boundary_points+2);
	  double distance = sqrt(pow(x-coords.at(0),2)+pow(y-coords.at(1),2));
	  double min_distance = distance;
	  interpolated_z = z;
	  interpolated_x = x;
	  interpolated_y = y;
	  for (int i = 3; i < _num_boundary_points; i+=3){
	    x = *(_boundary_points+i);
	    y = *(_boundary_points+i+1);
	    z = *(_boundary_points+i+2);
	    distance = sqrt(pow(x-coords.at(0),2)+pow(y-coords.at(1),2));
	    if (distance < min_distance){
	      min_distance = distance;
	      interpolated_x = x;
	      interpolated_y = y;
	      interpolated_z = z;
	    }
	  }
	}
      }

      std::vector<double> results = _tree->forward_transform_for_bottom(interpolated_x,interpolated_y,interpolated_z,_height_at_origin);
      DATA_TYPE interpolated_depth = results.at(0);

      double p1[3];
      double p2[3];
      double p3[3];

      if (_tree->Search_Results.size() < 9) {
      	if (interpolated_x <= 3999){
      	  xmin = interpolated_x-0.5;
      	  xmax = interpolated_x+1;
      	  ymin = interpolated_y-4;
      	  ymax = interpolated_y+4;
      	}
	
      	else {
      	  xmin = interpolated_x-1.5;
      	  xmax = interpolated_x;
      	  ymin = interpolated_y-4;
      	  ymax = interpolated_y+4;
      	}
	
      	_tree->traverse(xmin,ymin,xmax,ymax);

      }

      if (_tree->Search_Results.size() == 9){
      	for (int i = 0; i < 3; i++){
      	  p1[i] = _tree->Search_Results.at(i);
      	  p2[i] = _tree->Search_Results.at(i+3);
      	  p3[i] = _tree->Search_Results.at(i+6);
      	}
      	double v[3];
      	double w[3];
      	for (int i = 0; i < 3; i++){
      	  v[i] = p2[i]-p1[i];
      	  w[i] = p3[i]-p1[i];
      	}
      	double N[3];
      	N[0] = v[1]*w[2]-v[2]*w[1];
      	N[1] = v[2]*w[0]-v[0]*w[2];
      	N[2] = v[0]*w[1]-v[1]*w[0];
      	std::vector<double> normal = _tree->forward_transform_for_bottom(N[0],N[1],N[2],-200.0);
      	normal.at(1) = (4*atan(1)/180)*(90-normal.at(1));
      	normal.at(2) = (4*atan(1)/180)*normal.at(2);
	normal.at(1) = normal.at(1)/normal.at(0);
	normal.at(2) = normal.at(2)/normal.at(0);
	normal.at(0) = normal.at(0)/normal.at(0);
      	if (derivative){
      	  derivative[0] = normal.at(0);
      	  derivative[1] = normal.at(1);
      	  derivative[2] = normal.at(2);
      	}
      }
      
      else if (_tree->Search_Results.size() > 9){
      	std::vector<double> Search_Results = _tree->Search_Results;
      	std::vector<double> three_closest;
      	while (three_closest.size() < 9){
      	  double x = Search_Results.at(0);
      	  double y = Search_Results.at(1);
      	  double distance = sqrt(pow(x-interpolated_x,2)+pow(y-interpolated_y,2));
      	  double min_distance = distance;
      	  double point_index = 0;
      	  for (int i = 3; i < Search_Results.size(); i += 3){
      	    x = _tree->Search_Results.at(i);
      	    y = _tree->Search_Results.at(i+1);
      	    distance = sqrt(pow(x-interpolated_x,2)+pow(interpolated_y,2));
      	    if (distance < min_distance){
      	      min_distance = distance;
      	      point_index = i;
      	    }
      	  }
      	  three_closest.push_back(Search_Results.at(point_index));
      	  three_closest.push_back(Search_Results.at(point_index+1));
      	  three_closest.push_back(Search_Results.at(point_index+2));
      	  Search_Results.erase(Search_Results.begin()+(point_index-1),Search_Results.begin()+(point_index-1)+3);
      	}
      	for (int i = 0; i < 3; i++){
      	  p1[i] = three_closest.at(i);
      	  p2[i] = three_closest.at(i+3);
      	  p3[i] = three_closest.at(i+6);
      	}
      	double v[3];
      	double w[3];
      	for (int i = 0; i < 3; i++){
      	  v[i] = p2[i]-p1[i];
      	  w[i] = p3[i]-p1[i];
      	}
      	double N[3];
      	N[0] = v[1]*w[2]-v[2]*w[1];
      	N[1] = v[2]*w[0]-v[0]*w[2];
      	N[2] = v[0]*w[1]-v[1]*w[0];
      	std::vector<double> normal = _tree->forward_transform_for_bottom(N[0],N[1],N[2],-200.0);
      	normal.at(1) = (4*atan(1)/180)*(90-normal.at(1));
      	normal.at(2) = (4*atan(1)/180)*normal.at(2);
	normal.at(1) = normal.at(1)/normal.at(0);
	normal.at(2) = normal.at(2)/normal.at(0);
	normal.at(0) = normal.at(0)/normal.at(0);
      	if (derivative) {
      	  derivative[0] = normal.at(0);
      	  derivative[1] = normal.at(1);
      	  derivative[2] = normal.at(2);
      	}
      }
      return interpolated_depth;
    }


public:

     data_grid_quad(quadtree* tree, double origin_height, int num_boundary_points, double* boundary_points)
    {
	_tree = tree;
	_height_at_origin = origin_height;
	_num_boundary_points = num_boundary_points;
	_boundary_points = boundary_points;
    }

    /**
     * Destroys memory area for field data.
     */
    ~data_grid_quad()
    {
        /* for (unsigned n = 0; n < NUM_DIMS; ++n) { */
        /* 	if (_axis[n] != NULL) { */
        /* 		delete _axis[n]; */
        /* 	} */
        /* } */
        /* if (_data != NULL) { */
        /* 	delete[] _data; */
        /* } */
    }

}; // end data_grid_quad class
} // end of namespace types
} // end of namespace usml

#endif
