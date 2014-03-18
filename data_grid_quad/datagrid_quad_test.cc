#include <boost/test/unit_test.hpp>
#include <usml/types/types.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>

BOOST_AUTO_TEST_SUITE(datagridquad_test)

using namespace boost::unit_test;
using namespace usml::types;

BOOST_AUTO_TEST_CASE(datagridquad_test){

  cout << "=== datagridquad_test: interpolate_test ===" << endl;
  quadtree* tree = new quadtree(0.0,USML_TEST_DIR "/types/test/wedge_cartesian.txt");
  tree->populate_tree();
  
  const int MAX_CHARS_PER_LINE = 512;
  const int MAX_TOKENS_PER_LINE = 5;
  const char* DELIMITER = ",";
  int number_of_rows = 8001;
  int number_of_columns = 1;
  std::ifstream fin;
  std::vector<std::vector<double> > grid_coords;
  std::vector<double> list_of_points;

  fin.open(USML_TEST_DIR "/types/test/wedge_cartesian.txt");
  while (!fin.eof()){
    char buf[MAX_CHARS_PER_LINE];
    fin.getline(buf,MAX_CHARS_PER_LINE);
    int n = 0;
    const char* token[MAX_TOKENS_PER_LINE] = {};
    token[0] = strtok(buf,DELIMITER);
    if (token[0]){
      std::vector<double> temp_vec;
      for (n = 1; n < MAX_TOKENS_PER_LINE; n++){
  	token[n] = strtok(0,DELIMITER);
  	if (n > 1){
  	  temp_vec.push_back(atof(token[n]));
  	}
  	if (!token[n]){
  	  break;
  	}
      }
      grid_coords.push_back(temp_vec);
      temp_vec.erase(temp_vec.begin(),temp_vec.end());
    }
  }

  for (int i = 0; i < number_of_rows; i++){
    for (int j = 0; j < number_of_columns; j++){
      if (i == 0 || j == 0){
  	list_of_points.push_back(grid_coords.at(i*number_of_columns+j).at(0));
  	list_of_points.push_back(grid_coords.at(i*number_of_columns+j).at(1));
  	list_of_points.push_back(grid_coords.at(i*number_of_columns+j).at(2));
      }
      else if (i == number_of_rows-1 || j == number_of_columns-1){
  	list_of_points.push_back(grid_coords.at(i*number_of_columns+j).at(0));
  	list_of_points.push_back(grid_coords.at(i*number_of_columns+j).at(1));
  	list_of_points.push_back(grid_coords.at(i*number_of_columns+j).at(2));
      }
    }
  }

  BOOST_CHECK_EQUAL(list_of_points.size(),3*8001);
  double* boundary_points = new double[list_of_points.size()];
  std::copy(list_of_points.begin(),list_of_points.end(),boundary_points);
  grid_coords.erase(grid_coords.begin(),grid_coords.end());
  list_of_points.erase(list_of_points.begin(),list_of_points.end());
  
  data_grid_quad<double,2>* grid = new data_grid_quad<double,2>(tree,-200,3*8001,boundary_points);
  double location[3];

  cout << "Location 1: (x,y,z) = (0,0,0)" << endl;

  std::vector<double> interp_point = tree->forward_transform(0.0,0.0,0.0,-200);
  location[0] = interp_point.at(1);
  location[1] = interp_point.at(2);
  location[2] = interp_point.at(0);
  interp_point.erase(interp_point.begin(),interp_point.end());

  double interpolated_rho = grid->interpolate(location);
  BOOST_CHECK_CLOSE(interpolated_rho,6367208.77774376702,1e-3);

  cout << "Location 2: (x,y,z) = (0.75,0,0)" << endl;

  interp_point = tree->forward_transform(0.75,0.0,0.0,-200);
  location[0] = interp_point.at(1);
  location[1] = interp_point.at(2);
  location[2] = interp_point.at(0);
  interp_point.erase(interp_point.begin(),interp_point.end());
  double normal[3];

  interpolated_rho = grid->interpolate(location,normal);
  BOOST_CHECK_CLOSE(interpolated_rho,6367208.82770177815,1e-3);
  BOOST_CHECK_CLOSE(normal[0],1,1e-3);
  BOOST_CHECK_CLOSE(normal[1],1.64467e-07,1e-3);

  cout << "Location 3: (x,y,z) = (0,10,0)" << endl;

  interp_point = tree->forward_transform(0.0,10.0,0.0,-200);
  location[0] = interp_point.at(1);
  location[1] = interp_point.at(2);
  location[2] = interp_point.at(0);
  interp_point.erase(interp_point.begin(),interp_point.end());

  interpolated_rho = grid->interpolate(location);
  BOOST_CHECK_CLOSE(interpolated_rho,6367208.77774376702,1e-3);

  delete tree;
  delete grid;
  delete[] boundary_points;
}

BOOST_AUTO_TEST_SUITE_END()
