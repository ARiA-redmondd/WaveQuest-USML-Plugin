#include <iomanip>
#include <list>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <cstring>
#include <map>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

struct listofpoints{//old points
  double latitude;
  double longitude;
  double depth;
};

struct node{//old pointers
  double minX;
  double minY;
  double maxX;
  double maxY;
  bool leaf;
  std::list<listofpoints> points;
  long branch;
  double deviation;


  node* qA = NULL;
  node* qB = NULL;
  node* qC = NULL;
  node* qD = NULL;
  node* fatherBranch;

};

class quadtree {
public : 
  quadtree(double std_dev, const char* filename);
  void populate_tree();
  void insert(node* current, double minX, double minY, double maxX, double maxY, node* parent, long main_branch);
  void populatenodelist(node* current);
  void traverse_passive(double search_minx, double search_miny , double search_maxx, double search_maxy);
  void traverse(double search_minx, double search_miny , double search_maxx, double search_maxy);
  void traverse_passive_helper(node* current, double search_minx, double search_miny , double search_maxx, double search_maxy);
  void traverse_helper(node* current, double search_minx, double search_miny , double search_maxx, double search_maxy);
  bool standard_deviation(node* current);
  bool compare_ranges(double search_minx, double search_miny, double search_maxx, double search_maxy, node* current);
  void test_tree (node* current);
  void checkleaves(node* current);
  ~quadtree();
  void destroy_tree(node* current);	


  std::list<listofpoints> populate_vector(const char* filename);
  std::vector<double> ltp_to_ecef(double xt, double yt, double zt, double h);
  std::vector<double> ecef_to_ltp(double xe, double ye, double ze, double h);
  std::vector<double> ecef_to_wgs84_spolar(double x, double y, double z);
  std::vector<double> ecef_to_wgs84_spolar_for_bottom(double x, double y, double z);
  std::vector<double> wgs84_spolar_to_ecef(double h, double lat, double lng);
  std::vector<double> forward_transform(double x, double y, double z, double h);
  std::vector<double> forward_transform_for_bottom(double x, double y, double z, double height_mid_map);
  std::vector<double> inverse_transform(double latitude, double longitude, double height, double height_mid_map);

  node* root;
  double THRESHOLD;
  int numberofleaves;
  std::vector<double> Search_Results;

private : 
  long main_branch;
  double class_maxX;
  double class_maxY;
  double class_minX;
  double class_minY;
  long temp_num;
  std::list<listofpoints> data;
  long max_branch;
 
  int branch;
  int numberofpoints;
  
 
  int FINENESS_LEVEL;
  
  double midheight;
};
