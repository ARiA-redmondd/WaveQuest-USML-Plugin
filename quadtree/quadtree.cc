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
#include "quadtree.h"

quadtree::quadtree(double std_dev, const char* filename){
   root = new node;
   
   root->branch = 0;
   
   main_branch = 0;
   class_minX = 9999.99;
   class_minY = 9999.99;
   class_maxX = -9999.99;
   class_maxY = -9999.99;
   temp_num = 0;
   data = populate_vector(filename);

   root->minX = class_minX;
   root->maxX = class_maxX;
   root->minY = class_minY;
   root->maxY = class_maxY;


   max_branch = 0;
   branch = 0;
   numberofpoints = 0;
   numberofleaves = 0;
   THRESHOLD = std_dev;
   FINENESS_LEVEL = 0;
  }

  void quadtree::populate_tree(){

    double class_midX = (root->maxX+root->minX)/2.0;
    double class_midY = (root->maxY+root->minY)/2.0;
    root->branch = main_branch;
    root->qA = new node;
    root->qB = new node;
    root->qC = new node;
    root->qD = new node;

    while (!data.empty()){
      root->points.push_back(data.front());
      data.erase(data.begin());
    }
    insert(root->qA, root->minX, root->minY, class_midX, class_midY, root, main_branch);
    insert(root->qB, class_midX, root->minY, root->maxX, class_midY, root, main_branch);
    insert(root->qC, root->minX, class_midY, class_midX, root->maxY, root, main_branch);
    insert(root->qD, class_midX, class_midY, root->maxX, root->maxY, root, main_branch);

  }
  
   void quadtree::insert(node* current, double minX, double minY, double maxX, double maxY, node* parent, long main_branch){
    //current = new node;
   if(!parent->points.empty()){
    current->qA = new node;
    current->qB = new node;
    current->qC = new node;
    current->qD = new node;

    current->minX = minX;
    current->minY = minY;
    current->maxX = maxX;
    current->maxY = maxY;
    current->branch = 0;
    current->fatherBranch = parent;
    current->branch = main_branch;

    if(main_branch > max_branch){
      max_branch = main_branch;
    }

    populatenodelist(current);
    double midX = (current->maxX+current->minX)/2.0;
    double midY = (current->maxY+current->minY)/2.0;
    


    if(standard_deviation(current) && !current->points.empty()){
      current->leaf = true;
      delete current->qA;
      delete current->qB;
      delete current->qC;
      delete current->qD;
      current->qA = NULL;
      current->qB = NULL;
      current->qC = NULL;
      current->qD = NULL;
      numberofleaves++;
    }
    else {
      insert(current->qA, minX, minY, midX, midY, current, main_branch + 1 );
      insert(current->qB, midX, minY, maxX, midY, current, main_branch + 1 );
      insert(current->qC, minX, midY, midX, maxY, current, main_branch + 1 );
      insert(current->qD, midX, midY, maxX, maxY, current, main_branch + 1 );
    }
  }
 }
void quadtree::populatenodelist(node* current) {
    std::list<listofpoints> temp;
    double maX = root->minX;
    double miX = root->maxX;
    double maY = root->minY;
    double miY = root->maxY;
		
    
    while(!current->fatherBranch->points.empty()){
      double currentlat = current->fatherBranch->points.front().latitude;
      double currentlon = current->fatherBranch->points.front().longitude;
      
      if(currentlat <= current->maxX && currentlat >= current->minX && currentlon <= current->maxY && currentlon >= current->minY){
	current->points.push_back(current->fatherBranch->points.front());
	if(maX <= currentlat){
	  maX = currentlat;
	}
	if(miX >= currentlat){
	  miX = currentlat;
	}
	
	if(maY <= currentlon){
	  maY = currentlon;
	}
	if(miY >= currentlon){
	  miY = currentlon;
	}
      }
      else{
	temp.push_back(current->fatherBranch->points.front());
      }
      current->fatherBranch->points.pop_front();
    }


    current->fatherBranch->points = temp;
    temp.clear();
    current->maxX = maX;
    current->minX = miX;
    current->maxY = maY;
    current->minY = miY;
  }

  void quadtree::traverse_passive_helper(node* current, double search_minx, double search_miny , double search_maxx, double search_maxy){
    if(compare_ranges(search_minx, search_miny, search_maxx, search_maxy, current) && !current->points.empty()){
	Search_Results.push_back(current->points.front().latitude); //latitude
	Search_Results.push_back(current->points.front().longitude); //longitude
	Search_Results.push_back(current->points.front().depth); //depth
    }
    
    if(current->qA != NULL && (!((search_minx < current->qA->minX && search_maxx < current->qA->minX) || (search_miny < current->qA->minY && search_maxy < current->qA->minY))))
      {
	traverse_passive_helper(current->qA, search_minx, search_miny, search_maxx, search_maxy);
      }
    if(current->qB != NULL && (!((search_minx < current->qB->minX && search_maxx < current->qB->minX) || (search_miny < current->qB->minY && search_maxy < current->qB->minY))))
      {
		traverse_passive_helper(current->qB, search_minx, search_miny, search_maxx, search_maxy);
      }
    if(current->qC != NULL && (!((search_minx < current->qC->minX && search_maxx < current->qC->minX) || (search_miny < current->qC->minY && search_maxy < current->qC->minY))))
      {
	traverse_passive_helper(current->qC, search_minx, search_miny, search_maxx, search_maxy);
      }
    if(current->qD != NULL && (!((search_minx < current->qD->minX && search_maxx < current->qD->minX) || (search_miny < current->qD->minY && search_maxy < current->qD->minY))))
      {
	traverse_passive_helper(current->qD, search_minx, search_miny, search_maxx, search_maxy);
      }
  }

void quadtree::traverse_passive(double search_minx, double search_miny , double search_maxx, double search_maxy){
  if(!Search_Results.empty()){
  Search_Results.erase(Search_Results.begin(),Search_Results.end());
  }
  traverse_passive_helper(root,search_minx,search_miny ,search_maxx, search_maxy);
}

void quadtree::traverse(double search_minx, double search_miny , double search_maxx, double search_maxy){
  Search_Results.erase(Search_Results.begin(),Search_Results.end());
  traverse_helper(root,search_minx, search_miny , search_maxx, search_maxy);
}

  void quadtree::traverse_helper(node* current, double search_minx, double search_miny , double search_maxx, double search_maxy){
    if(compare_ranges(search_minx, search_miny, search_maxx, search_maxy, current) && !current->points.empty()){
      std::list<listofpoints> temp_list = current->points;
      //      for(int i = 0; i< current->points.size(); i++){
      while (!temp_list.empty()){
	Search_Results.push_back(temp_list.front().latitude); //latitude
	Search_Results.push_back(temp_list.front().longitude); //longitude
	Search_Results.push_back(temp_list.front().depth); //depth
	// Search_Results.push_back(current->points.at(i).latitude); //latitude
	// Search_Results.push_back(current->points.at(i).longitude); //longitude
	// Search_Results.push_back(current->points.at(i).depth); //depth
	temp_list.pop_front();
	temp_num++;
      }
      temp_list.clear();
    }
    
    
    if(current->qA != NULL && (!((search_minx < current->qA->minX && search_maxx < current->qA->minX) || (search_miny < current->qA->minY && search_maxy < current->qA->minY))))
      {
	traverse_helper(current->qA, search_minx, search_miny, search_maxx, search_maxy);
      }
    if(current->qB != NULL && (!((search_minx < current->qB->minX && search_maxx < current->qB->minX) || (search_miny < current->qB->minY && search_maxy < current->qB->minY))))
      {
	traverse_helper(current->qB, search_minx, search_miny, search_maxx, search_maxy);
      }
    if(current->qC != NULL && (!((search_minx < current->qC->minX && search_maxx < current->qC->minX) || (search_miny < current->qC->minY && search_maxy < current->qC->minY))))
      {
	traverse_helper(current->qC, search_minx, search_miny, search_maxx, search_maxy);
      }
    if(current->qD != NULL && (!((search_minx < current->qD->minX && search_maxx < current->qD->minX) || (search_miny < current->qD->minY && search_maxy < current->qD->minY))))
      {
	traverse_helper(current->qD, search_minx, search_miny, search_maxx, search_maxy);
      }
  }

bool quadtree::compare_ranges(double search_minx, double search_miny, double search_maxx, double search_maxy, node* current) {
		double current_minx = current->minX;
		double current_maxx = current->maxX;
		double current_miny = current->minY;
		double current_maxy = current->maxY;

		if(search_minx <= current_minx && search_maxx >= current_maxx && search_miny <= current_miny && search_maxy >= current_maxy){
			return true;
		}
		else if(search_maxx <= current_maxx && search_minx >= current_minx && search_miny >= current_miny && search_maxy <= current_maxy){
			return true;
		}
		else if(search_maxx <= current_maxx && search_maxy <= current_maxy && (search_maxx >= current_minx || search_maxy >= current_miny)){
			return true;
		}
		else if((search_minx >= current_minx || search_miny >= current_miny) && search_minx <= current_maxx && search_miny <= current_maxy){
			return true;
		}
		else{
			return false;
		}
}

  bool quadtree::standard_deviation(node* current){
    double temp = 0;
    double mean = 0.0;
    double num = 0;
    
    if(current->points.size() > 1){
      
      // for (int i = 0; i < current->points.size(); i++){
      // 	mean = mean + current->points.at(i).depth;
      // }
      
      std::list<listofpoints> temp_list = current->points;

      while (!temp_list.empty()){
	mean = mean + temp_list.front().depth;
	temp_list.pop_front();
      }

      mean = mean / current->points.size();

      temp_list = current->points;

      // for (int j = 0; j< current->points.size(); j++){
      // 	num = num + pow((current->points.at(j).depth - mean), 2);
	
      // }

      while (!temp_list.empty()){
	num = num + pow((temp_list.front().depth - mean), 2);
	temp_list.pop_front();
      }
      
      temp_list.clear();

      temp = sqrt((num/(current->points.size()-1)));
      current->deviation = temp;
      
      if(temp <= THRESHOLD){
	if(FINENESS_LEVEL <= 0){
	  return true;// if it meets the threshold
	}
	else{
	  
	  double midX = (current->maxX+current->minX)/2.0;
	  double midY = (current->maxY+current->minY)/2.0;
	  
	  // insert_std(current->qA, current->minX, current->minY, midX, midY, current, 1);
	  
	  // insert_std(current->qB, midX, current->minY, current->maxX, midY, current, 1);
	  
	  // insert_std(current->qC,current->minX, midY, midX, current->maxY, current, 1);
	  
	  // insert_std(current->qD, midX, midY, current->maxX, current->maxY, current, 1);
	  insert(current->qA, current->minX, current->minY, midX, midY, current, 1);
	  
	  insert(current->qB, midX, current->minY, current->maxX, midY, current, 1);
	  
	  insert(current->qC,current->minX, midY, midX, current->maxY, current, 1);
	  
	  insert(current->qD, midX, midY, current->maxX, current->maxY, current, 1);
	  
	  
	  return true;// if it meets the threshold
	}
      }
    }
    else{
      return true;//has 1 or less elements
    }
    
    return false; // has more than 1 listofpoints and doesn't meet standard deviation
  }

  std::list<listofpoints> quadtree::populate_vector(const char* filename) {
    std::list<listofpoints> Data;
    
    std::ifstream file;
    file.open(filename);

    char trash;
    double trash1, trash2, x, y, z;

    while (!file.eof()){
      file>>trash1>>trash>>trash2>>trash>>x>>trash>>y>>trash>>z;
	if(file.eof()){
		break;
	}	 
	 if (class_maxX < x){
	    class_maxX = x;
	  }
	  if (class_maxY < y){
	    class_maxY = y;
	  }
	  if (class_minX > x){
	    class_minX = x;
	  }
	  if (class_minY > y){
	    class_minY = y;
	  }
	  listofpoints tempPoint;
	  tempPoint.latitude = x;
	  tempPoint.longitude = y;
	  tempPoint.depth = z;
	  Data.push_back(tempPoint);
    }
file.close();
    return Data;
  }
  

std::vector<double> quadtree::ltp_to_ecef(double xt, double yt, double zt, double h){
  double x0, y0, z0;
  double xe, ye, ze;
  double a = 6378137.0;
  double e = 8.1819190842622*pow(10,-2);
  double N = a/sqrt(1-pow(e,2)/4);
  x0 = 3*(h+N)/4;
  y0 = sqrt(3)*(h+N)/4;
  z0 = (h+(1-pow(e,2))*N)/2;
  xe = x0-xt/2-sqrt(3)*yt/4+3*zt/4;
  ye = y0+sqrt(3)*xt/2-yt/4+sqrt(3)*zt/4;
  ze = z0+sqrt(3)*yt/2+zt/2;
  std::vector<double> result_vec;
  result_vec.push_back(xe);
  result_vec.push_back(ye);
  result_vec.push_back(ze);
  return result_vec;
}

std::vector<double> quadtree::ecef_to_ltp(double xe, double ye, double ze, double h){
  std::vector<double> result_vec;
  double x0, y0, z0;
  double xt, yt, zt;
  double a = 6378137.0;
  double e = 8.1819190842622*pow(10,-2);
  double N = a/sqrt(1-pow(e,2)/4);
  x0 = 3*(h+N)/4;
  y0 = sqrt(3)*(h+N)/4;
  z0 = (h+(1-pow(e,2))*N)/2;
  xt = -(xe-x0)/2+sqrt(3)*(ye-y0)/2;
  yt = -sqrt(3)*(xe-x0)/4-(ye-y0)/4+sqrt(3)*(ze-z0)/2;
  zt = 3*(xe-x0)/4+sqrt(3)*(ye-y0)/4+(ze-z0)/2;
  result_vec.push_back(xt);
  result_vec.push_back(yt);
  result_vec.push_back(zt);
  return result_vec;
}

std::vector<double> quadtree::ecef_to_wgs84_spolar(double x, double y, double z){
  std::vector<double> result_vec;
  double lat,lng,h;
  double p,Rn;
  double lat_next,lat_now, h_now;
  double a = 6378137.0;
  double e = 8.1819191*pow(10,-2);
  lng = (180/(4*atan(1)))*atan2(y,x);
  p = sqrt(pow(x,2)+pow(y,2));
  lat_now = atan2(p,z);
  for (int i = 0; i < 5; i++){
    Rn = a/sqrt(1-pow(e,2)*pow(sin(lat_now),2));
    h_now = p/cos(lat_now) - Rn;
    lat_next = atan((z/p)*pow(1-pow(e,2)*(Rn/(Rn+h_now)),-1));
    lat_now = lat_next;
  }
  lat = lat_now;
  Rn = a/sqrt(1-pow(e,2)*pow(sin(lat),2));
  h = p/cos(lat)-Rn;
  lat = (180/(4*atan(1)))*lat_now;
  result_vec.push_back(h);
  result_vec.push_back(lat);
  result_vec.push_back(lng);
  return result_vec;
}

  std::vector<double> quadtree::ecef_to_wgs84_spolar_for_bottom(double x, double y, double z){
    std::vector<double> result_vec;
    double lat,lng,h;
    double p,Rn,rho;
    double lat_next,lat_now, h_now;
    double a = 6378137.0;
    double e = 8.1819191*pow(10,-2);
    lng = (180/(4*atan(1)))*atan2(y,x);
    p = sqrt(pow(x,2)+pow(y,2));
    lat_now = atan2(p,z);
    for (int i = 0; i < 5; i++){
      Rn = a/sqrt(1-pow(e,2)*pow(sin(lat_now),2));
      h_now = p/cos(lat_now) - Rn;
      lat_next = atan((z/p)*pow(1-pow(e,2)*(Rn/(Rn+h_now)),-1));
      lat_now = lat_next;
    }
    lat = lat_now;
    Rn = a/sqrt(1-pow(e,2)*pow(sin(lat),2));
    h = p/cos(lat)-Rn;
    lat = (180/(4*atan(1)))*lat_now;
    const double f = 1.0 / 298.257223563;
    const double e2 = f * (2.0 - f);
    double sinT = sin((4*atan(1)/180)*lat);
    double w = sqrt(1.0 - e2 * sinT * sinT);
    double rm = a * (1-e2) / (w*w*w) ;
    double rv = a / w ;
    double radius_curvature = sqrt( rm * rv );
    rho = h+radius_curvature;
    result_vec.push_back(rho);
    result_vec.push_back(lat);
    result_vec.push_back(lng);
    return result_vec;
  }

std::vector<double> quadtree::wgs84_spolar_to_ecef(double h, double lat, double lng){
  std::vector<double> result_vec;
  double x,y,z;
  double a = 6378137.0;
  double e = 8.1819190842622*pow(10,-2);
  double Rn = a/sqrt(1-pow(e,2)*pow(sin((4*atan(1)/180)*lat),2));
  x = (Rn+h)*cos((4*atan(1)/180)*lat)*cos((4*atan(1)/180)*lng);
  y = (Rn+h)*cos((4*atan(1)/180)*lat)*sin((4*atan(1)/180)*lng);
  z = ((1-pow(e,2))*Rn+h)*sin((4*atan(1)/180)*lat);
  result_vec.push_back(x);
  result_vec.push_back(y);
  result_vec.push_back(z);
  return result_vec;
}
  
  std::vector<double> quadtree::forward_transform(double x, double y, double z, double h){
    std::vector<double> result1;
    std::vector<double> result2;
    result1 = ltp_to_ecef(x,y,z,h);
    result2 = ecef_to_wgs84_spolar(result1.at(0),result1.at(1),result1.at(2));
    result1.erase(result1.begin(),result1.end());
    return result2;
  }
  
  std::vector<double> quadtree::forward_transform_for_bottom(double x, double y, double z, double height_mid_map){
    std::vector<double> result1;
    std::vector<double> result2;
    result1 = ltp_to_ecef(x,y,z,height_mid_map);
    result2 = ecef_to_wgs84_spolar_for_bottom(result1.at(0),result1.at(1),result1.at(2));
    result1.erase(result1.begin(),result1.end());
    return result2;
  }
  
  std::vector<double> quadtree::inverse_transform(double latitude, double longitude, double height, double height_mid_map){
    std::vector<double> result1;
    std::vector<double> result2;
    // std::cout << latitude << std::endl;
    // std::cout << longitude << std::endl;
    // std::cout << height << std::endl;
    result1 = wgs84_spolar_to_ecef(height,latitude,longitude);
    result2 = ecef_to_ltp(result1.at(0),result1.at(1),result1.at(2),height_mid_map);
    result1.erase(result1.begin(),result1.end());
    return result2;
  }

void quadtree::destroy_tree(node* current){
  if(current->qA != NULL){
    destroy_tree(current->qA);
  }
  if(current->qB != NULL){
    destroy_tree(current->qB);
  }
  if(current->qC != NULL){
    destroy_tree(current->qC);
  }
  if(current->qD != NULL){
    destroy_tree(current->qD);
  }

 delete current;

}

quadtree::~quadtree(){
  destroy_tree(root->qA);
  destroy_tree(root->qB);
  destroy_tree(root->qC);
  destroy_tree(root->qD);
  delete root;
}
void quadtree::checkleaves(node* current){
    
    if(!current->leaf && current->points.size() > 0){
	 std::cout << "error \n";
     }
    if(current->qA != NULL)
      {
	checkleaves(current->qA);
      }
    if(current->qB != NULL)
      {
	checkleaves(current->qB);
      }
    if(current->qC != NULL)
      {
	checkleaves(current->qC);
      }
    if(current->qD != NULL)
      {
	checkleaves(current->qD);
      }

}
