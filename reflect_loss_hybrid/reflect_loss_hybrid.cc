#include <math.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <usml/ocean/reflect_loss_hybrid.h>
#include <exception>

using namespace usml::ocean;
using namespace std;
using namespace boost::math;

double reflect_loss_hybrid::interpolate(double val_shallow,double val_deep,double depth_shallow, double depth_inter, double depth_deep){
  double result;
  result = val_shallow + (val_deep-val_shallow)*(depth_inter-depth_shallow)/(depth_deep-depth_shallow);
  return result;
}

std::vector<double> reflect_loss_hybrid::surface_loss(int numf, double* freq, double grazing){
  std::vector<double> loss;
  double windspeed_knots = 1.94384*windspeed;
  double rms_disp = 30.48*((1.86*0.01*pow(windspeed_knots,2))/4.02);
  double g;
  for (int i = 0; i < numf; i++){
    double result;
    if (freq[i] <= 4000.0){
      g = 2*3.14159*0.01*freq[i]*rms_disp*sin(grazing)/c;
      g = pow(g,2);
      result = -20*log10(cyl_bessel_i(0,2*g)*exp(-2*g));
    }
    else if (freq[i] >= 10000.0){
      double temp;
      temp = (1.43*pow(10,-5)/sin(grazing))*pow(windspeed,3.75)*freq[i];
      result = min(30.0,temp);
    }
    else{
      double val_shallow;
      g = 2*3.14159*0.01*freq[i]*rms_disp*sin(grazing)/c;
      g = pow(g,2);
      val_shallow = -20*log10(cyl_bessel_i(0,2*g)*exp(-2*g));
      double val_deep;
      double temp;
      temp = (1.43*pow(10,-5)/sin(grazing))*pow(windspeed,3.75)*freq[i];
      val_deep = min(30.0,temp);
      result = interpolate(val_shallow,val_deep,4000.0,freq[i],10000.0);
    }
    loss.push_back(result);
  }
  return loss;
}

void reflect_loss_hybrid::reflect_loss(
    const wposition1& location,
    const seq_vector& frequencies, double angle,
    boost::numeric::ublas::vector<double>* amplitude, \
    boost::numeric::ublas::vector<double>* phase) {
  std::vector<double> loss;
  std::vector<double> freq_vec;
  double* freq = new double[frequencies.size()];
  for (seq_vector::const_iterator iter = frequencies.begin(); iter < frequencies.end(); iter++){
    freq_vec.push_back(*iter);
  }
  copy(freq_vec.begin(),freq_vec.end(),freq);
  loss = surface_loss(frequencies.size(),freq,angle);
  for (int i = 0; i < frequencies.size(); i++){
    (*amplitude)(i) = loss.at(i);
    if (phase){
      (*phase)(i) = 4*atan(1);
    }
  }
  loss.erase(loss.begin(),loss.end());
  freq_vec.erase(freq_vec.begin(),freq_vec.end());
  delete[] freq;
}



reflect_loss_hybrid::~reflect_loss_hybrid() {}
