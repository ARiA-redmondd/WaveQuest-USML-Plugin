#ifndef USML_OCEAN_REFLECT_LOSS_HYBRID_H
#define USML_OCEAN_REFLECT_LOSS_HYBRID_H

#include <usml/ocean/reflect_loss_model.h>
#include <vector>
#include <algorithm>
#include <fstream>

namespace usml{
  namespace ocean{
    class USML_DECLSPEC reflect_loss_hybrid : public reflect_loss_model{
    private:

	double windspeed;

	double c;

    public:

	double interpolate(double val_shallow, double val_deep, double depth_shallow, double depth_inter, double depth_deep);

	std::vector<double> surface_loss(int numf, double* freq, double grazing);

      reflect_loss_hybrid(double w_speed, double c_speed) : 
	windspeed(w_speed), c(c_speed) {}

      virtual void reflect_loss(const wposition1& location, const seq_vector& \
				frequencies, double angle, vector<double>* \
				amplitude, vector<double>* phase=NULL);

      virtual ~reflect_loss_hybrid();
    };
  }
}
#endif
