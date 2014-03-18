/**
 * @file boundary_grid.h
 * Creates a bottom model from a 1-D or 2-D data grid.
 */

#ifndef USML_OCEAN_BOUNDARY_GRID_QUAD_H
#define USML_OCEAN_BOUNDARY_GRID_QUAD_H

#include <usml/ocean/boundary_model.h>
#include <usml/ocean/reflect_loss_rayleigh.h>

namespace usml {
namespace ocean {

/**
 * Bottom model constructed from a 1-D or 2-D data grid.
 * The coordinate system for each kind of data set is:
 *
 *      - 1-D: Assumes that the bottom depth is a function of latitude and
 *             that the geodetic axes have been transformed to
 *             their spherical earth equivalents (theta).
 *      - 2-D: Assumes that the order of axes in the grid is
 *             (latitude, longitude) and that the geodetic
 *             axes have been transformed to their spherical earth
 *             equivalents (theta,phi).
 *
 * Uses the GRID_INTERP_PCHIP interpolation in both directions
 * to reduce sudden changes in surface normal direction.  Values outside of the
 * latitude/longitude axes defined by the data grid at limited to the values
 * at the grid edge.
 */
template< class DATA_TYPE, int NUM_DIMS > class boundary_grid_quad
    : public boundary_model
{
    //**************************************************
    // height model

protected:

    /** Boundary for all locations. */
    data_grid_quad<DATA_TYPE, NUM_DIMS>* _height;

public:

    /**
     * Compute the height of the boundary and it's surface normal at
     * a series of locations.
     *
     * @param location      Location at which to compute boundary.
     * @param rho           Surface height in spherical earth coords (output).
     * @param normal        Unit normal relative to location (output).
     * @param quick_interp  Determines if you want a fast nearest or pchip interp
     */
    virtual void height(const wposition& location, matrix<double>* rho,
        wvector* normal = NULL, bool quick_interp = false) {
    }

    /**
     * Compute the height of the boundary and it's surface normal at
     * a single location.  Often used during reflection processing.
     *
     * @param location      Location at which to compute boundary.
     * @param rho           Surface height in spherical earth coords (output).
     * @param normal        Unit normal relative to location (output).
     * @param quick_interp  Determines if you want a fast nearest or pchip interp
     */
    virtual void height(const wposition1& location, double* rho,
        wvector1* normal = NULL, bool quick_interp = false) {
            if (normal) {
                double loc[3] = { location.theta(), location.phi(), location.altitude() };
                double grad[3];
                *rho = this->_height->interpolate(loc, grad);
		normal->rho(grad[0]);
		normal->theta(grad[1]);
		normal->phi(grad[2]);
            } else {
	      double loc[3] = { location.latitude(), location.longitude(), location.altitude() };
                *rho = this->_height->interpolate(loc);
            }
    }

    //**************************************************
    // initialization

    /**
     * Initialize depth and reflection loss components for a boundary.
     *
     * @param height            Bottom depth (meters) as a function of position.
     *                          Assumes control of this grid and deletes
     *                          it when the class is destroyed.
     * @param reflect_loss      Reflection loss model.  Defaults to a
     *                          Rayleigh reflection for "sand" if NULL.
     *                          The boundary_model takes over ownship of this
     *                          reference and deletes it as part of its destructor.
     */
    boundary_grid_quad(data_grid_quad<DATA_TYPE, NUM_DIMS>* height,
        reflect_loss_model* reflect_loss = NULL) :
        boundary_model(reflect_loss), _height(height) {
        if (_reflect_loss_model == NULL) {
            _reflect_loss_model = new reflect_loss_rayleigh(
                reflect_loss_rayleigh::SAND);
        }
    }

    /**
     * Delete boundary grid.
     */
    virtual ~boundary_grid_quad() {
        delete _height;
    }

};

}  // end of namespace ocean
}  // end of namespace usml

#endif
