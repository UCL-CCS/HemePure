
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#include "lb/iolets/InOutLetFileWK.h"
#include "lb/iolets/BoundaryComms.h"
#include "configuration/SimConfig.h"
#include <fstream>
#include "log/Logger.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "util/utilityStructs.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {
      InOutLetFileWK::InOutLetFileWK() :
          InOutLetWK(), area(1.0), units(NULL)
      {
      }

      InOutLet* InOutLetFileWK::Clone() const
      {
        InOutLetFileWK* copy = new InOutLetFileWK(*this);
        return copy;
      }

      InOutLetFileWK::~InOutLetFileWK()
      {
      }
	      
      distribn_t InOutLetFileWK::GetScaleFactor(const LatticePosition& x) const
      {
	/* These absolute normal values can still be negative here,
	 * but are corrected below to become positive. */
	double abs_normal[3] = {normal.x, normal.y, normal.z};

	// prevent division by 0 errors if the normals are 0.0
	if (normal.x < 0.0000001) { abs_normal[0] = 0.0000001; }
	if (normal.y < 0.0000001) { abs_normal[1] = 0.0000001; }
	if (normal.z < 0.0000001) { abs_normal[2] = 0.0000001; }

	int xyz_directions[3] = { 1, 1, 1 };

	std::vector<int> xyz;
	xyz.push_back(0);
	xyz.push_back(0);
	xyz.push_back(0);

	double xyz_residual[3] = {0.0, 0.0, 0.0};
	/* The residual values increase by the normal values at every time step. When they hit >1.0, then
	 * xyz is incremented and a new grid point is attempted.
	 * In addition, the specific residual value is decreased by 1.0. */

	if (normal.x < 0.0)
	{
		xyz_directions[0] = -1;
		xyz[0] = floor(x.x);
		abs_normal[0] = -abs_normal[0];
		// start with a negative residual because we already moved partially in this direction
		xyz_residual[0] = -(x.x - floor(x.x));
	} else {
		xyz[0] = std::ceil(x.x);
		xyz_residual[0] = -(std::ceil(x.x) - x.x);
	}

	if (normal.y < 0.0)
	{
		xyz_directions[1] = -1;
		xyz[1] = floor(x.y);
		abs_normal[1] = -abs_normal[1];
		xyz_residual[1] = -(x.y - floor(x.y));
	} else {
		xyz[1] = std::ceil(x.y);
		xyz_residual[1] = -(std::ceil(x.y) - x.y);
	}

	if (normal.z < 0.0)
	{
		xyz_directions[2] = -1;
		xyz[2] = floor(x.z);
		abs_normal[2] = -abs_normal[2];
		xyz_residual[2] = -(x.z - floor(x.z));
	} else {
		xyz[2] = std::ceil(x.z);
		xyz_residual[2] = -(std::ceil(x.z) - x.z);
	}

	distribn_t scale_factor = 0;
	int iterations = 0;

	while (iterations < 3)
	{
		if (weights_table.count(xyz) > 0)
		{
        
			// Q = vLocal (0.5pi a**2)(a**2/(a**2 - r**2)
                        // where r is the distance from the centreline
	                // a is the radius of the circular iolet
                        scale_factor = 0.5 * area / weights_table.at(xyz);
		
			return scale_factor;
		}

		// propagate residuals to the move to the next grid point
		double xstep = (1.0 - xyz_residual[0]) / abs_normal[0];
		double ystep = (1.0 - xyz_residual[1]) / abs_normal[1];
		double zstep = (1.0 - xyz_residual[2]) / abs_normal[2];

		double all_step = 0.0;
		int xyz_change = 0;

		if(xstep < ystep) {
			if (xstep < zstep) {
				all_step = xstep;
				xyz_change = 0;
			} else {
				if (ystep < zstep) {
					all_step = ystep;
					xyz_change = 1;
				} else {
					all_step = zstep;
					xyz_change = 2;
				}
			}
		} else {
			if (ystep < zstep) {
				all_step = ystep;
				xyz_change = 1;
			} else {
				all_step = zstep;
				xyz_change = 2;
			}
		}

		xyz_residual[0] += abs_normal[0] * all_step;
		xyz_residual[1] += abs_normal[1] * all_step;
		xyz_residual[2] += abs_normal[2] * all_step;

		xyz[xyz_change] += xyz_directions[xyz_change];

		xyz_residual[xyz_change] -= 1.0;

		iterations++;
	}

	/* Lists the sites which should be in the wall, outside of the main inlet.
	 * If you are unsure, you can increase the log level of this, run HemeLB
	 * for 1 time step, and plot these points out. */
	log::Logger::Log<log::Trace, log::OnePerCore>("%f %f %f", x.x, x.y, x.z);
	return 0.0;

      }
    
    
	void InOutLetFileWK::Initialise(const util::UnitConverter* unitConverter)
	{
		log::Logger::Log<log::Warning, log::Singleton>(" --> initialising WK Outlet");
		units = unitConverter;

		useWeightsFromFile = false;
#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
		useWeightsFromFile = true;
#endif

		if (useWeightsFromFile) {
			// If the new velocity approximation is enabled, then we want to create a lookup table here.
			const std::string in_name = wkWeightsFilePath;
			util::check_file(in_name.c_str());

			// load and read file
			std::fstream myfile;
			myfile.open(in_name.c_str(), std::ios_base::in);
			log::Logger::Log<log::Warning, log::Singleton>(" ----> loading weights file: %s",in_name.c_str());

			std::string input_line;
			/* input files are in ASCII, in format:
			 * coord_x coord_y coord_z weights_value */
			while (myfile.good())
			{
				int x, y, z;
				double v;
				myfile >> x >> y >> z >> v;

				std::vector<int> xyz;
				xyz.push_back(x);
				xyz.push_back(y);
				xyz.push_back(z);
				weights_table[xyz] = v;

				log::Logger::Log<log::Trace, log::OnePerCore>("%lld %lld %lld %f",
						x,
						y,
						z,
						weights_table[xyz]);
			}
			myfile.close();
		}
	}

    
    }
  }
}
