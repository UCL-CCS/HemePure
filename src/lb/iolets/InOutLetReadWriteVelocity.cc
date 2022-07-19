
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/InOutLetReadWriteVelocity.h"
#include "lb/iolets/BoundaryComms.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "util/utilityFunctions.h"
#include "util/utilityStructs.h"
#include "configuration/SimConfig.h"
#include <sys/stat.h>
#include <fstream>
#include <algorithm>
#include <cmath>

namespace hemelb
{
	namespace lb
	{
		namespace iolets
		{
			InOutLetReadWriteVelocity::InOutLetReadWriteVelocity() :
				InOutLetVelocity(), units(NULL), couplingTimeStep(1)
			{
			}

			InOutLet* InOutLetReadWriteVelocity::Clone() const
			{
				InOutLet* copy = new InOutLetReadWriteVelocity(*this);
				return copy;
			}

			void InOutLetReadWriteVelocity::DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep)
      		{
        		if (comms->GetNumProcs() == 1) return;

        		comms->Receive(&maxVelocity);
        		comms->Send(&maxVelocityNew);
        		comms->WaitAllComms();
      		}

			void InOutLetReadWriteVelocity::DoPreStreamCoupling(const site_t& siteID,
																const LatticeTimeStep& timeStep,
                            		            				const LatticeVector& sitePos,
                    		                    				const LatticeDensity& density,
            		                            				const LatticeVelocity& velocity)
      		{
        		if (siteID == centreSiteID && couplingTimeStep == timeStep)
				{
					struct stat infile, outfile;
					bool updated = false, written = false;
					do
					{
						if(stat(velocityFilePath.c_str(), &infile) == 0)
						{
							std::fstream infile(velocityFilePath.c_str(), std::ios_base::in);
							log::Logger::Log<log::Debug, log::OnePerCore>("Reading iolet values from file:");

							double timeRead, valueRead;
							infile >> timeRead >> valueRead;
							infile.close();
							log::Logger::Log<log::Trace, log::OnePerCore>("Time: %f Value: %f", timeRead, valueRead);

							if (timeRead == couplingTimeStep)
							{
								couplingTimeStep += couplingFrequency;
								std::remove(velocityFilePath.c_str());
								maxVelocityNew = units->ConvertVelocityToLatticeUnits(valueRead * velocityConversionFactor);
								updated = true;
								printf("siteID %ld, timeStep %ld, couplingTimeStep %ld, maxV %lf, valueRead %lf\n", siteID, timeStep, couplingTimeStep, maxVelocityNew, valueRead);
							}
						}

						if(stat(pressureFilePath.c_str(), &outfile) != 0)
						{
							std::fstream outfile(pressureFilePath.c_str(), std::ios_base::out);
							log::Logger::Log<log::Debug, log::OnePerCore>("Writing iolet values to file:");

							LatticeTimeStep timeWrite = timeStep + couplingFrequency;
							PhysicalPressure valueWrite;
							valueWrite = units->ConvertPressureToPhysicalUnits(density * Cs2) * pressureConversionFactor;
							log::Logger::Log<log::Trace, log::OnePerCore>("Time: %f Value: %f", timeWrite, valueWrite);

							outfile << timeWrite << " " << valueWrite;
							outfile.close();
							written = true;
							printf("siteID %ld, timeStep %ld, couplingTimeStep %ld, press %lf, valueWrite %lf\n", siteID, timeStep, couplingTimeStep, density*Cs2, valueWrite);
						}
					}
					while (!(updated & written));
				}
      		}

			void InOutLetReadWriteVelocity::DoPostStreamCoupling(const site_t& siteID,
																 const LatticeTimeStep& timeStep,
                            		            				 const LatticeVector& sitePos)
      		{
        		if (siteID == centreSiteID)
        		{
          			maxVelocity = maxVelocityNew;
        		}
      		}

			LatticeVelocity InOutLetReadWriteVelocity::GetVelocity(const LatticePosition& x,
					const LatticeTimeStep t) const
			{

				if (!useWeightsFromFile)
				{
					// v(r) = vMax (1 - r**2 / a**2)
					// where r is the distance from the centreline
					LatticePosition displ = x - position;
					LatticeDistance z = displ.Dot(normal);
					Dimensionless rSqOverASq = (displ.GetMagnitudeSquared() - z * z) / (radius * radius);
					assert(rSqOverASq <= 1.0);

					// brackets to ensure that the scalar multiplies are done before vector * scalar
					return normal * (maxVelocity * (1. - rSqOverASq));
				} else {
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

					LatticeVelocity v_tot = 0;
					int iterations = 0;

					while (iterations < 3)
					{
						if (weights_table.count(xyz) > 0)
						{
							v_tot = normal * weights_table.at(xyz) * maxVelocity;
							return v_tot;
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
					return normal * 0.0;
				}

			}

			void InOutLetReadWriteVelocity::Initialise(const util::UnitConverter* unitConverter)
			{
				log::Logger::Log<log::Warning, log::Singleton>(" --> initialising vInlet");
				units = unitConverter;

				useWeightsFromFile = false;
#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
				useWeightsFromFile = true;
#endif

				if (useWeightsFromFile) {
					// If the new velocity approximation is enabled, then we want to create a lookup table here.
					const std::string in_name = velocityFilePath + ".weights.txt";
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
