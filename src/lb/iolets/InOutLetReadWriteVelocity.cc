
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
				InOutLetVelocity(), units(NULL), area(1), weights_sum(0), warmUpLength(0),
				maxVelocity(1), maxVelocityNew(1), couplingTimeStep(2), smoothingFactor(1)
			{
			}

			InOutLet* InOutLetReadWriteVelocity::Clone() const
			{
				InOutLet* copy = new InOutLetReadWriteVelocity(*this);
				return copy;
			}

			void InOutLetReadWriteVelocity::DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep)
      		{
				/**
				* Here the send and receive requests are placed. The message is received at or before the wait
				* barrier set by BoundaryValues::FinishReceive().
				*/
        		comms->Receive(&maxVelocity);
        		comms->Send(&maxVelocityNew);

				// Here the reductions are blocking communications; they have to be made non-blocking
				const BoundaryCommunicator& bcComm = comms->GetCommunicator();
				densityAvg = bcComm.Reduce(densitySum, MPI_SUM, bcComm.GetBCProcRank());
				siteCount = bcComm.Reduce(siteCount, MPI_SUM, bcComm.GetBCProcRank());
				if (siteCount != 0)
				{
					densityAvg = densityAvg / siteCount;
				}

				densitySum = 0.0;
				siteCount = 0;
      		}

			void InOutLetReadWriteVelocity::DoPreStreamCoupling(const site_t& siteID,
																const LatticeTimeStep& timeStep,
                            		            				const LatticeVector& sitePos,
                    		                    				const LatticeDensity& density,
            		                            				const LatticeVelocity& velocity)
      		{
				if (siteID == centreSiteID && timeStep == warmUpLength + couplingTimeStep)
				{
					double couplingTime = (double)startTime + (double)(couplingTimeStep - 2) * units->GetTimeStepLength();
					printf("Looking for a flow rate value at time %e\n", couplingTime);
					bool updated = false;
					while (!updated)
					{
						struct stat infile;
						if (stat(flowRateFilePath.c_str(), &infile) == 0)
						{
							std::fstream infile(flowRateFilePath.c_str(), std::ios_base::in);
							log::Logger::Log<log::Debug, log::OnePerCore>("Reading flow rate value from file: %s", flowRateFilePath.c_str());

							double timeRead, valueRead;
							infile >> std::scientific >> timeRead >> std::scientific >> valueRead;
							infile.close();
							log::Logger::Log<log::Debug, log::OnePerCore>("timeRead: %e, valueRead: %e", timeRead, valueRead);

							//printf("timeRead %e, couplingTime %e, diff %e\n", timeRead, couplingTime, std::abs(timeRead - couplingTime));
							if (std::abs(timeRead - couplingTime) < 0.5 * units->GetTimeStepLength())
							{
								//std::remove(flowRateFilePath.c_str());
								maxVelocityNew = units->ConvertVelocityToLatticeUnits(ConvertFlowRateToVelocity(valueRead * flowRateConversionFactor));
								updated = true;
								printf("timeStep %lu, timeRead %e, valueRead %e, maxV_phy %.15lf, maxV %.15lf\n", timeStep, timeRead, valueRead, ConvertFlowRateToVelocity(valueRead * flowRateConversionFactor), maxVelocityNew);
							}
							// Enhance stability by applying a low-pass filter, which basically implements an exponential moving average.
							maxVelocityNew = smoothingFactor * maxVelocityNew + (1.0 - smoothingFactor) * maxVelocity;
						}
						else
						{
							printf("Looking for a flow rate file at %s\n", flowRateFilePath.c_str());
						}
					};

					bool written = false;
					while (!written)
					{
						struct stat outfile;
						if (stat(pressureFilePath.c_str(), &outfile) == 0)
						{
							std::fstream outfile(pressureFilePath.c_str(), std::ios_base::out);
							log::Logger::Log<log::Debug, log::OnePerCore>("Writing pressure value to file: %s", pressureFilePath.c_str());

							double timeWrite = (double)startTime + (double)(couplingTimeStep + couplingFrequency - 2) * units->GetTimeStepLength();
							double valueWrite = units->ConvertPressureToPhysicalUnits(densityAvg * Cs2) * pressureConversionFactor;
							log::Logger::Log<log::Debug, log::OnePerCore>("timeWrite: %e, valueWrite: %e", timeWrite, valueWrite);

							outfile << std::scientific << timeWrite << " " << std::scientific << valueWrite;
							outfile.close();
							written = true;
							printf("timeStep %lu, timeWrite %e, densityAvg %.15lf, valueWrite %e\n", timeStep, timeWrite, densityAvg, valueWrite);
						}
						else
						{
							printf("Looking for a pressure file at %s\n", pressureFilePath.c_str());
						}
					};

					couplingTimeStep += couplingFrequency;
				}
				densitySum += density;
				siteCount ++;
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
				// Get the max velocity
				LatticeSpeed max = maxVelocity;

				// If we're in the warm-up phase, scale down the imposed velocity
				if (t <= warmUpLength)
				{
					max *= t / LatticeSpeed(warmUpLength + 1);
				}

				if (!useWeightsFromFile)
				{
					// v(r) = vMax (1 - r**2 / a**2)
					// where r is the distance from the centreline
					LatticePosition displ = x - position;
					LatticeDistance z = displ.Dot(normal);
					LatticeDistance rSq = displ.GetMagnitudeSquared() - z * z;
					Dimensionless rSqOverASq = rSq / (radius * radius);
					if (rSqOverASq > 1.0)
					{
						log::Logger::Log<log::Error, log::OnePerCore>(
							"An IOLET site with r = %lf lies outside the IOLET radius %lf.",
							std::sqrt(rSq), radius);
						std::exit(16);
					}

					// brackets to ensure that the scalar multiplies are done before vector * scalar
					return normal * (max * (1. - rSqOverASq));
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
							v_tot = normal * weights_table.at(xyz) * max;
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
					const std::string in_name = weightsFilePath;
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
						weights_sum += v;

						log::Logger::Log<log::Trace, log::OnePerCore>("%lld %lld %lld %f",
								x,
								y,
								z,
								weights_table[xyz]);
					}
					weights_sum = units->ConvertAreaToPhysicalUnits(weights_sum);
					myfile.close();
				}
				else
				{
					weights_sum = 0.5 * area;
				}

				struct stat infile;
				if (stat(flowRateFilePath.c_str(), &infile) == 0)
				{
					std::fstream infile(flowRateFilePath.c_str(), std::ios_base::in);
					log::Logger::Log<log::Debug, log::Singleton>("Reading flow rate value from file: %s", flowRateFilePath.c_str());

					double timeRead, valueRead;
					infile >> std::scientific >> timeRead >> std::scientific >> valueRead;
					infile.close();
					log::Logger::Log<log::Debug, log::Singleton>("timeRead: %e, valueRead: %e", timeRead, valueRead);

					startTime = timeRead;
					maxVelocity = units->ConvertVelocityToLatticeUnits(ConvertFlowRateToVelocity(valueRead * flowRateConversionFactor));
					maxVelocityNew = maxVelocity;
					printf("Initialisation: timeRead %e, valueRead %e, maxV_phy %.15lf, maxV %.15lf\n", timeRead, valueRead, ConvertFlowRateToVelocity(valueRead * flowRateConversionFactor), maxVelocityNew);
				}
				else
				{
					log::Logger::Log<log::Error, log::Singleton>("Missing flow rate file at %s", flowRateFilePath.c_str());
					std::exit(17);
				}
			}

		}
	}
}
