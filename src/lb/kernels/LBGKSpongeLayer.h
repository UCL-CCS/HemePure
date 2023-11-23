
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_LBGK_SPONGELAYER_H
#define HEMELB_LB_KERNELS_LBGK_SPONGELAYER_H

#include <cstdlib>
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
	namespace lb
	{
		namespace kernels
		{
			/**
			 * LBGKSpongeLayer: This class implements the LBGK single-relaxation time kernel with a viscous sponge layer.
			 */
			template<class LatticeType>
				class LBGKSpongeLayer : public BaseKernel<LBGKSpongeLayer<LatticeType>, LatticeType>
			{
				public:
					LBGKSpongeLayer(InitParams& initParams) :
						tau0(initParams.lbmParams->GetTau()), vRatio(initParams.lbmParams->ViscosityRatio),
						state(initParams.state)
					{
						InitState(initParams);
					}

					inline void DoCalculateDensityMomentumFeq(HydroVars<LBGKSpongeLayer<LatticeType> >& hydroVars, site_t index)
					{
						LatticeType::CalculateDensityMomentumFEq(hydroVars.f,
								hydroVars.density,
								hydroVars.momentum.x,
								hydroVars.momentum.y,
								hydroVars.momentum.z,
								hydroVars.velocity.x,
								hydroVars.velocity.y,
								hydroVars.velocity.z,
								hydroVars.f_eq.f);

						for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
						{
							hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
						}

						CalculateTau(hydroVars, index);
					}

					inline void DoCalculateFeq(HydroVars<LBGKSpongeLayer>& hydroVars, site_t index)
					{
						LatticeType::CalculateFeq(hydroVars.density,
								hydroVars.momentum.x,
								hydroVars.momentum.y,
								hydroVars.momentum.z,
								hydroVars.f_eq.f);

						for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
						{
							hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
						}

						CalculateTau(hydroVars, index);
					}

					inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<LBGKSpongeLayer>& hydroVars)
					{
						distribn_t omega = - 1.0 / hydroVars.tau;

						for (Direction direction = 0; direction < LatticeType::NUMVECTORS; ++direction)
						{
							hydroVars.SetFPostCollision(direction,
									hydroVars.f[direction] + hydroVars.f_neq.f[direction] * omega);
						}
					}

				private:
          			/**
           			*  Helper method to set/update member variables. Called from the constructor and Reset()
           			*
           			*  @param initParams struct used to store variables required for initialisation of various operators
           			*/
          			void InitState(const kernels::InitParams& initParams)
          			{
            			vTau.resize(initParams.latDat->GetLocalFluidSiteCount());
						// Width of a sponge layer (in number of sites)
						const LatticeDistance width = initParams.lbmParams->SpongeLayerWidth;
						const LatticeDistance widthSq = width * width;

			            for (site_t i = 0; i < vTau.size(); i++)
            			{
							distribn_t vRatioTot = 1.0;
              				const LatticeVector& siteLocation = initParams.latDat->GiveMeGlobalSiteCoords(i);
              				for (int j = 0; j < initParams.outletPositions.size(); j++)
              				{
                				const LatticeDistance distSq = (siteLocation - initParams.outletPositions[j]).GetMagnitudeSquared();
                				if (distSq <= widthSq)
                				{
									// Quadratic function
									vRatioTot *= vRatio - ((vRatio - 1.0) / widthSq) * distSq;

									// Sinusoidal function
									//vRatioTot *= (0.5 * (vRatio - 1.0)) * (1.0 + cos((PI / widthSq) * distSq)) + 1.0;
                				}
              				}
							// Note that viscosity is proportional to (tau - 0.5)
							vTau[i] = vRatioTot * (tau0 - 0.5) + 0.5;
            			}
          			}

					/**
					* Calculate the relaxation time of the collision and temporaily store it in hydroVars.
					* The sponge layer is maintained for a certain time and then dissolved.
					*/
					inline void CalculateTau(HydroVars<LBGKSpongeLayer<LatticeType> >& hydroVars, site_t index)
					{
						LatticeTimeStep timeStep = state->GetTimeStep();
						if (timeStep <= (LatticeTimeStep)vRatio)
						{
							hydroVars.tau = vTau[index];
						}
						else if (timeStep < 2 * (LatticeTimeStep)vRatio)
						{
							// Linear decay from vTau to tau0
							hydroVars.tau = (tau0 - vTau[index]) / vRatio * timeStep + (2.0 * vTau[index] - tau0);
						}
						else
						{
							hydroVars.tau = tau0;
						}
					}

					// Normal relaxation time
					const distribn_t tau0;
					// Ratio of the maximum viscosity in the sponge layer to the normal viscosity
					const Dimensionless vRatio;
					// Pointer to the simulation state which provides the current time step.
					SimulationState* state;
          			// Vector containing the viscous relaxation time for each site in the domain.
          			std::vector<distribn_t> vTau;
			};

		}
	}
}

#endif /* HEMELB_LB_KERNELS_LBGK_SPONGELAYER_H */
