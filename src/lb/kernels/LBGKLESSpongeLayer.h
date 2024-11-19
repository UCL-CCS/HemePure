
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_LBGK_LES_SPONGELAYER_H
#define HEMELB_LB_KERNELS_LBGK_LES_SPONGELAYER_H

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
			 * LBGKLESSpongeLayer: This class implements the LBGK single-relaxation time kernel with a viscous sponge layer.
			 */
			template<class LatticeType>
				class LBGKLESSpongeLayer : public BaseKernel<LBGKLESSpongeLayer<LatticeType>, LatticeType>
			{
				public:
					LBGKLESSpongeLayer(InitParams& initParams) :
						tau0(initParams.lbmParams->GetTau()), vRatio(initParams.lbmParams->ViscosityRatio),
						lifetime(initParams.lbmParams->SpongeLayerLifetime), state(initParams.state)
					{
						InitState(initParams);
					}

					inline void DoCalculateDensityMomentumFeq(HydroVars<LBGKLESSpongeLayer<LatticeType> >& hydroVars, site_t index)
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

					inline void DoCalculateFeq(HydroVars<LBGKLESSpongeLayer>& hydroVars, site_t index)
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

					inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<LBGKLESSpongeLayer>& hydroVars)
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
								// const int dist = (siteLocation - initParams.outletPositions[j]).GetByDirection(util::Direction::Direction::X);
								// const LatticeDistance distSq = dist * dist;
								const LatticeDistance dist = std::sqrt(distSq);
                				if (distSq <= widthSq)
                				{
									// Quadratic function
									vRatioTot *= 1.0 + (vRatio - 1.0) * (dist / width - 1.0) * (dist / width - 1.0);

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
					inline void CalculateTau(HydroVars<LBGKLESSpongeLayer<LatticeType> >& hydroVars, site_t index)
					{
						LatticeTimeStep timeStep = state->GetTimeStep();
						double tau_les = compute_tau_smagorinsky(hydroVars);
						if (timeStep <= lifetime / 2)
						{
							if (vTau[index] == tau0)
							{
								hydroVars.tau = tau_les;
							}
							else
							{
								hydroVars.tau = vTau[index];
							}
						}
						else if (timeStep < lifetime)
						{
							// Linear decay from vTau to tau0
							// hydroVars.tau = (tau0 - vTau[index]) * 2.0 / lifetime * timeStep + (2.0 * vTau[index] - tau0);
							hydroVars.tau = (tau_les - vTau[index]) * 2.0 / lifetime * timeStep + (2.0 * vTau[index] - tau_les);
						}
						else
						{
							// hydroVars.tau = tau0;
							hydroVars.tau = tau_les;
						}
					}

					double compute_tau_smagorinsky(HydroVars<LBGKLESSpongeLayer<LatticeType> >& hydroVars) const
					{
						// calculate tau using smagorinsky local correction
						double dx = 1.0;
						double dt = 1.0;
						double C = dx / dt;
						double rho1 = 1.0;
						double localTau;
						double C_smag = 0.1;
						// Compute non-equilibrium values
						for (unsigned int ii = 0; ii < LatticeType::NUMVECTORS; ++ii)
						{
							hydroVars.f_neq.f[ii] = hydroVars.f[ii] - hydroVars.f_eq.f[ii];
						}
						double Q_12 = 0.0;
						// Calculate diagonal and upper diagonal of the non equilibrium stress tensor
						for (int i = 0; i < 3; ++i) {
							for (int j = 0; j < 3; ++j) {
								double qij = 0.0;
								for (int v = 0; v < LatticeType::NUMVECTORS; ++v) {
									qij += LatticeType::discreteVelocityVectors[i][v] * LatticeType::discreteVelocityVectors[j][v] * hydroVars.f_neq.f[v];
								}
								Q_12 += qij * qij;
							}
						}
						Q_12 = sqrt(Q_12);
						// eq 36 Koda 2015, csmag is smagorinsky constant here is c_smag^2 in the paper is c_smag
						localTau = 1. / 2. *
									(tau0 + sqrt((tau0 * rho1 * C)*(tau0 * rho1 * C) + 
									18.0 * 1.4142135623730950488016887242097 * rho1 * C_smag * C_smag * Q_12) / (rho1 * C));

						return localTau;
					}

					// Normal relaxation time
					const distribn_t tau0;
					// Ratio of the maximum viscosity in the sponge layer to the normal viscosity
					const Dimensionless vRatio;
					// Lifetime of the sponge layer
					const LatticeTimeStep lifetime;
					// Pointer to the simulation state which provides the current time step.
					SimulationState* state;
          			// Vector containing the viscous relaxation time for each site in the domain.
          			std::vector<distribn_t> vTau;
			};

		}
	}
}

#endif /* HEMELB_LB_KERNELS_LBGK_LES_SPONGELAYER_H */
