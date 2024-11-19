
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_KERNELS_TRT_SPONGELAYER_H
#define HEMELB_LB_KERNELS_TRT_SPONGELAYER_H

#include <cstdlib>
#include "lb/HFunction.h"
#include "util/utilityFunctions.h"
#include "lb/kernels/BaseKernel.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      /**
       *  TRTSpongeLayer: This class implements a two-relaxation time kernel.
       */
      template<class LatticeType>
      class TRTSpongeLayer : public BaseKernel<TRTSpongeLayer<LatticeType>, LatticeType>
      {
          // Store the directions as pairs of opposites
          // (Note that the zero vector is it's own opposite)
          typedef std::pair<Direction, Direction> Opposites;
          typedef std::vector<Opposites> OppList;
          OppList directionPairs;
          Direction iZero;

        public:
          TRTSpongeLayer(InitParams& initParams):
            tau0(initParams.lbmParams->GetTau()), vRatio(initParams.lbmParams->ViscosityRatio),
						lifetime(initParams.lbmParams->SpongeLayerLifetime), state(initParams.state)
          {
            for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i)
            {
              Direction iBar = LatticeType::INVERSEDIRECTIONS[i];
              if (i == iBar)
                iZero = i;

              if (iBar >= i)
                directionPairs.push_back(std::make_pair(i, iBar));
            }
            InitState(initParams);
          }

          inline void DoCalculateDensityMomentumFeq(HydroVars<TRTSpongeLayer<LatticeType> >& hydroVars, site_t index)
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

          inline void DoCalculateFeq(HydroVars<TRTSpongeLayer>& hydroVars, site_t index)
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

          inline void DoCollide(const LbmParameters* const lbmParams, HydroVars<TRTSpongeLayer>& hydroVars)
          {
            // Note HemeLB defines omega = -1/ tau
            // Magic number determines the other relaxation time
            // Lambda = (tau_plus - 1/2) (tau_minus - 1/2)
            // Choose such that HWBB walls are always in the right place.
            const distribn_t Lambda = lbmParams->GetRelaxationParameter();
            const distribn_t tau_plus = hydroVars.tau;
            const distribn_t omega_plus = lbmParams->GetOmega();
            const distribn_t tau_minus = 0.5 + Lambda / (tau_plus - 0.5);
            const distribn_t omega_minus =  -1.0 / tau_minus;

            // Special case the null velocity.
            hydroVars.SetFPostCollision(iZero,
                                        hydroVars.f[iZero] + omega_plus * hydroVars.f_neq.f[iZero]);

            // Now deal with the non-zero
            for (OppList::const_iterator oppIt = directionPairs.begin();
                oppIt != directionPairs.end();
                ++oppIt)
            {
              const Direction& i = oppIt->first;
              const Direction& iBar = oppIt->second;

              distribn_t sym = 0.5 * omega_plus * (hydroVars.f_neq.f[i] + hydroVars.f_neq.f[iBar]);
              distribn_t asym = 0.5 * omega_minus * (hydroVars.f_neq.f[i] - hydroVars.f_neq.f[iBar]);
              hydroVars.SetFPostCollision(i, hydroVars.f[i] + sym + asym);
              hydroVars.SetFPostCollision(iBar, hydroVars.f[iBar] + sym - asym);
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
					inline void CalculateTau(HydroVars<TRTSpongeLayer<LatticeType> >& hydroVars, site_t index)
					{
						LatticeTimeStep timeStep = state->GetTimeStep();
						if (timeStep <= lifetime / 2)
						{
							hydroVars.tau = vTau[index];
						}
						else if (timeStep < lifetime)
						{
							// Linear decay from vTau to tau0
							hydroVars.tau = (tau0 - vTau[index]) * 2.0 / lifetime * timeStep + (2.0 * vTau[index] - tau0);
						}
						else
						{
							hydroVars.tau = tau0;
						}
					}

          // Normal symmetric relaxation time (tau_plus)
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

#endif /* HEMELB_LB_KERNELS_TRT_SPONGELAYER_H */
