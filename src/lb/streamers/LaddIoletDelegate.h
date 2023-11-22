
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_LADDIOLETDELEGATE_H
#define HEMELB_LB_STREAMERS_LADDIOLETDELEGATE_H

#include "lb/streamers/SimpleBounceBackDelegate.h"

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      template<typename CollisionImpl>
      class LaddIoletDelegate : public SimpleBounceBackDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          LaddIoletDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
              SimpleBounceBackDelegate<CollisionType>(delegatorCollider, initParams),
                  bValues(initParams.boundaryObject)
          {
            // Adapted from YangPressureDelegate
					  for (int i = 0; i < bValues->GetLocalIoletCount(); ++i)
					  {
						  SortDirectionsCloseToIOletNormal(bValues->GetLocalIolet(i));
					  }
          }

          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latticeData,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& ii)
          {
            // Translating from Ladd, J. Fluid Mech. "Numerical simulations
            // of particulate suspensions via a discretized Boltzmann
            // equation. Part 1. Theoretical foundation", 1994
            // Eq (3.2) -- simple bounce-back -- becomes:
            //   f_i'(r, t+1) = f_i(r, t*)
            // Eq (3.3) --- modified BB -- becomes:
            //   f_i'(r, t+1) = f_i(r, t*) - 2 a1_i \rho u . c_i
            // where u is the velocity of the boundary half way along the
            // link and a1_i = w_1 / cs2

            if (!site.HasIolet(ii)) return;

            int boundaryId = site.GetIoletId();
            iolets::InOutLetVelocity* iolet =
                dynamic_cast<iolets::InOutLetVelocity*>(bValues->GetIolets()[boundaryId]);
            LatticePosition sitePos(site.GetGlobalSiteCoords());
            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[ii];

            // Couple with an external system if there is
            if (unstreamed == ii)
            {
						  iolet->DoPreStreamCoupling(site.GetIndex(), bValues->GetTimeStep(), sitePos,
							  							           hydroVars.density, hydroVars.velocity);
            }

            LatticePosition halfWay(sitePos);
            halfWay.x += 0.5 * LatticeType::CX[ii];
            halfWay.y += 0.5 * LatticeType::CY[ii];
            halfWay.z += 0.5 * LatticeType::CZ[ii];

            LatticeVelocity wallMom(iolet->GetVelocity(halfWay, bValues->GetTimeStep()));
            //TODO: Add site.GetGlobalSiteCoords() as a first argument?

            if (CollisionType::CKernel::LatticeType::IsLatticeCompressible())
            {
              wallMom *= hydroVars.density;
            }

            distribn_t correction = 2. * LatticeType::EQMWEIGHTS[ii]
                * (wallMom.x * LatticeType::CX[ii] + wallMom.y * LatticeType::CY[ii]
                    + wallMom.z * LatticeType::CZ[ii]) / Cs2;

            * (latticeData->GetFNew(SimpleBounceBackDelegate<CollisionImpl>::GetBBIndex(site.GetIndex(),
                                                                                        ii))) =
                hydroVars.GetFPostCollision()[ii] - correction;
          }

          inline void PostStepLink(geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							const Direction& ii)
					{
            if (!site.HasIolet(ii)) return;

						int boundaryId = site.GetIoletId();
						iolets::InOutLet* localIOlet = bValues->GetIolets()[boundaryId];
            Direction unstreamed = LatticeType::INVERSEDIRECTIONS[ii];

						// Finalise the coupling with the external system
            if (unstreamed == ii)
            {
						  localIOlet->DoPostStreamCoupling(site.GetIndex(), bValues->GetTimeStep(), site.GetGlobalSiteCoords());
            }
					}

        protected:
          // Copied from YangPressureDelegate
					void SortDirectionsCloseToIOletNormal(iolets::InOutLet* localIOlet)
					{
						const LatticePosition& ioletNormal = localIOlet->GetNormal();
						std::array<Direction, LatticeType::NUMVECTORS> dirs;
        		std::array<LatticeDistance, LatticeType::NUMVECTORS> dist;

	        	for (Direction k = 0; k < LatticeType::NUMVECTORS; ++k)
						{
							const LatticePosition ck = LatticePosition(LatticeType::CXD[k],
					          							                       LatticeType::CYD[k],
		                													           LatticeType::CZD[k]);
							const LatticePosition unitCk = ck.GetNormalised();
							dist[k] = LatticePosition(unitCk - ioletNormal).GetMagnitudeSquared();
        			dirs[k] = k;
      			}

						// Sort dirs by comparing any two elements of dist.
      			std::sort(dirs.begin(), dirs.end(), [&dist](Direction i, Direction j) {return dist[i] < dist[j];});

						// Store the results in the iolet object.
						localIOlet->SetDirectionsCloseToNormal(dirs.begin(), dirs.end());
      		}

        private:
          iolets::BoundaryValues* bValues;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_LADDIOLETDELEGATE_H */
