
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKDELEGATE_H
#define HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
	namespace lb
	{
		namespace streamers
		{
			template<typename CollisionImpl>
				class GrinbergKarniadakisWKDelegate : public BaseStreamerDelegate<CollisionImpl>
			{
				public:
					typedef CollisionImpl CollisionType;
					typedef typename CollisionType::CKernel::LatticeType LatticeType;

					GrinbergKarniadakisWKDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
						collider(delegatorCollider), iolet(*initParams.boundaryObject), nashDelegate(delegatorCollider, initParams)
				{
				}

					inline void StreamLink(const LbmParameters* lbmParams,
							geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();
						iolets::InOutLetWK* wkIolet = dynamic_cast<iolets::InOutLetWK*>(iolet.GetIolets()[boundaryId]);

						if (site.GetIndex() == wkIolet->GetCentreSiteID())
						{
							LatticePressure pressure = wkIolet->GetPressure(0); // 0 is dummy
							distribn_t R0 = wkIolet->GetResistance();
							distribn_t C0 = wkIolet->GetCapacitance();

							LatticePosition sitePos(site.GetGlobalSiteCoords());
							distribn_t scaleFactor = wkIolet->GetScaleFactor(sitePos);

							// Calculate the velocity at the ghost site, as the component normal to the iolet.
							util::Vector3D<Dimensionless> ioletNormal = wkIolet->GetNormal();
							distribn_t component = hydroVars.velocity.Dot(ioletNormal);

							// Explicit integration scheme
							//LatticePressure pressureNew = (1.0/C0)*scaleFactor*std::abs(component) + (1.0 - 1.0/(R0*C0))*pressure;

							// Semi-implicit integration scheme
							LatticePressure pressureNew = R0/(1.0 + R0*C0) * (scaleFactor*std::abs(component) + C0*pressure);

							wkIolet->SetDensityNew(pressureNew / Cs2);
						}

						nashDelegate.StreamLink(lbmParams, latticeData, site, hydroVars, direction);
					}
					
					inline void PostStepLink(geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();
						iolets::InOutLetWK* wkIolet = dynamic_cast<iolets::InOutLetWK*>(iolet.GetIolets()[boundaryId]);
						wkIolet->SetDensity(wkIolet->GetDensityNew(0)); // 0 is dummy
					}

				protected:
					CollisionType& collider;
					iolets::BoundaryValues& iolet;

				private:
					NashZerothOrderPressureDelegate<CollisionType> nashDelegate;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKDELEGATE_H
