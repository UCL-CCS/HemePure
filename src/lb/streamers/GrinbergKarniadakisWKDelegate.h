
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
						collider(delegatorCollider), iolet(*initParams.boundaryObject)
				{
				}

					inline void StreamLink(const LbmParameters* lbmParams,
							geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();

#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
						iolets::InOutLetFileWK* wkIolet = dynamic_cast<iolets::InOutLetFileWK*>(iolet.GetIolets()[boundaryId]);
#else
						iolets::InOutLetWK* wkIolet = dynamic_cast<iolets::InOutLetWK*>(iolet.GetIolets()[boundaryId]);
#endif
						LatticeDensity ghostDensity = wkIolet->GetDensity(0.0);
						
						// Calculate the velocity at the ghost site, as the component normal to the iolet.
						util::Vector3D<float> ioletNormal = wkIolet->GetNormal();

						// Note that the division by density compensates for the fact that v_x etc have momentum
						// not velocity.
						distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal);

						if (site.GetIndex() == wkIolet->GetCentreSiteID())
						{
							LatticePressure pressure = wkIolet->GetPressure(0);
							distribn_t R0 = wkIolet->GetResistance();
							distribn_t C0 = wkIolet->GetCapacitance();

							LatticePosition sitePos(site.GetGlobalSiteCoords());
							distribn_t scaleFactor = wkIolet->GetScaleFactor(sitePos);

							// Explicit integration scheme
							//LatticeDensity ghostDensityNew = ((1.0/C0)*scaleFactor*std::abs(component) + (1.0 - 1.0/(R0*C0))*pressure)/Cs2;

							// Semi-implicit integration scheme
							LatticeDensity ghostDensityNew = ((R0/(1.0 + R0*C0))*(scaleFactor*std::abs(component) + C0*pressure))/Cs2;

							wkIolet->SetDensityNew(ghostDensityNew);
						}

						// TODO it's ugly that we have to do this.
						// TODO having to give 0 as an argument is also ugly.
						// TODO it's ugly that we have to give hydroVars a nonsense distribution vector
						// that doesn't get used.
						kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(site);

						ghostHydrovars.density = ghostDensity;
						ghostHydrovars.momentum = ioletNormal * component * ghostDensity;

						collider.kernel.CalculateFeq(ghostHydrovars, 0);

						Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

						*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed)
							= ghostHydrovars.GetFEq()[unstreamed];
					}
					
					inline void PostStepLink(geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();

#ifdef HEMELB_USE_VELOCITY_WEIGHTS_FILE
						iolets::InOutLetFileWK* wkIolet = dynamic_cast<iolets::InOutLetFileWK*>(iolet.GetIolets()[boundaryId]);
#else
						iolets::InOutLetWK* wkIolet = dynamic_cast<iolets::InOutLetWK*>(iolet.GetIolets()[boundaryId]);
#endif
						wkIolet->SetDensity(wkIolet->GetDensityNew(0.0));
					}

				protected:
					CollisionType& collider;
					iolets::BoundaryValues& iolet;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKDELEGATE_H
