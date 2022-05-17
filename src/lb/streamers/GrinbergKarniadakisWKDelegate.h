
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
						iolets::InOutLet* localIOlet = iolet.GetIolets()[boundaryId];
						localIOlet->DoPreStreamCoupling(site.GetIndex(), site.GetGlobalSiteCoords(),
														hydroVars.density, hydroVars.velocity);
						nashDelegate.StreamLink(lbmParams, latticeData, site, hydroVars, direction);
					}
					
					inline void PostStepLink(geometry::LatticeData* const latticeData,
							const geometry::Site<geometry::LatticeData>& site,
							const Direction& direction)
					{
						int boundaryId = site.GetIoletId();
						iolets::InOutLet* localIOlet = iolet.GetIolets()[boundaryId];
						localIOlet->DoPostStreamCoupling(site.GetIndex(), site.GetGlobalSiteCoords());
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
