
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
                                                iolets::InOutLetFileWK* wkIolet = dynamic_cast<iolets::InOutLetFileWK*>(iolet.GetLocalIolet(boundaryId));
#else                                               
					       	iolets::InOutLetWK* wkIolet = dynamic_cast<iolets::InOutLetWK*>(iolet.GetLocalIolet(boundaryId));
#endif						
						distribn_t Ptm1;
						double R0 = wkIolet->GetRwk();	
					        
						Ptm1 = (wkIolet->GetDensity(iolet.GetTimeStep())-1.0)*Cs2;
						
						LatticePosition sitePos(site.GetGlobalSiteCoords());

				        	LatticePosition halfWay(sitePos);
						halfWay.x += 0.5 * LatticeType::CX[direction];
						halfWay.y += 0.5 * LatticeType::CY[direction];
						halfWay.z += 0.5 * LatticeType::CZ[direction];

						distribn_t scaleFactor = wkIolet->GetQtScaleFactor(halfWay); 
						distribn_t distance = wkIolet->GetDistance(halfWay); 

						// Calculate the velocity at the ghost site, as the component normal to the iolet.
						util::Vector3D<float> ioletNormal = iolet.GetLocalIolet(boundaryId)->GetNormal();

						// Note that the division by density compensates for the fact that v_x etc have momentum
						// not velocity.
						distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal); 						
						distribn_t ghostDensity = wkIolet->GetDensity(0.0);

						if (distance < 0.5)
						{

						// Set the density at the "ghost" site to be the density prediceted by the 2-element WK.
						ghostDensity = ((R0/(1.0 + 0.2/lbmParams->GetTimeStep()))*(scaleFactor*std::abs(component) + 0.2*Ptm1/(lbmParams->GetTimeStep()*R0)))/Cs2 + 1.0;
						
						wkIolet->SetDensity(ghostDensity);
						}

if (iolet.GetTimeStep()%50==0&&distance<0.5) //(sitePos.x==7 && sitePos.y==13 && sitePos.z==240)||(sitePos.x==193&&sitePos.y==13&&sitePos.z==131)))
{
	std::cout << "x, scalefactor, velocity, R, Ptm1, ghostDensity :" << halfWay.x << "," << scaleFactor << "," << std::abs(component) << "," << R0 << "," << Ptm1 << "," << ghostDensity << std::endl;
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
				protected:
					CollisionType& collider;
					iolets::BoundaryValues& iolet;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_GRINBERGKARNIADAKISWKDELEGATE_H
