// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_LB_HPP
#define HEMELB_LB_LB_HPP

#include "io/writers/xdr/XdrMemWriter.h"
#include "lb/lb.h"
#include "util/unique.h"
#include "lb/InitialCondition.h"
#include "lb/InitialCondition.hpp"

namespace hemelb
{
	namespace lb
	{

		template<class LatticeType>
			hemelb::lb::LbmParameters* LBM<LatticeType>::GetLbmParams()
			{
				return &mParams;
			}

		template<class LatticeType>
			lb::MacroscopicPropertyCache& LBM<LatticeType>::GetPropertyCache()
			{
				return propertyCache;
			}

		template<class LatticeType>
			LBM<LatticeType>::LBM(configuration::SimConfig *iSimulationConfig,
					net::Net* net,
					geometry::LatticeData* latDat,
					SimulationState* simState,
					reporting::Timers &atimings,
					geometry::neighbouring::NeighbouringDataManager *neighbouringDataManager) :
				mSimConfig(iSimulationConfig), mNet(net), mLatDat(latDat), mState(simState),
				mParams(iSimulationConfig->GetTimeStepLength(), iSimulationConfig->GetVoxelSize()), timings(atimings),
				propertyCache(*simState, *latDat), neighbouringDataManager(neighbouringDataManager)
		{
			ReadParameters();
		}

		template<class LatticeType>
			void LBM<LatticeType>::InitInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state)
			{
				initParams.siteRanges.resize(2);

				initParams.siteRanges[0].first = 0;
				initParams.siteRanges[1].first = mLatDat->GetMidDomainSiteCount();
				state = 0;
				initParams.siteRanges[0].second = initParams.siteRanges[0].first + mLatDat->GetMidDomainCollisionCount(state);
				initParams.siteRanges[1].second = initParams.siteRanges[1].first + mLatDat->GetDomainEdgeCollisionCount(state);

				initParams.siteCount = mLatDat->GetMidDomainCollisionCount(state) + mLatDat->GetDomainEdgeCollisionCount(state);
			}

		template<class LatticeType>
			void LBM<LatticeType>:: AdvanceInitParamsSiteRanges(kernels::InitParams& initParams, unsigned& state)
			{
				initParams.siteRanges[0].first += mLatDat->GetMidDomainCollisionCount(state);
				initParams.siteRanges[1].first += mLatDat->GetDomainEdgeCollisionCount(state);
				++state;
				initParams.siteRanges[0].second = initParams.siteRanges[0].first + mLatDat->GetMidDomainCollisionCount(state);
				initParams.siteRanges[1].second = initParams.siteRanges[1].first + mLatDat->GetDomainEdgeCollisionCount(state);

				initParams.siteCount = mLatDat->GetMidDomainCollisionCount(state) + mLatDat->GetDomainEdgeCollisionCount(state);
			}

		template<class LatticeType>
			void LBM<LatticeType>::InitCollisions()
			{
				/**
				 * Ensure the boundary objects have all info necessary.
				 */
				PrepareBoundaryObjects();

				// TODO Note that the convergence checking is not yet implemented in the
				// new boundary condition hierarchy system.
				// It'd be nice to do this with something like
				// MidFluidCollision = new ConvergenceCheckingWrapper(new WhateverMidFluidCollision());

				kernels::InitParams initParams = kernels::InitParams();
				initParams.latDat = mLatDat;
				initParams.lbmParams = &mParams;
				initParams.neighbouringDataManager = neighbouringDataManager;
				for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
				{
					initParams.outletPositions.push_back(mOutletValues->GetLocalIolet(outlet)->GetPosition());
				}

				unsigned collId;
				InitInitParamsSiteRanges(initParams, collId);
				mMidFluidCollision = new tMidFluidCollision(initParams);

				AdvanceInitParamsSiteRanges(initParams, collId);
				mWallCollision = new tWallCollision(initParams);

				AdvanceInitParamsSiteRanges(initParams, collId);
				initParams.boundaryObject = mInletValues;
				mInletCollision = new tInletCollision(initParams);

				AdvanceInitParamsSiteRanges(initParams, collId);
				initParams.boundaryObject = mOutletValues;
				mOutletCollision = new tOutletCollision(initParams);

				AdvanceInitParamsSiteRanges(initParams, collId);
				initParams.boundaryObject = mInletValues;
				mInletWallCollision = new tInletWallCollision(initParams);

				AdvanceInitParamsSiteRanges(initParams, collId);
				initParams.boundaryObject = mOutletValues;
				mOutletWallCollision = new tOutletWallCollision(initParams);
			}

		template<class LatticeType>
			void LBM<LatticeType>::Initialise(iolets::BoundaryValues* iInletValues,
					iolets::BoundaryValues* iOutletValues,
					const util::UnitConverter* iUnits)
			{
				mInletValues = iInletValues;
				mOutletValues = iOutletValues;
				mUnits = iUnits;

				InitCollisions();

				//JM SetInitialConditions();
			}

		template<class LatticeType>
			void LBM<LatticeType>::PrepareBoundaryObjects()
			{
				// First, iterate through all of the inlet and outlet objects, finding out the minimum density seen in the simulation.
				distribn_t minDensity = std::numeric_limits<distribn_t>::max();

				for (unsigned inlet = 0; inlet < mInletValues->GetLocalIoletCount(); ++inlet)
				{
					minDensity = std::min(minDensity, mInletValues->GetLocalIolet(inlet)->GetDensityMin());
				}

				for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
				{
					minDensity = std::min(minDensity, mOutletValues->GetLocalIolet(outlet)->GetDensityMin());
				}

				// Now go through them again, informing them of the minimum density.
				for (unsigned inlet = 0; inlet < mInletValues->GetLocalIoletCount(); ++inlet)
				{
					mInletValues->GetLocalIolet(inlet)->SetMinimumSimulationDensity(minDensity);
				}

				for (unsigned outlet = 0; outlet < mOutletValues->GetLocalIoletCount(); ++outlet)
				{
					mOutletValues->GetLocalIolet(outlet)->SetMinimumSimulationDensity(minDensity);
				}
			}

		template<class LatticeType>
			
			    void LBM<LatticeType>::SetInitialConditions(const net::IOCommunicator& ioComms)
			    {
			      auto icond = InitialCondition::FromConfig(mSimConfig->GetInitialCondition());
			      icond.SetFs<LatticeType>(mLatDat, ioComms);
			      icond.SetTime(mState);
			    }
		/** JM Method before trying to bring Checkpointing in	
			void LBM<LatticeType>::SetInitialConditions()
			{
				distribn_t density = mUnits->ConvertPressureToLatticeUnits(mSimConfig->GetInitialPressure()) / Cs2;

				for (site_t i = 0; i < mLatDat->GetLocalFluidSiteCount(); i++)
				{
					distribn_t f_eq[LatticeType::NUMVECTORS];

					LatticeType::CalculateFeq(density, 0.0, 0.0, 0.0, f_eq);

					distribn_t* f_old_p = mLatDat->GetFOld(i * LatticeType::NUMVECTORS);
					distribn_t* f_new_p = mLatDat->GetFNew(i * LatticeType::NUMVECTORS);

					for (unsigned int l = 0; l < LatticeType::NUMVECTORS; l++)
					{
						f_new_p[l] = f_old_p[l] = f_eq[l];
					}
				}
			}
		**/

		template<class LatticeType>
			void LBM<LatticeType>::RequestComms()
			{
				timings[hemelb::reporting::Timers::lb].Start();

				// Delegate to the lattice data object to post the asynchronous sends and receives
				// (via the Net object).
				// NOTE that this doesn't actually *perform* the sends and receives, it asks the Net
				// to include them in the ISends and IRecvs that happen later.
				mLatDat->SendAndReceive(mNet);

				timings[hemelb::reporting::Timers::lb].Stop();
			}

		template<class LatticeType>
			void LBM<LatticeType>::PreSend()
			{
				timings[hemelb::reporting::Timers::lb].Start();
				timings[hemelb::reporting::Timers::lb_calc].Start();

				/**
				 * In the PreSend phase, we do LB on all the sites that need to have results sent to
				 * neighbouring ranks ('domainEdge' sites). In site id terms, this means we start at the
				 * end of the sites whose neighbours all lie on this rank ('midDomain'), then progress
				 * through the sites of each type in turn.
				 */
				site_t offset = mLatDat->GetMidDomainSiteCount();

				StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetDomainEdgeCollisionCount(0));
				offset += mLatDat->GetDomainEdgeCollisionCount(0);

				StreamAndCollide(mWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(1));
				offset += mLatDat->GetDomainEdgeCollisionCount(1);

				mInletValues->FinishReceive();
				StreamAndCollide(mInletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(2));
				offset += mLatDat->GetDomainEdgeCollisionCount(2);

				mOutletValues->FinishReceive();
				StreamAndCollide(mOutletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(3));
				offset += mLatDat->GetDomainEdgeCollisionCount(3);

				StreamAndCollide(mInletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(4));
				offset += mLatDat->GetDomainEdgeCollisionCount(4);

				StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(5));

				timings[hemelb::reporting::Timers::lb_calc].Stop();
				timings[hemelb::reporting::Timers::lb].Stop();
			}

		template<class LatticeType>
			void LBM<LatticeType>::PreReceive()
			{
				timings[hemelb::reporting::Timers::lb].Start();
				timings[hemelb::reporting::Timers::lb_calc].Start();

				/**
				 * In the PreReceive phase, we perform LB for all the sites whose neighbours lie on this
				 * rank ('midDomain' rather than 'domainEdge' sites). Ideally this phase is the longest bit (maximising time for the asynchronous sends
				 * and receives to complete).
				 *
				 * In site id terms, this means starting at the first site and progressing through the
				 * midDomain sites, one type at a time.
				 */
				site_t offset = 0;

				StreamAndCollide(mMidFluidCollision, offset, mLatDat->GetMidDomainCollisionCount(0));
				offset += mLatDat->GetMidDomainCollisionCount(0);

				StreamAndCollide(mWallCollision, offset, mLatDat->GetMidDomainCollisionCount(1));
				offset += mLatDat->GetMidDomainCollisionCount(1);

				StreamAndCollide(mInletCollision, offset, mLatDat->GetMidDomainCollisionCount(2));
				offset += mLatDat->GetMidDomainCollisionCount(2);

				StreamAndCollide(mOutletCollision, offset, mLatDat->GetMidDomainCollisionCount(3));
				offset += mLatDat->GetMidDomainCollisionCount(3);

				StreamAndCollide(mInletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(4));
				offset += mLatDat->GetMidDomainCollisionCount(4);

				StreamAndCollide(mOutletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(5));

				timings[hemelb::reporting::Timers::lb_calc].Stop();
				timings[hemelb::reporting::Timers::lb].Stop();
			}

		template<class LatticeType>
			void LBM<LatticeType>::PostReceive()
			{
				timings[hemelb::reporting::Timers::lb].Start();

				// Copy the distribution functions received from the neighbouring
				// processors into the destination buffer "f_new".
				// This is done here, after receiving the sent distributions from neighbours.
				mLatDat->CopyReceived();

				// Do any cleanup steps necessary on boundary nodes
				site_t offset = mLatDat->GetMidDomainSiteCount();

				timings[hemelb::reporting::Timers::lb_calc].Start();

				//TODO yup, this is horrible. If you read this, please improve the following code.
				PostStep(mMidFluidCollision, offset, mLatDat->GetDomainEdgeCollisionCount(0));
				offset += mLatDat->GetDomainEdgeCollisionCount(0);

				PostStep(mWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(1));
				offset += mLatDat->GetDomainEdgeCollisionCount(1);

				PostStep(mInletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(2));
				offset += mLatDat->GetDomainEdgeCollisionCount(2);

				PostStep(mOutletCollision, offset, mLatDat->GetDomainEdgeCollisionCount(3));
				offset += mLatDat->GetDomainEdgeCollisionCount(3);

				PostStep(mInletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(4));
				offset += mLatDat->GetDomainEdgeCollisionCount(4);

				PostStep(mOutletWallCollision, offset, mLatDat->GetDomainEdgeCollisionCount(5));

				offset = 0;

				PostStep(mMidFluidCollision, offset, mLatDat->GetMidDomainCollisionCount(0));
				offset += mLatDat->GetMidDomainCollisionCount(0);

				PostStep(mWallCollision, offset, mLatDat->GetMidDomainCollisionCount(1));
				offset += mLatDat->GetMidDomainCollisionCount(1);

				PostStep(mInletCollision, offset, mLatDat->GetMidDomainCollisionCount(2));
				offset += mLatDat->GetMidDomainCollisionCount(2);

				PostStep(mOutletCollision, offset, mLatDat->GetMidDomainCollisionCount(3));
				offset += mLatDat->GetMidDomainCollisionCount(3);

				PostStep(mInletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(4));
				offset += mLatDat->GetMidDomainCollisionCount(4);

				PostStep(mOutletWallCollision, offset, mLatDat->GetMidDomainCollisionCount(5));

				timings[hemelb::reporting::Timers::lb_calc].Stop();
				timings[hemelb::reporting::Timers::lb].Stop();
			}

		template<class LatticeType>
			void LBM<LatticeType>::EndIteration()
			{
				timings[hemelb::reporting::Timers::lb].Start();
				timings[hemelb::reporting::Timers::lb_calc].Start();

				// Swap f_old and f_new ready for the next timestep.
				mLatDat->SwapOldAndNew();

				timings[hemelb::reporting::Timers::lb_calc].Stop();
				timings[hemelb::reporting::Timers::lb].Stop();
			}

		template<class LatticeType>
			LBM<LatticeType>::~LBM()
			{
				// Delete the collision and stream objects we've been using
				delete mMidFluidCollision;
				delete mWallCollision;
				delete mInletCollision;
				delete mOutletCollision;
				delete mInletWallCollision;
				delete mOutletWallCollision;
			}

		template<class LatticeType>
			void LBM<LatticeType>::ReadParameters()
			{
				std::vector<lb::iolets::InOutLet*> inlets = mSimConfig->GetInlets();
				std::vector<lb::iolets::InOutLet*> outlets = mSimConfig->GetOutlets();
				inletCount = inlets.size();
				outletCount = outlets.size();
				mParams.StressType = mSimConfig->GetStressType();
				mParams.SetRelaxationParameter(mSimConfig->GetRelaxationParameter());
				mParams.ElasticWallStiffness = mSimConfig->GetElasticWallStiffness();
				mParams.BoundaryVelocityRatio = mSimConfig->GetBoundaryVelocityRatio();
				mParams.ViscosityRatio = mSimConfig->GetViscosityRatio();
				mParams.SpongeLayerWidth = mSimConfig->GetSpongeLayerWidth();
			}
	}
}

#endif /* HEMELB_LB_LB_HPP */
