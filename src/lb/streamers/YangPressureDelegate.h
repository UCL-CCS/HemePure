
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_YANGPRESSUREDELEGATE_H
#define HEMELB_LB_STREAMERS_YANGPRESSUREDELEGATE_H

#include "util/utilityFunctions.h"
#include "lb/streamers/BaseStreamerDelegate.h"

namespace hemelb
{
	namespace lb
	{
		namespace streamers
		{
			/**
			 * This class implements the pressure boundary condition described by Zhaoxia Yang in
			 * "Pressure condition for lattice Boltzmann methods on domains with curved boundaries",
			 * Computers and Mathematics with Applications 59 (2010) 2168-2177.
			 *
			 * The following assumptions are applied to simplify the implementation.
			 * 1. The iolet plane is straight rather than curved.
			 * 2. The collision kernel adopts a single relaxation time.
			 */
			template <typename CollisionImpl>
			class YangPressureDelegate : public BaseStreamerDelegate<CollisionImpl>
			{
			public:
				typedef CollisionImpl CollisionType;
				typedef typename CollisionType::CKernel::LatticeType LatticeType;

				YangPressureDelegate(CollisionType &delegatorCollider, kernels::InitParams &initParams) :
					collider(delegatorCollider),
					neighbouringLatticeData(initParams.latDat->GetNeighbouringData()),
					iolet(*initParams.boundaryObject)
				{
					// Find the lattice directions described by the mapping G and P.
					for (int i = 0; i < iolet.GetLocalIoletCount(); ++i)
					{
						SortDirectionsCloseToIOletNormal(iolet.GetLocalIolet(i));
					}

					// Loop over each site this streamer is responsible for, as specified in siteRanges.
					for (std::vector<std::pair<site_t, site_t>>::iterator rangeIt =
							 initParams.siteRanges.begin();
						 rangeIt != initParams.siteRanges.end(); ++rangeIt)
					{
						// Part of the loop mentioned
						for (site_t localIndex = rangeIt->first; localIndex < rangeIt->second; ++localIndex)
						{
							geometry::Site<const geometry::LatticeData> localSite =
								initParams.latDat->GetSite(localIndex);
							const LatticeVector &localSiteLocation = localSite.GetGlobalSiteCoords();
							const iolets::InOutLet* localIOlet = iolet.GetIolets()[localSite.GetIoletId()];

							// Ensure the data of the two fluid sites used for extrapolation is available.
							for (Direction i = 1; i < LatticeType::NUMVECTORS; ++i)
							{
								/**
								 * Taking one step from the current site in the current direction, i, gives an
								 * outer-wall (solid) node. Taking a further step in the direction mapped
								 * by P gives the first fluid site. Taking one more step in this direction
								 * gives the second fluid site.
								 */

								// Skip if there is no iolet in this direction.
								if (!localSite.HasIolet(i)) continue;

								const LatticeVector ci = LatticeVector(LatticeType::CX[i],
																	   LatticeType::CY[i],
																	   LatticeType::CZ[i]);
								const LatticeVector outerWallSiteLocation = localSiteLocation + ci;

								Direction j = 0, dirP = localIOlet->GetDirectionCloseToNormal(0);
								do
								{
									const LatticeVector cp = LatticeVector(LatticeType::CX[dirP],
																		   LatticeType::CY[dirP],
																		   LatticeType::CZ[dirP]);

									const LatticeVector firstFluidSiteLocation = outerWallSiteLocation + cp;
									proc_t firstFluidSiteHomeProc =
										initParams.latDat->GetProcIdFromGlobalCoords(firstFluidSiteLocation);

									const LatticeVector secondFluidSiteLocation = firstFluidSiteLocation + cp;
									proc_t secondFluidSiteHomeProc =
										initParams.latDat->GetProcIdFromGlobalCoords(secondFluidSiteLocation);

									if (firstFluidSiteHomeProc != SITE_OR_BLOCK_SOLID &&
										secondFluidSiteHomeProc != SITE_OR_BLOCK_SOLID)
									{
										if (firstFluidSiteHomeProc != initParams.latDat->GetLocalRank())
										{
											// Request data from another rank.
											initParams.neighbouringDataManager->RegisterNeededSite(
												initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(firstFluidSiteLocation));
										}
										if (secondFluidSiteHomeProc != initParams.latDat->GetLocalRank())
										{
											// Request data from another rank.
											initParams.neighbouringDataManager->RegisterNeededSite(
												initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(secondFluidSiteLocation));
										}
										break;
									}

									j++;
									dirP = localIOlet->GetDirectionCloseToNormal(j);
								}
								while (dirP != 0);

								/**
								 * Since the directions are sorted, further searches cannot be successful.
								 * Try again by taking two steps instead of one step from the outer-wall node.
								 */
								if (dirP == 0)
								{
									j = 0;
									dirP = localIOlet->GetDirectionCloseToNormal(0);
									do
									{
										const LatticeVector cp = LatticeVector(LatticeType::CX[dirP],
																			   LatticeType::CY[dirP],
																			   LatticeType::CZ[dirP]);

										// The only difference is the factor 2 in the next line.
										const LatticeVector firstFluidSiteLocation = outerWallSiteLocation + cp * 2;
										proc_t firstFluidSiteHomeProc =
											initParams.latDat->GetProcIdFromGlobalCoords(firstFluidSiteLocation);

										const LatticeVector secondFluidSiteLocation = firstFluidSiteLocation + cp;
										proc_t secondFluidSiteHomeProc =
											initParams.latDat->GetProcIdFromGlobalCoords(secondFluidSiteLocation);

										if (firstFluidSiteHomeProc != SITE_OR_BLOCK_SOLID &&
											secondFluidSiteHomeProc != SITE_OR_BLOCK_SOLID)
										{
											if (firstFluidSiteHomeProc != initParams.latDat->GetLocalRank())
											{
												// Request data from another rank.
												initParams.neighbouringDataManager->RegisterNeededSite(
													initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(firstFluidSiteLocation));
											}
											if (secondFluidSiteHomeProc != initParams.latDat->GetLocalRank())
											{
												// Request data from another rank.
												initParams.neighbouringDataManager->RegisterNeededSite(
													initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(secondFluidSiteLocation));
											}
											break;
										}

										j++;
										dirP = localIOlet->GetDirectionCloseToNormal(j);
									}
									while (dirP != 0);

									// If the search fails again, a finer grid should be used.
									if (dirP == 0)
									{
										hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::OnePerCore>(
											"A higher resolution is required at the site location [%ld, %ld, %ld]",
											localSiteLocation.x, localSiteLocation.y, localSiteLocation.z);
										std::exit(15);
									}
								}
							}
						}
					}
				}

				inline void
				StreamLink(const LbmParameters *lbmParams,
						   geometry::LatticeData *const latticeData,
						   const geometry::Site<geometry::LatticeData> &site,
						   kernels::HydroVars<typename CollisionType::CKernel> &hydroVars,
						   const Direction &direction)
				{
					iolets::InOutLet* localIOlet = iolet.GetIolets()[site.GetIoletId()];
					const LatticePosition& ioletNormal = localIOlet->GetNormal();
					Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];
					Direction dirG = localIOlet->GetDirectionCloseToNormal(0);

					// Couple with an external system if there is
					if (unstreamed == dirG)
					{
						localIOlet->DoPreStreamCoupling(site.GetIndex(), iolet.GetTimeStep(), site.GetGlobalSiteCoords(),
														hydroVars.density, hydroVars.velocity);
					}

					// Obtain cp, wallDistance, and the old distributions at the first and second fluid nodes.
					// Here wall is referred to as the iolet plane.
					bool special; // whether the outer-wall node has no neighbouring fluid site
					LatticeVector cp; // the lattice vector given by the mapping P
					LatticeDistance wallDistance; // distance between the wall and the first fluid site
					const distribn_t *firstFluidFOld, *secondFluidFOld;
					GetFluidSitesData(site, direction, latticeData, special, cp, wallDistance, firstFluidFOld, secondFluidFOld);

					// Calculate the densities, momenta, and equilibrium distributions and store them in HydroVars.
					kernels::HydroVars<typename CollisionType::CKernel> hVfirstFluid(firstFluidFOld);
					collider.kernel.CalculateDensityMomentumFeq(hVfirstFluid, 0); // the second argument is dummy
					kernels::HydroVars<typename CollisionType::CKernel> hVsecondFluid(secondFluidFOld);
					collider.kernel.CalculateDensityMomentumFeq(hVsecondFluid, 0); // the second argument is dummy

					// Calculate the non-equilibrium distributions on the wall (equation 16).
					distribn_t fNeqWall[LatticeType::NUMVECTORS];
					for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i)
					{
						fNeqWall[i] = (1.0 + wallDistance) * hVfirstFluid.GetFNeq().f[i]
									  - wallDistance * hVsecondFluid.GetFNeq().f[i];
					}

					LatticeVector cg = LatticeVector(LatticeType::CX[dirG], LatticeType::CY[dirG], LatticeType::CZ[dirG]);
					distribn_t fNew;

					if (unstreamed == dirG && cg == cp)
					{
						// Calculate the density on the wall (equation 13).
						// The equilibrium density, 1, is absorbed in the iolet pressure.
						// Also note that h = 1 in lattice units.
						const distribn_t visc = Cs2 * (hydroVars.tau - 0.5); // kinematic viscosity
						LatticeDensity densityWall = 3.0 * (localIOlet->GetPressure(iolet.GetTimeStep())
															+ visc * stress(-ioletNormal, hydroVars.tau, fNeqWall));

						// Interpolate the density at the first fluid site accounting for the BC (equation 12).
						LatticeDensity densityBC = (densityWall + wallDistance * hVsecondFluid.density) / (1.0 + wallDistance);

						// Calculate the post-collision distributions at the second fluid site.
						collider.Collide(lbmParams, hVsecondFluid);

						// Calculate the distribution of the unstreamed direction (equation 11).
						fNew = site.GetFOld<LatticeType>()[unstreamed]
								+ site.GetFOld<LatticeType>()[direction]
								+ 2.0 * LatticeType::EQMWEIGHTS[unstreamed] * (densityBC - hydroVars.density)
								- hVsecondFluid.GetFPostCollision()[direction];
					}
					else
					{
						// Construct a HydroVars structure at the outer-wall node.
						distribn_t fOuterWall[LatticeType::NUMVECTORS];
						kernels::HydroVars<typename CollisionType::CKernel> hVouterWall(fOuterWall);

						// Calculate the required quantities at the outer-wall node.
						const LatticePosition sqBracket = momentumCorrection(cp, -ioletNormal, hydroVars.tau, fNeqWall);
						distribn_t fNeqOuterWall;
						if (special)
						{
							// Extrapolate the density and momentum (equation B.4).
							hVouterWall.density = 3.0 * hVfirstFluid.density - 2.0 * hVsecondFluid.density;
							hVouterWall.momentum = (hVfirstFluid.momentum * (6.0 * wallDistance - 3.0)
													+ hVsecondFluid.momentum * (4.0 - 4.0 * wallDistance)
													- sqBracket * 3.0
												   ) / (1.0 + 2.0 * wallDistance);

							// Extrapolate the non-equilibrium distribution of the unstreamed direction (equation B.4).
							fNeqOuterWall = 3.0 * hVfirstFluid.GetFNeq().f[unstreamed]
											- 2.0 * hVsecondFluid.GetFNeq().f[unstreamed];
						}
						else
						{
							// Extrapolate the density and momentum (equation 18).
							hVouterWall.density = 2.0 * hVfirstFluid.density - hVsecondFluid.density;
							hVouterWall.momentum = (hVfirstFluid.momentum * 4.0 * wallDistance
													+ hVsecondFluid.momentum * (1.0 - 2.0 * wallDistance)
													- sqBracket
												   ) / (1.0 + 2.0 * wallDistance);

							// Extrapolate the non-equilibrium distribution of the unstreamed direction (equation 20).
							fNeqOuterWall = 2.0 * hVfirstFluid.GetFNeq().f[unstreamed]
											- hVsecondFluid.GetFNeq().f[unstreamed];
						}

						// Calculate the equilibrium distributions at the outer-wall node.
						LatticeType::CalculateFeq(hVouterWall.density,
												  hVouterWall.momentum.x,
												  hVouterWall.momentum.y,
												  hVouterWall.momentum.z,
												  hVouterWall.GetFEqPtr());

						// Calculate the distribution of the unstreamed direction (equation 17).
						// Assumption 2 is applied here: A equals 1/tau times the identity matrix.
						fNew = hVouterWall.GetFEq().f[unstreamed] + fNeqOuterWall * (1.0 - 1.0 / hydroVars.tau);
					}

					*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed) = fNew;
				}

				inline void PostStepLink(geometry::LatticeData* const latticeData,
						const geometry::Site<geometry::LatticeData>& site,
						const Direction& direction)
				{
					int boundaryId = site.GetIoletId();
					iolets::InOutLet* localIOlet = iolet.GetIolets()[boundaryId];
					Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];
					Direction dirG = localIOlet->GetDirectionCloseToNormal(0);

					// Finalise the coupling with the external system
					if (unstreamed == dirG)
					{
						localIOlet->DoPostStreamCoupling(site.GetIndex(), iolet.GetTimeStep(), site.GetGlobalSiteCoords());
					}
				}

			protected:
				/**
				 * Sort the lattice directions in an ascending order of how well they allign with
				 * the iolet normal vector. See equation A.1 in the paper. Due to assumption 1,
				 * the directions are the same for all sites linked to the same iolet.
				 */
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

				/**
				 * Obtain the data at the two fluid sites. The logic used in the constructor is applied.
				 */
				void GetFluidSitesData(const geometry::Site<geometry::LatticeData> &site,
									   const Direction &i,
									   geometry::LatticeData *const latDat,
									   bool &special,
									   LatticeVector &cp,
									   LatticeDistance &firstFluidWallDistance,
									   const distribn_t* &firstFluidFOld,
									   const distribn_t* &secondFluidFOld)
				{
					const iolets::InOutLet* localIOlet = iolet.GetIolets()[site.GetIoletId()];

					// The parameter i is assumed to be the stream direction.
					const LatticeVector ci = LatticeVector(LatticeType::CX[i],
														   LatticeType::CY[i],
														   LatticeType::CZ[i]);
					const LatticeVector outerWallSiteLocation = site.GetGlobalSiteCoords() + ci;

					special = false;
					Direction j = 0, dirP = localIOlet->GetDirectionCloseToNormal(0);
					do
					{
						cp = LatticeVector(LatticeType::CX[dirP], LatticeType::CY[dirP], LatticeType::CZ[dirP]);

						const LatticeVector firstFluidSiteLocation = outerWallSiteLocation + cp;
						proc_t firstFluidSiteHomeProc = latDat->GetProcIdFromGlobalCoords(firstFluidSiteLocation);

						const LatticeVector secondFluidSiteLocation = firstFluidSiteLocation + cp;
						proc_t secondFluidSiteHomeProc = latDat->GetProcIdFromGlobalCoords(secondFluidSiteLocation);

						if (firstFluidSiteHomeProc != SITE_OR_BLOCK_SOLID &&
							secondFluidSiteHomeProc != SITE_OR_BLOCK_SOLID)
						{
							const Direction oppP = LatticeType::INVERSEDIRECTIONS[dirP];
							if (firstFluidSiteHomeProc == latDat->GetLocalRank())
							{
								geometry::Site<geometry::LatticeData> firstFluidSite =
									latDat->GetSite(latDat->GetContiguousSiteId(firstFluidSiteLocation));
								firstFluidFOld = firstFluidSite.GetFOld<LatticeType>();
								firstFluidWallDistance = firstFluidSite.GetWallDistance<LatticeType>(oppP);
							}
							else
							{
								const geometry::neighbouring::ConstNeighbouringSite
									firstFluidSite = neighbouringLatticeData.GetSite(
										latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(firstFluidSiteLocation));
								firstFluidFOld = firstFluidSite.GetFOld<LatticeType>();
								firstFluidWallDistance = firstFluidSite.GetWallDistance<LatticeType>(oppP);
							}

							if (secondFluidSiteHomeProc == latDat->GetLocalRank())
							{
								geometry::Site<geometry::LatticeData> secondFluidSite =
									latDat->GetSite(latDat->GetContiguousSiteId(secondFluidSiteLocation));
								secondFluidFOld = secondFluidSite.GetFOld<LatticeType>();
							}
							else
							{
								const geometry::neighbouring::ConstNeighbouringSite
									secondFluidSite = neighbouringLatticeData.GetSite(
										latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(secondFluidSiteLocation));
								secondFluidFOld = secondFluidSite.GetFOld<LatticeType>();
							}
							break;
						}

						j++;
						dirP = localIOlet->GetDirectionCloseToNormal(j);
					}
					while (dirP != 0);

					if (dirP == 0)
					{
						special = true;
						j = 0;
						dirP = localIOlet->GetDirectionCloseToNormal(0);
						do
						{
							cp = LatticeVector(LatticeType::CX[dirP], LatticeType::CY[dirP], LatticeType::CZ[dirP]);

							// The only difference is the factor 2 in the next line.
							const LatticeVector firstFluidSiteLocation = outerWallSiteLocation + cp * 2;
							proc_t firstFluidSiteHomeProc = latDat->GetProcIdFromGlobalCoords(firstFluidSiteLocation);

							const LatticeVector secondFluidSiteLocation = firstFluidSiteLocation + cp;
							proc_t secondFluidSiteHomeProc = latDat->GetProcIdFromGlobalCoords(secondFluidSiteLocation);

							if (firstFluidSiteHomeProc != SITE_OR_BLOCK_SOLID &&
								secondFluidSiteHomeProc != SITE_OR_BLOCK_SOLID)
							{
								const Direction oppP = LatticeType::INVERSEDIRECTIONS[dirP];
								if (firstFluidSiteHomeProc == latDat->GetLocalRank())
								{
									geometry::Site<geometry::LatticeData> firstFluidSite =
										latDat->GetSite(latDat->GetContiguousSiteId(firstFluidSiteLocation));
									firstFluidFOld = firstFluidSite.GetFOld<LatticeType>();
									firstFluidWallDistance = firstFluidSite.GetWallDistance<LatticeType>(oppP);
								}
								else
								{
									const geometry::neighbouring::ConstNeighbouringSite
										firstFluidSite = neighbouringLatticeData.GetSite(
											latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(firstFluidSiteLocation));
									firstFluidFOld = firstFluidSite.GetFOld<LatticeType>();
									firstFluidWallDistance = firstFluidSite.GetWallDistance<LatticeType>(oppP);
								}

								if (secondFluidSiteHomeProc == latDat->GetLocalRank())
								{
									geometry::Site<geometry::LatticeData> secondFluidSite =
										latDat->GetSite(latDat->GetContiguousSiteId(secondFluidSiteLocation));
									secondFluidFOld = secondFluidSite.GetFOld<LatticeType>();
								}
								else
								{
									const geometry::neighbouring::ConstNeighbouringSite
										secondFluidSite = neighbouringLatticeData.GetSite(
											latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(secondFluidSiteLocation));
									secondFluidFOld = secondFluidSite.GetFOld<LatticeType>();
								}
								break;
							}

							j++;
							dirP = localIOlet->GetDirectionCloseToNormal(j);
						}
						while (dirP != 0);
					}
				}

				/**
				 * This function calculates the value of B_k(s) in equation 15.
				 */
				inline distribn_t B_k(const Direction k, const LatticePosition s)
				{
					const LatticePosition ck = LatticePosition(LatticeType::CXD[k],
															   LatticeType::CYD[k],
															   LatticeType::CZD[k]);
					distribn_t sc2 = s.Dot(ck);
					sc2 *= sc2;
					const distribn_t s2 = s.GetMagnitudeSquared();
					const distribn_t c2 = ck.GetMagnitudeSquared();
					return (sc2 - s2 * c2 / 3.0) / (2.0 * Cs2);
				}

				/**
				 * This function calculates the stress term in equation 14.
				 * Assumption 2 is applied here: the collision matrix A equals 1/tau
				 * times the identity matrix, where tau is the single relaxation time.
				 */
				inline distribn_t stress(const LatticePosition vec,
										 const distribn_t tau,
										 const distribn_t *const fNeq)
				{
					distribn_t stress = 0;
					for (Direction k = 0; k < LatticeType::NUMVECTORS; ++k)
					{
						stress += B_k(k, vec) * (fNeq[k] + fNeq[LatticeType::INVERSEDIRECTIONS[k]]);
					}
					// Note that h = 1 in lattice units.
					stress *= -1.0 / (2.0 * tau);
					return stress;
				}

				inline LatticePosition tangentialVector(const LatticePosition vec, const LatticePosition normal)
				{
					distribn_t normalComp = vec.Dot(normal);
					return (vec - normal * normalComp).Normalise();
				}

				/**
				 * This function calculates the value of gamma_1 in equation 19.
				 */
				inline distribn_t gamma1(const LatticePosition tanVec,
										 const distribn_t tanComp,
										 const distribn_t tau,
										 const distribn_t *const fNeq)
				{
					return 2.0 * tanComp * stress(tanVec, tau, fNeq);
				}

				/**
				 * This function calculates the value of gamma_2 in equation 19.
				 */
				inline distribn_t gamma2(const LatticePosition cp,
										 const distribn_t tau,
										 const distribn_t *const fNeq)
				{
					return 2.0 * stress(cp, tau, fNeq);
				}

				/**
				 * This function calculates the square bracket term in equation 18.
				 */
				inline LatticePosition momentumCorrection(const LatticePosition cp,
														  const LatticePosition normal,
														  const distribn_t tau,
														  const distribn_t *const fNeq)
				{
					LatticePosition tanVec = tangentialVector(cp, normal);
					distribn_t tanComp = cp.Dot(tanVec);
					distribn_t normalComp = cp.Dot(normal);
					distribn_t g1 = gamma1(tanVec, tanComp, tau, fNeq);
					distribn_t g2 = gamma2(cp, tau, fNeq);
					return tanVec * g1 + normal * (g2 - tanComp * g1) / normalComp;
				}

				// Attributes of this class
				CollisionType &collider;
				const geometry::neighbouring::NeighbouringLatticeData &neighbouringLatticeData;
				iolets::BoundaryValues &iolet;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_YANGPRESSUREDELEGATE_H
