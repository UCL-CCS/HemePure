
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
			 * 2. All solid nodes beyond the iolet have two neighbouring fluid sites
			 *    in a certain lattice direction.
			 * 3. The collision kernel adopts a single relaxation time.
			 * 4. There is only one inlet and one outlet in each rank. (To be relaxed.)
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

					std::cout << "rank " << initParams.latDat->GetLocalRank() << " siteCount " << initParams.siteCount << " localIoletCount " << initParams.boundaryObject->GetLocalIoletCount() << " ioletType " << initParams.boundaryObject->GetIoletType() << std::endl;

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

							/**
							 * Sort the lattice directions in an ascending order based on their allignment to
							 * the iolet normal vector. Due to assumption 1, these directions are the same for
							 * all sites linked to the same iolet. Therefore, only the first site is concerned.
							 * (Currently assumption 4 is applied, but it will be relaxed.)
							 */
							if (localIndex == rangeIt->first)
							{
								SortDirectionsAlignToIOletNormal(localSite.GetIoletId());
								dirG = dirs[0]; // the direction described by the mapping G
								cg = LatticeVector(LatticeType::CX[dirG], LatticeType::CY[dirG], LatticeType::CZ[dirG]);
							}

							// Ensure the data of the two fluid sites is available.
							for (Direction i = 1; i < LatticeType::NUMVECTORS; ++i)
							{
								/**
								 * Taking one step from the current site in the current direction, i, gives an
								 * outer-neighbour (solid) node. Taking a further step in the direction mapped
								 * by P gives the first fluid site. Taking one more step in this direction
								 * gives the second fluid site.
								 */

								// Skip if there is no iolet in this direction
								if (!localSite.HasIolet(i)) continue;

								const LatticeVector ci = LatticeVector(LatticeType::CX[i],
																	   LatticeType::CY[i],
																	   LatticeType::CZ[i]);
								const LatticeVector solidSiteLocation = localSiteLocation + ci;

								Direction j = 0, dirP = dirs[0];
								do
								{
									const LatticeVector cp = LatticeVector(LatticeType::CX[dirP],
																		   LatticeType::CY[dirP],
																		   LatticeType::CZ[dirP]);

									const LatticeVector firstFluidSiteLocation = solidSiteLocation + cp;
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
											// Request data from another rank
											initParams.neighbouringDataManager->RegisterNeededSite(
												initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(firstFluidSiteLocation));
										}
										if (secondFluidSiteHomeProc != initParams.latDat->GetLocalRank())
										{
											// Request data from another rank
											initParams.neighbouringDataManager->RegisterNeededSite(
												initParams.latDat->GetGlobalNoncontiguousSiteIdFromGlobalCoords(secondFluidSiteLocation));
										}
										//printf("site [%ld %ld %ld], dirG %u, dirP %u\n", localSiteLocation.x, localSiteLocation.y, localSiteLocation.z, dirG, dirP);
										break;
									}

									j++;
									dirP = dirs[j];
								}
								while (dirP != 0);

								// Since dirs are sorted, having dirP == 0 implies the mapping P is empty,
								// violating assumption 2. In this case, a finer grid should be used.
								if (dirP == 0)
								{
									hemelb::log::Logger::Log<hemelb::log::Error, hemelb::log::OnePerCore>(
										"A higher resolution is required for site [%ld, %ld, %ld] of index %ld",
										localSiteLocation.x, localSiteLocation.y, localSiteLocation.z, localIndex);
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
					if (!site.HasIolet(direction)) return;

					int boundaryId = site.GetIoletId();
					LatticePosition ioletNormal = iolet.GetIolets()[boundaryId]->GetNormal();
					distribn_t fNew;

					// Obtain cp, wallDistance, and the old distributions at the first and second fluid nodes.
					// Here wall is referred to as the iolet plane.
					LatticeVector cp;			  // the lattice vector given by the mapping P
					LatticeDistance wallDistance; // distance between the wall and the first fluid site
					const distribn_t *firstFluidFOld, *secondFluidFOld;
					GetFluidSitesData(site, direction, latticeData, cp, wallDistance, firstFluidFOld, secondFluidFOld);
					//printf("site [%ld %ld %ld], dir %u, cp [%ld %ld %ld], wallDistance %lf, cg [%ld %ld %ld]\n",
					//	site.GetGlobalSiteCoords().x, site.GetGlobalSiteCoords().y, site.GetGlobalSiteCoords().z,
					//	direction, cp.x, cp.y, cp.z, wallDistance, cg.x, cg.y, cg.z);

					// Calculate the densities, momenta, and equilibrium distributions and store them in HydroVars
					kernels::HydroVars<typename CollisionType::CKernel> hVfirstFluid(firstFluidFOld);
					collider.kernel.CalculateDensityMomentumFeq(hVfirstFluid, 0); // the second argument is dummy
					kernels::HydroVars<typename CollisionType::CKernel> hVsecondFluid(secondFluidFOld);
					collider.kernel.CalculateDensityMomentumFeq(hVsecondFluid, 0); // the second argument is dummy

					// Calculate the non-equilibrium distributions on the wall (equation 16)
					distribn_t fNeqWall[LatticeType::NUMVECTORS];
					for (Direction i = 0; i < LatticeType::NUMVECTORS; ++i)
					{
						fNeqWall[i] = (1.0 + wallDistance) * hVfirstFluid.GetFNeq().f[i]
									  - wallDistance * hVsecondFluid.GetFNeq().f[i];
					}

					Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];
					if (unstreamed == dirG && cg == cp)
					{
						// Calculate the density on the wall (equation 13).
						// The equilibrium density, 1, is absorbed in the iolet pressure.
						// Also note that h = 1 in lattice units.
						iolets::InOutLetCosine *cosIolet =
							dynamic_cast<iolets::InOutLetCosine *>(iolet.GetIolets()[boundaryId]);
						const distribn_t visc = Cs2 * (hydroVars.tau - 0.5); // kinematic viscosity
						LatticeDensity densityWall = 3.0 * (cosIolet->GetPressure(iolet.GetTimeStep()) + visc *
															stress(-ioletNormal, hydroVars.tau, fNeqWall));

						// Interpolate the density at the first fluid site accounting for the BC (equation 12)
						LatticeDensity densityBC = (densityWall + wallDistance * hVsecondFluid.density) / (1.0 + wallDistance);

						// Calculate the post-collision distributions at the second fluid site
						collider.Collide(lbmParams, hVsecondFluid);

						// Calculate the distribution of the unstreamed direction (equation 11)
						fNew = site.GetFOld<LatticeType>()[unstreamed]
								+ site.GetFOld<LatticeType>()[direction]
								+ 2.0 * LatticeType::EQMWEIGHTS[unstreamed] * (densityBC - hydroVars.density)
								- hVsecondFluid.GetFPostCollision()[direction];
					}
					else
					{
						// Construct a HydroVars structure at the outer-neighbour node
						distribn_t fOuterWall[LatticeType::NUMVECTORS];
						kernels::HydroVars<typename CollisionType::CKernel> hVouterWall(fOuterWall);

						// Extrapolate the density and momentum at the outer-neighbour node (equation 18)
						hVouterWall.density = 2.0 * hVfirstFluid.density - hVsecondFluid.density;
						hVouterWall.momentum = (hVfirstFluid.momentum * 4.0 * wallDistance +
												hVsecondFluid.momentum * (1.0 - 2.0 * wallDistance)
												- momentumCorrection(cp, -ioletNormal, hydroVars.tau, fNeqWall)
											   ) / (1.0 + 2.0 * wallDistance);

						// Calculate the equilibrium distributions
						LatticeType::CalculateFeq(hVouterWall.density,
												  hVouterWall.momentum.x,
												  hVouterWall.momentum.y,
												  hVouterWall.momentum.z,
												  hVouterWall.GetFEqPtr());

						// Extrapolate the non-equilibrium distribution of the unstreamed direction (equation 20)
						distribn_t fNeqOuterWall = 2.0 * hVfirstFluid.GetFNeq().f[unstreamed]
												   - hVsecondFluid.GetFNeq().f[unstreamed];

						// Calculate the distribution of the unstreamed direction (equation 17).
						// Assumption 3 is applied here: A equals 1/tau times the identity matrix.
						fNew = hVouterWall.GetFEq().f[unstreamed] + fNeqOuterWall * (1.0 - 1.0 / hydroVars.tau);
					}


					// Below is the Nash BC

					// Set the density at the "ghost" site to be the density of the iolet.
					distribn_t ghostDensity = iolet.GetBoundaryDensity(boundaryId);

					// Calculate the velocity at the ghost site, as the component normal to the iolet.
					// util::Vector3D<float> ioletNormal = iolet.GetIolets()[boundaryId]->GetNormal();

					// Note that the division by density compensates for the fact that v_x etc have momentum
					// not velocity.
					distribn_t component = (hydroVars.momentum / hydroVars.density).Dot(ioletNormal);

					// TODO it's ugly that we have to do this.
					// TODO having to give 0 as an argument is also ugly.
					// TODO it's ugly that we have to give hydroVars a nonsense distribution vector
					// that doesn't get used.
					kernels::HydroVars<typename CollisionType::CKernel> ghostHydrovars(site);

					ghostHydrovars.density = ghostDensity;
					ghostHydrovars.momentum = ioletNormal * component * ghostDensity;

					collider.kernel.CalculateFeq(ghostHydrovars, 0);

					// Direction unstreamed = LatticeType::INVERSEDIRECTIONS[direction];

					distribn_t fOld = ghostHydrovars.GetFEq()[unstreamed];
					//*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed) = ghostHydrovars.GetFEq()[unstreamed];
					if (unstreamed == dirG && cg == cp)
					{
						//printf("unstreamed [%ld %ld %ld]\n", site.GetGlobalSiteCoords().x, site.GetGlobalSiteCoords().y, site.GetGlobalSiteCoords().z);
						*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed) = fNew;
					}
					else
					{
						*latticeData->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS + unstreamed) = fNew;
					}

					if (iolet.GetTimeStep() % 50 == 0)
					{
						//printf("direction %u, Nash %.15lf, Yang %.15lf, dif %.5lf\n", direction, fOld, fNew, (fNew - fOld)/fOld);
					}
				}

			protected:
				/**
				 * Sort the lattice directions in an ascending order based on their allignment
				 * to the iolet normal vector. See equation A.1 in the paper.
				 */
				inline void SortDirectionsAlignToIOletNormal(const int boundaryId)
				{
        			// Assumption 1 is implied here
					LatticePosition ioletNormal = iolet.GetIolets()[boundaryId]->GetNormal();

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
					// Sort dirs by comparing any two elements of dist
        			std::sort(dirs.begin(), dirs.end(), [&dist](Direction i, Direction j) {return dist[i] < dist[j];});
					for (Direction k = 0; k < LatticeType::NUMVECTORS; ++k)
					{
						printf("%u %lf ", dirs[k], dist[dirs[k]]);
					}
					printf("\n");
      			}

				/**
				 * Obtain the data at the two fluid sites. The logic used in the constructor is applied.
				 */
				inline void GetFluidSitesData(const geometry::Site<geometry::LatticeData> &site,
											  const Direction &i,
											  geometry::LatticeData *const latDat,
											  LatticeVector &cp,
											  LatticeDistance &firstFluidWallDistance,
											  const distribn_t* &firstFluidFOld,
											  const distribn_t* &secondFluidFOld)
				{
					// The parameter i is assumed to be the stream direction
					const LatticeVector ci = LatticeVector(LatticeType::CX[i],
														   LatticeType::CY[i],
														   LatticeType::CZ[i]);
					const LatticeVector solidSiteLocation = site.GetGlobalSiteCoords() + ci;

					Direction j = 0, dirP = dirs[0];
					do
					{
						cp = LatticeVector(LatticeType::CX[dirP], LatticeType::CY[dirP], LatticeType::CZ[dirP]);

						const LatticeVector firstFluidSiteLocation = solidSiteLocation + cp;
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
						dirP = dirs[j];
					}
					while (dirP != 0);
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
				 * Assumption 3 is applied here: the collision matrix A equals 1/tau
				 * times the identity matrix, where tau is the single relaxation time.
				 */
				inline distribn_t stress(const LatticePosition vec,
										 const distribn_t tau,
										 const distribn_t *const fNeq)
				{
					distribn_t stress = 0;
					for (Direction k = 0; k < LatticeType::NUMVECTORS; ++k)
					{
						stress += B_k(k, vec) * fNeq[k];
					}
					// Note that A_k = A_k* = 1/tau and h = 1 in lattice units
					stress *= -1.0 / tau;
					return stress;
				}

				inline LatticePosition tangentialVector(const LatticePosition vec, const LatticePosition normal)
				{
					distribn_t normalComp = vec.Dot(normal);
					return (vec - vec * normalComp).Normalise();
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
				std::array<Direction, LatticeType::NUMVECTORS> dirs;
				Direction dirG;
				LatticeVector cg;
			};
		}
	}
}

#endif // HEMELB_LB_STREAMERS_YANGPRESSUREDELEGATE_H
