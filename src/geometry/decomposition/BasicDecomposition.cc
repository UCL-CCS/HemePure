
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/decomposition/BasicDecomposition.h"
#include "net/mpi.h"

namespace hemelb
{
	namespace geometry
	{
		namespace decomposition
		{
			BasicDecomposition::BasicDecomposition(const Geometry& geometry,
					const lb::lattices::LatticeInfo& latticeInfo,
					const net::MpiCommunicator& communicator,
					const std::unordered_map<site_t, std::pair<uint16_t, uint16_t> >& blockInformation,
					const std::unordered_map<site_t, uint16_t>& blockWeights) :
				geometry(geometry), latticeInfo(latticeInfo), communicator(communicator),
				blockInformation(blockInformation), blockWeights(blockWeights)
			{
			}

			// Very dumb decomposition. For testing only!
			void BasicDecomposition::DecomposeDumb(
					std::unordered_map<site_t, proc_t>& procAssignedToEachBlock,
					sitedata_t nonEmptyBlocks)
			{
				site_t solidBlock = 0;
				for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
				{
					if (blockInformation.find(block) != blockInformation.end()) {
						procAssignedToEachBlock[block] = solidBlock/
							(double) nonEmptyBlocks*communicator.Size();
						solidBlock++;
					}
				}
			}

			void BasicDecomposition::Decompose(std::unordered_map<site_t, proc_t>& procAssignedToEachBlock)
			{
				// Keep a count of the number of non-empty blocks that haven't yet been assigned
				// a processor.
				site_t unvisitedFluidBlockCount = 0;
				for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
				{
					if (blockInformation.find(block) != blockInformation.end())
					{
#ifdef HEMELB_USE_GMYPLUS
						unvisitedFluidBlockCount += blockWeights.at(block);
#else
						unvisitedFluidBlockCount++;
#endif
					}
				}

				DivideBlocks(procAssignedToEachBlock,
						unvisitedFluidBlockCount,
						geometry,
						communicator.Size(),
						blockInformation);
			}

			//void BasicDecomposition::Validate(std::vector<proc_t>& procAssignedToEachBlock)
			//{
			//	log::Logger::Log<log::Debug, log::OnePerCore>("Validating procForEachBlock");

			//	std::vector<proc_t> procForEachBlockRecv = communicator.AllReduce(procAssignedToEachBlock, MPI_MAX);

			//	for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
			//	{
			//		if (procAssignedToEachBlock[block] != procForEachBlockRecv[block])
			//		{
			//			log::Logger::Log<log::Critical, log::OnePerCore>("At least one other proc thought block %li should be on proc %li but we locally had it as %li",
			//					block,
			//					procAssignedToEachBlock[block],
			//					procForEachBlockRecv[block]);
			//		}
			//	}
			//}

			void BasicDecomposition::DivideBlocks(std::unordered_map<site_t, proc_t>& unitForEachBlock,
					site_t unassignedBlocks,
					const Geometry& geometry,
					const proc_t unitCount,
					const std::unordered_map<site_t, std::pair<uint16_t, uint16_t> >& blockInformation)
			{
				// Initialise the unit being assigned to, and the approximate number of blocks
				// required on each unit.
				proc_t currentUnit = 0;

				site_t targetBlocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (communicator.Size()));

				// Create a set to monitor whether each block has been assigned yet.
				std::unordered_set<site_t> blockAssigned;

				// Create lists of the current edge of blocks on the current proc and the edge being expanded into
				std::vector<BlockLocation> currentEdge;
				std::vector<BlockLocation> expandedEdge;

				// Domain Decomposition. Pick a site. Set it to the rank we are
				// looking at. Find its neighbours and put those on the same
				// rank, then find the next-nearest neighbours, etc. until we
				// have a completely joined region, or there are enough fluid
				// sites on the rank.  In the former case, start again at
				// another site. In the latter case, move on to the next rank.
				// Do this until all sites are assigned to a rank. There is a
				// high chance of of all sites on a rank being joined.

				site_t blockNumber = -1;
				site_t blocksOnCurrentProc = 0;

				// Iterate over all blocks.
				for (site_t blockCoordI = 0; blockCoordI < geometry.GetBlockDimensions().x; blockCoordI++)
				{
					for (site_t blockCoordJ = 0; blockCoordJ < geometry.GetBlockDimensions().y; blockCoordJ++)
					{
						for (site_t blockCoordK = 0; blockCoordK < geometry.GetBlockDimensions().z; blockCoordK++)
						{
							// Block number is the number of the block we're currently on.
							blockNumber++;

							if (blockInformation.find(blockNumber) == blockInformation.end())
							{
								continue;
							}

							// Alternatively, if this block has already been assigned, move on.
							if (blockAssigned.find(blockNumber) != blockAssigned.end())
							{
								continue;
							}

							// Assign this block to the current unit.
							blockAssigned.insert(blockNumber);
							unitForEachBlock[blockNumber] = currentUnit;

#ifdef HEMELB_USE_GMYPLUS
							blocksOnCurrentProc += blockWeights.at(blockNumber);
#else
							blocksOnCurrentProc++;
#endif

							// Record the location of this initial block.
							currentEdge.clear();
							BlockLocation lNew(blockCoordI, blockCoordJ, blockCoordK);
							currentEdge.push_back(lNew);

							// The subdomain can grow.
							bool isRegionGrowing = true;

							int blockAssignedCount = blockAssigned.size();

							// While the region can grow (i.e. it is not bounded by solids or visited
							// sites), and we need more sites on this particular rank.
							while (blocksOnCurrentProc < targetBlocksPerUnit && isRegionGrowing)
							{
								expandedEdge.clear();

								// Sites added to the edge of the mClusters during the iteration.
								isRegionGrowing = Expand(expandedEdge,
										blockAssigned,
										unitForEachBlock,
										blocksOnCurrentProc,
										currentEdge,
										currentUnit,
										targetBlocksPerUnit);

								// When the new layer of edge sites has been found, swap the buffers for
								// the current and new layers of edge sites.
								currentEdge.swap(expandedEdge);
							}

							//// If we tried to expand and did not get far.
							//// Prevents creation of small isolated regions.
							//if (((blockAssigned.size() - blockAssignedCount) < 15) && (currentUnit > 0))
							//{
							//	for (std::unordered_map<site_t, proc_t>::const_iterator iter = unitForEachBlock.begin();
							//			iter != unitForEachBlock.end(); ++iter)
							//	{
							//		if (iter->second == currentUnit)
							//			// Allocate these blocks to previous unit.
							//			unitForEachBlock.at(iter->first) = currentUnit - 1;
							//	}
							//}

							// If we have enough sites, we have finished.
							if (blocksOnCurrentProc >= targetBlocksPerUnit)
							{
								++currentUnit;

								unassignedBlocks -= blocksOnCurrentProc;
								targetBlocksPerUnit = (site_t) ceil((double) unassignedBlocks / (double) (unitCount - currentUnit));

								blocksOnCurrentProc = 0;
							}
							// If not, we have to start growing a different region for the same rank:
							// region expansions could get trapped.
						} // Block co-ord k
					} // Block co-ord j
				} // Block co-ord i

				//blockNumber = -1;
				//// Check which ranks own the neighbouring blocks of each block.
				//for (site_t blockCoordI = 0; blockCoordI < geometry.GetBlockDimensions().x; blockCoordI++)
				//{
				//	for (site_t blockCoordJ = 0; blockCoordJ < geometry.GetBlockDimensions().y; blockCoordJ++)
				//	{
				//		for (site_t blockCoordK = 0; blockCoordK < geometry.GetBlockDimensions().z; blockCoordK++)
				//		{
				//			// Block number is the number of the block we're currently on.
				//			if (unitForEachBlock.find(++blockNumber) == unitForEachBlock.end())
				//			{
				//				continue;
				//			}

				//			// Create map to record rank of neighbouring blocks.
				//			std::unordered_map<proc_t, int> KFC;
				//			for (site_t neighI = util::NumericalFunctions::max<site_t>(0, blockCoordI - 1); (neighI
				//						<= (blockCoordI + 1)) && (neighI < geometry.GetBlockDimensions().x); ++neighI)
				//			{
				//				for (site_t neighJ = util::NumericalFunctions::max<site_t>(0, blockCoordJ - 1); (neighJ
				//							<= (blockCoordJ + 1)) && (neighJ < geometry.GetBlockDimensions().y); ++neighJ)
				//				{
				//					for (site_t neighK = util::NumericalFunctions::max<site_t>(0, blockCoordK - 1); (neighK
				//								<= (blockCoordK + 1)) && (neighK < geometry.GetBlockDimensions().z); ++neighK)
				//					{
				//						site_t neighBlockId = geometry.GetBlockIdFromBlockCoordinates(neighI, neighJ, neighK);
				//						if (unitForEachBlock.find(neighBlockId) != unitForEachBlock.end())
				//						{
				//							proc_t neighRank = unitForEachBlock.at(neighBlockId);
				//							if (KFC.find(neighRank) == KFC.end())
				//								KFC[neighRank] = 0;
				//							KFC.at(neighRank)++;
				//						}
				//					}
				//				}
				//			}

				//			// If no neighbouring blocks are owned by the rank owning the source block...
				//			// ... or if the number of blocks owned by rank owning source block is < 8.
				//			if (
				//					KFC.find(unitForEachBlock.at(blockNumber)) == KFC.end() ||
				//					KFC.at(  unitForEachBlock.at(blockNumber)) < 8)
				//			{
				//				int maxValue = 0; proc_t maxValueProc = unitForEachBlock.at(blockNumber);
				//				for (std::unordered_map<proc_t, int>::const_iterator iter = KFC.begin();
				//						iter != KFC.end(); ++iter)
				//				{
				//					if (iter->second > maxValue)
				//					{
				//						maxValue     = iter->second;
				//						maxValueProc = iter->first;
				//					}
				//				}
				//				unitForEachBlock.at(blockNumber) = maxValueProc;
				//			}
				//		}
				//	}
				//}

				// To measure quality of distribution.
				std::vector<sitedata_t> totalBlockWeights(communicator.Size(), 0);

				// Iterate over all blocks (again).
				for (site_t blockNumber = 0; blockNumber < geometry.GetBlockCount(); ++blockNumber)
				{
					// Weight of all blocks on partition.
					if (unitForEachBlock.find(blockNumber) != unitForEachBlock.end())
#ifdef HEMELB_USE_GMYPLUS
						totalBlockWeights[unitForEachBlock.at(blockNumber)] += blockWeights.at(blockNumber);
#else
						totalBlockWeights[unitForEachBlock.at(blockNumber)]++;
#endif
				}

				sitedata_t max = *std::max_element(
						totalBlockWeights.begin(), totalBlockWeights.end());
				sitedata_t min = *std::min_element(
						totalBlockWeights.begin(), totalBlockWeights.end());
				double average = std::accumulate(
						totalBlockWeights.begin(), totalBlockWeights.end(), 0)/(double)totalBlockWeights.size();

				double d_min = (double)min / average;
				double d_max = (double)max / average;

				double d_ratio = (d_max - d_min) / (d_max + d_min);

				communicator.Barrier();
				if (communicator.Rank() == 0)
					log::Logger::Log<log::Info, log::OnePerCore>("----> load distribution: %f", d_ratio);//totalBlockWeights[communicator.Rank()]/average);
			}

			bool BasicDecomposition::Expand(std::vector<BlockLocation>& expansionBlocks,
					std::unordered_set<site_t>& blockAssigned,
					std::unordered_map<site_t, proc_t>& unitForEachBlock,
					site_t &blocksOnCurrentUnit,
					const std::vector<BlockLocation>& edgeBlocks,
					const proc_t currentUnit,
					const site_t blocksPerUnit)
			{
				bool regionExpanded = false;

				// For sites on the edge of the domain (sites_a), deal with the neighbours.
				for (unsigned int edgeBlockId = 0; (edgeBlockId < edgeBlocks.size()) && (blocksOnCurrentUnit < blocksPerUnit); edgeBlockId++)
				{
					const BlockLocation& edgeBlockCoords = edgeBlocks[edgeBlockId];

					for (Direction direction = 1; direction < latticeInfo.GetNumVectors() && blocksOnCurrentUnit < blocksPerUnit; direction++)
					{
						// Record neighbour location.
						BlockLocation neighbourCoords = edgeBlockCoords + latticeInfo.GetVector(direction);

						//{
						//	int i = 0;
						//	char hostname[256];
						//	gethostname(hostname, sizeof(hostname));
						//	printf("PID %d (%d) on %s ready for attach\n", getpid(), communicator.Rank(), hostname);
						//	fflush(stdout);
						//	while (0 == i)
						//		sleep(5);
						//}

						// Move on if neighbour is outside the bounding box.
						if (!geometry.AreBlockCoordinatesValid(neighbourCoords))
						{
							continue;
						}

						site_t neighBlockId = geometry.GetBlockIdFromBlockCoordinates(neighbourCoords.x,
								neighbourCoords.y,
								neighbourCoords.z);

						// Don't use this block if it has no fluid sites, or if it has already been assigned to a processor.
						if (blockInformation.find(neighBlockId) == blockInformation.end() ||
								blockAssigned.find(neighBlockId) != blockAssigned.end())
						{
							continue;
						}

						// Set the rank for a neighbour and update the fluid site counters.
						blockAssigned.insert(neighBlockId);
						unitForEachBlock[neighBlockId] = currentUnit;
#ifdef HEMELB_USE_GMYPLUS
						blocksOnCurrentUnit += blockWeights.at(neighBlockId);
#else
						blocksOnCurrentUnit++;
#endif

						// Neighbour was found, so the region can grow.
						regionExpanded = true;

						// Record the location of the neighbour.
						expansionBlocks.push_back(neighbourCoords);
					}
				}
				return regionExpanded;
			}
		} /* namespace decomposition */
	} /* namespace geometry */
} /* namespace hemelb */
