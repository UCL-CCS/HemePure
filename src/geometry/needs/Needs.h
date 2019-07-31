
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_GEOMETRY_NEEDS_NEEDS_H
#define HEMELB_GEOMETRY_NEEDS_NEEDS_H
#include <vector>
#include "net/net.h"
#include "net/IOCommunicator.h"

// PURE
#include <unordered_set>
#include <unordered_map>

namespace hemelb
{
	namespace geometry
	{
		/*
		   JH Note: Class created to make the new needs communication strategy in #142 testable.
		 */
		/***
		 *  Class defining HemeLB needs communication
		 Used by geometry reader to know where to send which blocks.
		 */
		class Needs
		{
			public:
				/***
				 * Constructor for Needs manager.
				 * @param BlockCount Count of blocks
				 * @param readBlock Which cores need which blocks, as an array of booleans.
				 * @param readingGroupSize Number sof cores to use for reading blocks
				 * @param net Instance of Net communication class to use.
				 */
				Needs(const site_t blockCount,
						const std::unordered_set<site_t>& readBlock,
						const proc_t readingGroupSize,
						const proc_t readingGroupSpacing,
						net::InterfaceDelegationNet &net);

				/***
				 * Which processors need a given block?
				 * @param block Block number to query
				 * @return Vector of ranks in the decomposition topology which need this block
				 */
				const std::vector<proc_t> & CoreNeedingBlock(const site_t &block) const
				{
					if (coresWantingBlocksBuffer.find(block) != coresWantingBlocksBuffer.end())
						 return coresWantingBlocksBuffer.at(block);
					else return coresWantingBlocksBuffer.at(-1);
				}

				void CoreNeedingBlock_erase(const site_t &block)
				{
					coresWantingBlocksBuffer.erase(block);
				}

				const proc_t CoreReadingBlock(const proc_t readingCore) const
				{
					return readingCores[readingCore];
				}

				const bool Reader() const
				{
					return reader;
				}

				/***
				 * Which core should be responsible for reading a given block? This core does not necessarily
				 * require information about the block
				 *
				 * @param blockNumber Block number to query
				 * @return Rank in the decomposition topology, for core which should read the block.
				 */
				proc_t GetReadingCoreForBlock(const site_t blockNumber) const;
			private:
				bool reader;
				const net::MpiCommunicator & communicator;
				const proc_t readingGroupSize;
				const proc_t readingGroupSpacing;
				std::unordered_map<site_t, std::vector<proc_t> > coresWantingBlocksBuffer;
				std::vector<proc_t> readingCores;
		};
	}
}
#endif // HEMELB_GEOMETRY_NEEDS_NEEDS_H
