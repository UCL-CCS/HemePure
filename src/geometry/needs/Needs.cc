
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "geometry/needs/Needs.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
	namespace geometry
	{
		Needs::Needs(const site_t blockCount,
				const std::unordered_set<site_t>& readBlock,
				const proc_t readingGroupSize,
				const proc_t readingGroupSpacing,
				net::InterfaceDelegationNet & net) :
			communicator(net.GetCommunicator()),
			readingGroupSize(readingGroupSize), readingGroupSpacing(readingGroupSpacing)
		{

			// Map containing list of reading cores
			proc_t i = 1; proc_t readingCore = 0;
			while (i <= readingGroupSize) {
				readingCores.push_back(readingCore);
				readingCore += readingGroupSpacing; i++;
			}

			// Check if core will read
			reader = false;
			for (std::vector<proc_t>::iterator it = readingCores.begin(); it != readingCores.end(); ++it)
				if (*it == communicator.Rank())
					reader = true;

			// Compile the blocks needed here into an array of indices, instead of an array of bools
			std::vector<std::vector<site_t> > blocksNeededHere(readingGroupSize);
			for (site_t block = 0; block < blockCount; ++block)
			{
				if (readBlock.find(block) != readBlock.end())
				{
					blocksNeededHere[GetReadingCoreForBlock(block)].push_back(block);
				}
			}

			// Share the counts of needed blocks
			int blocksNeededSize[readingGroupSize];
			std::vector<int> blocksNeededSizes(communicator.Size());
			for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
			{
				blocksNeededSize[readingCore] = blocksNeededHere[readingCore].size();
				net.RequestGatherSend(blocksNeededSize[readingCore], readingCores[readingCore]);
			}
			if (reader)
			{
				net.RequestGatherReceive(blocksNeededSizes);
			}
			net.Dispatch();

			// Communicate the arrays of needed blocks
			std::vector<site_t> blocksNeededOn;
			for (proc_t readingCore = 0; readingCore < readingGroupSize; readingCore++)
			{
				net.RequestGatherVSend(blocksNeededHere[readingCore], readingCores[readingCore]);
			}
			if (reader)
			{
				net.RequestGatherVReceive(blocksNeededOn, blocksNeededSizes);
			}
			net.Dispatch();

			coresWantingBlocksBuffer[-1].push_back(-1);
			if (reader)
			{
				int needsPassed = 0;
				// Transpose the blocks needed on cores matrix
				for (proc_t sendingCore = 0; sendingCore < communicator.Size(); sendingCore++)
				{
					for (int needForThisSendingCore = 0; needForThisSendingCore < blocksNeededSizes[sendingCore];
							++needForThisSendingCore)
					{
						coresWantingBlocksBuffer[blocksNeededOn[needsPassed]].push_back(sendingCore);

						++needsPassed;
					}
				} //for sendingCore
			} //if a reading core

			//{
			//	int i = 0;
			//	char hostname[256];
			//	gethostname(hostname, sizeof(hostname));
			//	printf("PID %d (%d) on %s ready for attach\n", getpid(), communicator.Rank(), hostname);
			//	fflush(stdout);
			//	while (0 == i)
			//		sleep(5);
			//}
		}

		proc_t Needs::GetReadingCoreForBlock(const site_t blockNumber) const
		{
			return proc_t(blockNumber % readingGroupSize);
		}
	} //namespace
} //namespace
