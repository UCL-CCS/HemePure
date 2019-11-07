
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <cmath>
#include <list>
#include <map>
#include <algorithm>
#include <zlib.h>

#include <iostream>
#include <fstream>

#include "io/formats/geometry.h"
#include "io/writers/xdr/XdrMemReader.h"
#include "geometry/decomposition/BasicDecomposition.h"
#include "geometry/decomposition/OptimisedDecomposition.h"
#include "geometry/GeometryReader.h"
#include "lb/lattices/D3Q27.h"
#include "net/net.h"
#include "net/IOCommunicator.h"
#include "log/Logger.h"
#include "util/utilityFunctions.h"
#include "constants.h"

namespace hemelb
{
	namespace geometry
	{

		GeometryReader::GeometryReader(const lb::lattices::LatticeInfo& latticeInfo,
				reporting::Timers &atimings, const net::IOCommunicator& ioComm) :
			latticeInfo(latticeInfo), hemeLbComms(ioComm), timings(atimings)
		{
			// This rank should participate in the domain decomposition if
			//  - we're not on core 0 (the only core that might ever not participate)
			//  - there's only one processor (so core 0 has to participate)

			// Create our own group, without the root node if we're not running with it.
			if (ioComm.Size() > 1)
			{
				participateInTopology = !ioComm.OnIORank();

				std::vector<int> lExclusions(1);
				lExclusions[0] = 0;
				net::MpiGroup computeGroup = hemeLbComms.Group().Exclude(lExclusions);
				// Create a communicator just for the domain decomposition.
				computeComms = ioComm.Create(computeGroup);
			}
			else
			{
				participateInTopology = true;
				computeComms = ioComm;
			}
		}

		GeometryReader::~GeometryReader()
		{
		}

		Geometry GeometryReader::LoadAndDecompose(const std::string& dataFilePath)
		{
			timings[hemelb::reporting::Timers::fileRead].Start();

			// Create hints about how we'll read the file. See Chapter 13, page 400 of the MPI 2.2 spec.
			MPI_Info fileInfo;
			HEMELB_MPI_CALL(MPI_Info_create, (&fileInfo));
			std::string accessStyle = "access_style";
			std::string accessStyleValue = "sequential";
			std::string buffering = "collective_buffering";
			std::string bufferingValue = "true";

			HEMELB_MPI_CALL(MPI_Info_set, (fileInfo,
						const_cast<char*> (accessStyle.c_str()),
						const_cast<char*> (accessStyleValue.c_str()))
					);
			HEMELB_MPI_CALL(MPI_Info_set, (fileInfo,
						const_cast<char*> (buffering.c_str()),
						const_cast<char*> (bufferingValue.c_str()))
					);

			// Open the file.
			file = net::MpiFile::Open(hemeLbComms, dataFilePath, MPI_MODE_RDONLY, fileInfo);

			log::Logger::Log<log::Info, log::Singleton>("----> opened data file %s", dataFilePath.c_str());
			// TODO: Why is there this fflush?
			fflush(NULL);

			// Set the view to the file.
			file.SetView(0, MPI_CHAR, MPI_CHAR, "native", fileInfo);

			// The processor assigned to each block.
			std::unordered_map<site_t, proc_t>* principalProcForEachBlock =
				new std::unordered_map<site_t, proc_t>();

#ifdef HEMELB_USE_GMYPLUS
			log::Logger::Log<log::Info, log::Singleton>("----> accepting *.gmy+");
#endif
#ifdef HEMELB_USE_MPI_CALL
			log::Logger::Log<log::Info, log::Singleton>("----> using standard MPI calls in ReadInBlock()");
#endif

			log::Logger::Log<log::Info, log::Singleton>("----> reading preamble");
			Geometry geometry = ReadPreamble();

			log::Logger::Log<log::Info, log::Singleton>("----> reading header (start)");
			ReadHeader(geometry.GetBlockCount());
			log::Logger::Log<log::Info, log::Singleton>("----> reading header (end)");

			// Close the file - only the ranks participating in the topology need to read it again.
			file.Close();

			sitedata_t siteCount = 0;
			for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
			{
				if (fluidSitesOnEachBlock.find(block) != fluidSitesOnEachBlock.end())
					siteCount += fluidSitesOnEachBlock.at(block);
			}

#ifndef HEMELB_USE_PARMETIS
			fluidSitesOnEachBlock.clear();
#endif

			log::Logger::Log<log::Info, log::Singleton>(
					"----> non-empty blocks: %lu",
					nonEmptyBlocks);
			log::Logger::Log<log::Info, log::Singleton>(
					"          total blocks: %lu",
					geometry.GetBlockCount());
			log::Logger::Log<log::Info, log::Singleton>(
					"                 ratio: %f",
					nonEmptyBlocks/(double)geometry.GetBlockCount());
			log::Logger::Log<log::Info, log::Singleton>(
					"                 sites: %lu",
					siteCount);

			log::Logger::Log<log::Info, log::Singleton>(
					"----> blockInformation.size(): %lu",
					blockInformation.size());
			log::Logger::Log<log::Info, log::Singleton>(
					" fluidSitesOnEachBlock.size(): %lu",
					fluidSitesOnEachBlock.size());
			log::Logger::Log<log::Info, log::Singleton>(
					"          blockWeights.size(): %lu",
					blockWeights.size());

			log::Logger::Log<log::Info, log::Singleton>("----> is blockInformation.size() == nonEmptyBlocks? %s", blockInformation.size() == nonEmptyBlocks ? "yes" : "no");

#ifndef HEMELB_USE_PARMETIS
			log::Logger::Log<log::Info, log::Singleton>("----> not optimising decomposition");
#endif

			log::Logger::Log<log::Info, log::Singleton>("----> basic decomposition (start)");
			timings[hemelb::reporting::Timers::initialDecomposition].Start();
			if (participateInTopology)
			{
				// Get an initial base-level decomposition of the domain macro-blocks over processors.
				// This will later be improved upon by ParMETIS.
				decomposition::BasicDecomposition basicDecomposer(geometry,
						latticeInfo,
						computeComms,
						blockInformation,
						blockWeights);
		//		basicDecomposer.DecomposeBlock(
		//				*principalProcForEachBlock);
		//		basicDecomposer.DecomposeDumb(
		//				*principalProcForEachBlock, nonEmptyBlocks);
				basicDecomposer.Decompose(
						*principalProcForEachBlock);
		//		basicDecomposer.Validate(
		//				*principalProcForEachBlock);
			}
			timings[hemelb::reporting::Timers::initialDecomposition].Stop();
			log::Logger::Log<log::Info, log::Singleton>("----> basic decomposition (end)");

			// The processor assigned to each block we know about.
			std::unordered_map<site_t, proc_t>* principalProcForEachBlockFiltered =
				new std::unordered_map<site_t, proc_t>();

			log::Logger::Log<log::Info, log::Singleton>("----> read blocks (start)");
			// Perform the initial read-in.
			if (participateInTopology)
			{
				// Reopen in the file just between the nodes in the topology decomposition. Read in blocks
				// local to this node.
				file = net::MpiFile::Open(computeComms, dataFilePath, MPI_MODE_RDONLY, fileInfo);

				ReadInBlocksWithHalo(geometry,
						*principalProcForEachBlock,
						*principalProcForEachBlockFiltered,
						computeComms.Rank());
			}
			timings[hemelb::reporting::Timers::fileRead].Stop();
			log::Logger::Log<log::Info, log::Singleton>("----> read blocks (end)");

			timings[hemelb::reporting::Timers::domainDecomposition].Start();
			// Having done an initial decomposition of the geometry, and read in the data, we optimise the
			// domain decomposition.
			if (participateInTopology)
			{
				OptimiseDomainDecomposition(geometry, *principalProcForEachBlock, *principalProcForEachBlockFiltered);
				file.Close();
			}
			// Finish up - close the file, set the timings, deallocate memory.
			HEMELB_MPI_CALL(MPI_Info_free, (&fileInfo));
			timings[hemelb::reporting::Timers::domainDecomposition].Stop();

			delete principalProcForEachBlock;
			delete principalProcForEachBlockFiltered;

			return geometry;
		}

		std::vector<char> GeometryReader::ReadOnAllTasks(sitedata_t nBytes)
		{
			std::vector<char> buffer(nBytes);
			const net::MpiCommunicator& comm = file.GetCommunicator();
			if (comm.Rank() == HEADER_READING_RANK)
			{
				file.Read(buffer);
			}
			comm.Broadcast(buffer, HEADER_READING_RANK);
			return buffer;
		}

		/**
		 * Read in the section at the beginning of the config file.
		 */
		Geometry GeometryReader::ReadPreamble()
		{
			const unsigned preambleBytes = io::formats::geometry::PreambleLength;
			std::vector<char> preambleBuffer = ReadOnAllTasks(preambleBytes);

			// Create an Xdr translator based on the read-in data.
			io::writers::xdr::XdrReader preambleReader = io::writers::xdr::XdrMemReader(&preambleBuffer[0],
					preambleBytes);

			unsigned hlbMagicNumber, gmyMagicNumber, version;
			// Read in housekeeping values.
			preambleReader.readUnsignedInt(hlbMagicNumber);
			preambleReader.readUnsignedInt(gmyMagicNumber);
			preambleReader.readUnsignedInt(version);

			// Check the value of the HemeLB magic number.
			if (hlbMagicNumber != io::formats::HemeLbMagicNumber)
			{
				throw Exception() << "This file does not start with the HemeLB magic number."
					<< " Expected: " << unsigned(io::formats::HemeLbMagicNumber)
					<< " Actual: " << hlbMagicNumber;
			}

			// Check the value of the geometry file magic number.
			if (gmyMagicNumber != io::formats::geometry::MagicNumber)
			{
				throw Exception() << "This file does not have the geometry magic number."
					<< " Expected: " << unsigned(io::formats::geometry::MagicNumber)
					<< " Actual: " << gmyMagicNumber;
			}

			if (version != io::formats::geometry::VersionNumber)
			{
				throw Exception() << "Version number incorrect."
					<< " Supported: " << unsigned(io::formats::geometry::VersionNumber)
					<< " Input: " << version;
			}

			// Variables we'll read.
			// We use temporary vars here, as they must be the same size as the type in the file
			// regardless of the internal type used.
			unsigned int blocksX, blocksY, blocksZ, blockSize;
			double voxelSize;
			util::Vector3D<double> origin;

			// Read in the values.
			preambleReader.readUnsignedInt(blocksX);
			preambleReader.readUnsignedInt(blocksY);
			preambleReader.readUnsignedInt(blocksZ);
			preambleReader.readUnsignedInt(blockSize);
			preambleReader.readDouble(voxelSize);
			for (unsigned int i = 0; i < 3; ++i)
			{
				preambleReader.readDouble(origin[i]);
			}

			// Read the padding unsigned int.
			unsigned paddingValue;
			preambleReader.readUnsignedInt(paddingValue);

			return Geometry(util::Vector3D<site_t>(blocksX, blocksY, blocksZ),
					blockSize);
		}

		/**
		 * Read the header section, with minimal information about each block.
		 */
		void GeometryReader::ReadHeader(site_t blockCount)
		{
			sitedata_t maxLength = 100000; //std::numeric_limits<signed>::max()
			sitedata_t maxBytes  = GetHeaderLength(maxLength);

			std::vector<char> partialBuffer(maxBytes);

			const net::MpiCommunicator& comm = file.GetCommunicator();

			sitedata_t last;
			sitedata_t offset = 0;
			sitedata_t chunks = blockCount/maxLength;

			nonEmptyBlocks = 0;
			while (offset < chunks+1)
			{
				// The last block to be read this iteration
				last = (offset+1)*maxLength;

				// The last iteration
				if (offset == chunks)
				{
					partialBuffer.resize(GetHeaderLength(blockCount%maxLength));
					partialBuffer.shrink_to_fit();
					last = blockCount;
				}

				if (comm.Rank() == HEADER_READING_RANK)
				{
					file.ReadAt(io::formats::geometry::PreambleLength+
							offset*maxBytes, partialBuffer);
				}
				comm.Broadcast(partialBuffer, HEADER_READING_RANK);

				// Create a Xdr translation object to translate from binary
				hemelb::io::writers::xdr::XdrReader preambleReader =
					hemelb::io::writers::xdr::XdrMemReader(&partialBuffer[0], (unsigned int) partialBuffer.capacity());

				// Read in the data
				for (sitedata_t block = offset*maxLength; block < last; block++)
				{
					unsigned int sites, weights, bytes, uncompressedBytes;
					preambleReader.readUnsignedInt(sites);
#ifdef HEMELB_USE_GMYPLUS
					preambleReader.readUnsignedInt(weights);
#endif
					preambleReader.readUnsignedInt(bytes);
					preambleReader.readUnsignedInt(uncompressedBytes);

					if (sites > 0)
					{
						if (std::pow(2,16)-1 < sites)
							throw Exception() << "Too large (sites)!";
						if (std::pow(2,16)-1 < weights)
							throw Exception() << "Too large (weights)!";
						if (std::pow(2,16)-1 < bytes)
							throw Exception() << "Too large (bytes)!";
						if (std::pow(2,16)-1 < uncompressedBytes)
							throw Exception() << "Too large (uncompressedBytes)!";

						// Essential block information
						blockInformation[block].first  = bytes;
						blockInformation[block].second = uncompressedBytes;

						if (blockInformation.at(block).first != bytes)
							throw Exception() << "blockInformation.at(block).first != bytes!";
						if (blockInformation.at(block).second != uncompressedBytes)
							throw Exception() << "blockInformation.at(block).second != uncompressedBytes!";

						// The number of fluid sites on this block.
#ifdef HEMELB_USE_PARMETIS
						fluidSitesOnEachBlock[block] = sites;
#else
						if (comm.Rank() == HEADER_READING_RANK)
							fluidSitesOnEachBlock[block] = sites;
#endif

#ifdef HEMELB_USE_GMYPLUS
						// 'Computational weight' of this block.
						blockWeights[block] = weights;
#endif

						// Count the number of blocks containing fluid sites.
						nonEmptyBlocks++;
					}
				} offset++;
			}
		}

		/**
		 * Read in the necessary blocks from the file.
		 */
		void GeometryReader::ReadInBlocksWithHalo(Geometry& geometry,
				std::unordered_map<site_t, proc_t>& unitForEachBlock,
				std::unordered_map<site_t, proc_t>& unitForEachBlockFiltered,
				const proc_t localRank)
		{
			// Create a list of which blocks to read in.
			timings[hemelb::reporting::Timers::readBlocksPrelim].Start();

			// Populate the list of blocks to read (including a halo one block wide around all
			// local blocks).
			log::Logger::Log<log::Debug, log::OnePerCore>("----> determining blocks to read");
			std::unordered_set<site_t> readBlock = DecideWhichBlocksToReadIncludingHalo(geometry,
					unitForEachBlock,
					unitForEachBlockFiltered,
					localRank);

#ifndef HEMELB_USE_PARMETIS
			unitForEachBlock.clear();
#endif

			if (READING_GROUP_SPACING*(READING_GROUP_SIZE-1) > (computeComms.Size()-1))
				throw Exception() << "Bad reading core configuration!";

			// Next we spread round the lists of which blocks each core needs access to.
			log::Logger::Log<log::Debug, log::OnePerCore>("----> informing reading cores of block needs");
			net::Net net = net::Net(computeComms);
			Needs needs(geometry.GetBlockCount(),
					readBlock,
					util::NumericalFunctions::min(
					READING_GROUP_SIZE, computeComms.Size()),
					READING_GROUP_SPACING,
					net);
			timings[hemelb::reporting::Timers::readBlocksPrelim].Stop();

			// Set the initial offset to the first block, which will be updated as we progress
			// through the blocks.
			MPI_Offset offset = io::formats::geometry::PreambleLength
				+ GetHeaderLength(geometry.GetBlockCount());

			// Iterate over each block... and trim the fat.
#ifndef HEMELB_USE_PARMETIS
			if (!needs.Reader())
			{
				for (site_t nextBlockToRead = 0; nextBlockToRead < geometry.GetBlockCount(); ++nextBlockToRead)
				{
					if (readBlock.find(nextBlockToRead) == readBlock.end())
						blockInformation.erase(nextBlockToRead);
				}
			}
			computeComms.Barrier();
#endif
			log::Logger::Log<log::Info, log::OnePerCore>(
					"----> blockInformation.size(): %lu",
					blockInformation.size());

			log::Logger::Log<log::Debug, log::OnePerCore>("----> ReadInBlocks() (start)");
			timings[hemelb::reporting::Timers::readBlocksAll].Start();
			// Iterate over each block.
			for (site_t nextBlockToRead = 0; nextBlockToRead < geometry.GetBlockCount(); ++nextBlockToRead)
			{
				// Read in the block on all cores (nothing will be done if this core doesn't need the block).
				ReadInBlock(offset,
						geometry,
						needs.CoreNeedingBlock(nextBlockToRead), // this is a list of cores needing block
						needs.CoreReadingBlock(GetReadingCoreForBlock(nextBlockToRead)),
						nextBlockToRead,
						(readBlock.find(nextBlockToRead) != readBlock.end()));

#ifndef HEMELB_USE_PARMETIS
				// Don't need this anymore...
				readBlock.erase(nextBlockToRead);

				// or this nasty thing.
				needs.CoreNeedingBlock_erase(nextBlockToRead);
				if (needs.CoreNeedingBlock(nextBlockToRead).size() != 1)
					throw Exception() << "needs.ProcessorsNeedingBlock(nextBlockToRead) not clear!";
#endif

				// Update the offset to be ready for the next block.
				if (blockInformation.find(nextBlockToRead) != blockInformation.end())
				{
					offset += blockInformation.at(nextBlockToRead).first;
#ifndef HEMELB_USE_PARMETIS
					blockInformation.erase(nextBlockToRead);
#endif
				}
			}
			timings[hemelb::reporting::Timers::readBlocksAll].Stop();
			log::Logger::Log<log::Debug, log::OnePerCore>("----> ReadInBlocks() (end)");
		}

		void GeometryReader::ReadInBlock(MPI_Offset offsetSoFar, Geometry& geometry,
				const std::vector<proc_t>& procsWantingThisBlock,
				const proc_t readingCore, const site_t blockNumber, const bool neededOnThisRank)
		{
			// Easy case if there are no sites on the block.
			if (blockInformation.find(blockNumber) == blockInformation.end())
			{
				return;
			}

			net::Net net = net::Net(computeComms);

			std::vector<char> compressedBlockData;
			if (readingCore == computeComms.Rank())
			{
				log::Logger::Log<log::Debug, log::OnePerCore>("------> blockNumber = %li with:\n"
						"bytes             = %u,\n"
						"uncompressedBytes = %u and offsetSoFar = %li",
						blockNumber,
						blockInformation.at(blockNumber).first,
						blockInformation.at(blockNumber).second,
						offsetSoFar);

				timings[hemelb::reporting::Timers::readBlock].Start();
				// Read the data.
				compressedBlockData.resize(blockInformation.at(blockNumber).first);
				file.ReadAt(offsetSoFar, compressedBlockData);

				// Spread it...
				// unless procsWantingBlocksBuffer (procsWantingThisBlock) is empty.
				if (procsWantingThisBlock.front() != -1) {
					for (std::vector<proc_t>::const_iterator receiver = procsWantingThisBlock.begin(); receiver
							!= procsWantingThisBlock.end(); receiver++)
					{
						if (*receiver != computeComms.Rank())
						{
#ifdef HEMELB_USE_MPI_CALL
							MPI_Send(&compressedBlockData[0],
									compressedBlockData.size(),
									MPI_CHAR, *receiver, 0, computeComms);
#else
							net.RequestSendV(compressedBlockData, *receiver);
#endif
						}
					}
					timings[hemelb::reporting::Timers::readBlock].Stop();
				}
			}
			else if (neededOnThisRank)
			{
				compressedBlockData.resize(blockInformation.at(blockNumber).first);
#ifdef HEMELB_USE_MPI_CALL
				MPI_Recv(&compressedBlockData[0],
						compressedBlockData.size(),
						MPI_CHAR, readingCore, 0, computeComms, MPI_STATUS_IGNORE);
#else
				net.RequestReceiveV(compressedBlockData, readingCore);
#endif
			}
			else
			{
				return;
			}
			timings[hemelb::reporting::Timers::readNet].Start();
#ifndef HEMELB_USE_MPI_CALL
			net.Dispatch();
#endif
			timings[hemelb::reporting::Timers::readNet].Stop();

			timings[hemelb::reporting::Timers::readParse].Start();
			if (neededOnThisRank)
			{
				// Create an Xdr interpreter.
				std::vector<char> blockData = DecompressBlockData(compressedBlockData,
						blockInformation.at(blockNumber).second);
				io::writers::xdr::XdrMemReader lReader(&blockData.front(), blockData.size());

				ParseBlock(geometry, blockNumber, lReader);
			}
			//else if (geometry.Blocks.find(blockNumber) != geometry.Blocks.end())
			//{
			//	geometry.Blocks[blockNumber].Sites = std::vector<GeometrySite>(0, GeometrySite(false));
			//}
			timings[hemelb::reporting::Timers::readParse].Stop();
		}

		std::vector<char> GeometryReader::DecompressBlockData(const std::vector<char>& compressed,
				const unsigned int uncompressedBytes)
		{
			timings[hemelb::reporting::Timers::unzip].Start();
			// For zlib return codes.
			int ret;

			// Set up the buffer for decompressed data.
			std::vector<char> uncompressed(uncompressedBytes);

			// Set up the inflator.
			z_stream stream;
			stream.zalloc = Z_NULL;
			stream.zfree = Z_NULL;
			stream.opaque = Z_NULL;
			stream.avail_in = compressed.size();
			stream.next_in = reinterpret_cast<unsigned char*> (const_cast<char*> (&compressed.front()));

			ret = inflateInit(&stream);
			if (ret != Z_OK)
				throw Exception() << "Decompression error for block!";

			stream.avail_out = uncompressed.size();
			stream.next_out = reinterpret_cast<unsigned char*> (&uncompressed.front());

			ret = inflate(&stream, Z_FINISH);
			if (ret != Z_STREAM_END)
				throw Exception() << "Decompression error for block!";

			uncompressed.resize(uncompressed.size() - stream.avail_out);
			ret = inflateEnd(&stream);
			if (ret != Z_OK)
				throw Exception() << "Decompression error for block!";

			timings[hemelb::reporting::Timers::unzip].Stop();
			return uncompressed;
		}

		void GeometryReader::ParseBlock(Geometry& geometry, const site_t block,
				io::writers::xdr::XdrReader& reader)
		{
			// We start by clearing the sites on the block. We read the blocks twice (once before
			// optimisation and once after), so there can be sites on the block from the previous read.
			geometry.Blocks[block].Sites.clear();

			for (site_t localSiteIndex = 0; localSiteIndex < geometry.GetSitesPerBlock(); ++localSiteIndex)
			{
				geometry.Blocks.at(block).Sites.push_back(ParseSite(reader));
			}
		}

		GeometrySite GeometryReader::ParseSite(io::writers::xdr::XdrReader& reader)
		{
			// Read the fluid property.
			unsigned isFluid;
			bool success = reader.readUnsignedInt(isFluid);

			if (!success)
			{
				log::Logger::Log<log::Error, log::OnePerCore>("Error reading site type!");
			}

			/// @todo #598 use constant in hemelb::io::formats::geometry
			GeometrySite readInSite(isFluid != 0);

			// If solid, there's nothing more to do.
			if (!readInSite.isFluid)
			{
				return readInSite;
			}

			const io::formats::geometry::DisplacementVector& neighbourhood =
				io::formats::geometry::Get().GetNeighbourhood();
			// Prepare the links array to have enough space.
			readInSite.links.resize(latticeInfo.GetNumVectors() - 1);

			bool isGmyWallSite = false;

			// For each link direction...
			for (Direction readDirection = 0; readDirection < neighbourhood.size(); readDirection++)
			{
				// read the type of the intersection and create a link...
				unsigned intersectionType;
				reader.readUnsignedInt(intersectionType);

				GeometrySiteLink link;
				link.type = (GeometrySiteLink::IntersectionType) intersectionType;

				// walls have a floating-point distance to the wall...
				if (link.type == GeometrySiteLink::WALL_INTERSECTION)
				{
					isGmyWallSite = true;
					float distance;
					reader.readFloat(distance);
					link.distanceToIntersection = distance;
				}

				// inlets and outlets (which together with none make up the other intersection types)
				// have an iolet id and a distance float...
				else if (link.type != GeometrySiteLink::NO_INTERSECTION)
				{
					float distance;
					unsigned ioletId;
					reader.readUnsignedInt(ioletId);
					reader.readFloat(distance);

					link.ioletId = ioletId;
					link.distanceToIntersection = distance;
				}

				// now, attempt to match the direction read from the local neighbourhood to one in the
				// lattice being used for simulation. If a match is found, assign the link to the read
				// site.
				for (Direction usedLatticeDirection = 1; usedLatticeDirection < latticeInfo.GetNumVectors(); usedLatticeDirection++)
				{
					if (latticeInfo.GetVector(usedLatticeDirection) == neighbourhood[readDirection])
					{
						// If this link direction is necessary to the lattice in use, keep the link data.
						readInSite.links[usedLatticeDirection - 1] = link;
						break;
					}
				}
			}

			unsigned normalAvailable;
			reader.readUnsignedInt(normalAvailable);
			readInSite.wallNormalAvailable = (normalAvailable
					== io::formats::geometry::WALL_NORMAL_AVAILABLE);

			if (readInSite.wallNormalAvailable != isGmyWallSite)
			{
				std::string msg = isGmyWallSite
					? "wall fluid site without"
					: "bulk fluid site with";
				throw Exception() << "Malformed GMY file, "
					<< msg << " a defined wall normal currently not allowed.";
			}

			if (readInSite.wallNormalAvailable)
			{
				reader.readFloat(readInSite.wallNormal[0]);
				reader.readFloat(readInSite.wallNormal[1]);
				reader.readFloat(readInSite.wallNormal[2]);
			}

			return readInSite;
		}

		proc_t GeometryReader::GetReadingCoreForBlock(site_t blockNumber)
		{
			return proc_t(blockNumber % util::NumericalFunctions::min(READING_GROUP_SIZE,
						computeComms.Size()));
		}

		std::unordered_set<site_t> GeometryReader::DecideWhichBlocksToReadIncludingHalo(
				const Geometry& geometry,
				const std::unordered_map<site_t, proc_t>& unitForEachBlock,
				std::unordered_map<site_t, proc_t>& unitForEachBlockFiltered,
				proc_t localRank)
		{
			std::unordered_set<site_t> shouldReadBlock;

			// Read a block in if it has fluid sites and is to live on the current processor. Also read
			// in any neighbours with fluid sites.
			for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; ++blockI)
			{
				for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; ++blockJ)
				{
					for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; ++blockK)
					{
						site_t lBlockId = geometry.GetBlockIdFromBlockCoordinates(blockI, blockJ, blockK);

						if (unitForEachBlock.find(lBlockId) == unitForEachBlock.end())
						{
							continue;
						}

						if (unitForEachBlock.at(lBlockId) != localRank)
						{
							continue;
						}

						// Read in all neighbouring blocks.
						for (site_t neighI = util::NumericalFunctions::max<site_t>(0, blockI - 1); (neighI
									<= (blockI + 1)) && (neighI < geometry.GetBlockDimensions().x); ++neighI)
						{
							for (site_t neighJ = util::NumericalFunctions::max<site_t>(0, blockJ - 1); (neighJ
										<= (blockJ + 1)) && (neighJ < geometry.GetBlockDimensions().y); ++neighJ)
							{
								for (site_t neighK = util::NumericalFunctions::max<site_t>(0, blockK - 1); (neighK
											<= (blockK + 1)) && (neighK < geometry.GetBlockDimensions().z); ++neighK)
								{
									site_t lNeighId = geometry.GetBlockIdFromBlockCoordinates(neighI, neighJ, neighK);
									if (unitForEachBlock.find(lNeighId) != unitForEachBlock.end())
										unitForEachBlockFiltered[lNeighId] = unitForEachBlock.at(lNeighId);
									shouldReadBlock.insert(lNeighId);
								}
							}
						}
					}
				}
			}
			return shouldReadBlock;
		}

		void GeometryReader::OptimiseDomainDecomposition(Geometry& geometry,
				std::unordered_map<site_t, proc_t>& procForEachBlock,
				std::unordered_map<site_t, proc_t>& procForEachBlockFiltered)
		{
			decomposition::OptimisedDecomposition optimiser(timings,
					computeComms,
					geometry,
					latticeInfo,
					procForEachBlock,
					fluidSitesOnEachBlock);

			site_t geometrySize = geometry.Blocks.size();

//#ifdef HEMELB_USE_ZOLTAN
#ifdef HEMELB_USE_PARMETIS
			timings[hemelb::reporting::Timers::reRead].Start();
			log::Logger::Log<log::Debug, log::OnePerCore>("----> rereading blocks");
			// Reread the blocks based on the optimised decomposition.
			RereadBlocks(geometry,
					optimiser.GetMovesCountPerCore(),
					optimiser.GetMovesList(),
					procForEachBlock,
					procForEachBlockFiltered);
			timings[hemelb::reporting::Timers::reRead].Stop();

			log::Logger::Log<log::Info, log::OnePerCore>(
					"----> geometry.Block.size()/geometrySize: %f",
					geometry.Blocks.size()/(double)geometrySize);
#endif

			timings[hemelb::reporting::Timers::moves].Start();
			// Implement the decomposition now that we have read the necessary data.
			log::Logger::Log<log::Debug, log::OnePerCore>("----> implementing moves");
			ImplementMoves(geometry,
					procForEachBlock,
					procForEachBlockFiltered,
					optimiser.GetMovesCountPerCore(),
					optimiser.GetMovesList());
			timings[hemelb::reporting::Timers::moves].Stop();

			ShowDecomposition(geometry, procForEachBlockFiltered);
		}

		void GeometryReader::ShowDecomposition(Geometry& geometry,
				const std::unordered_map<site_t, proc_t>& procForEachBlockFiltered) const
		{
			// Open file for writing.
			std::ofstream myfile;
			myfile.open("./decomposition/globalSiteCoords_" + std::to_string(computeComms.Rank()) + ".txt");

			const int blockSize = geometry.GetBlockSize();
			for (site_t blockI = 0; blockI < geometry.GetBlockDimensions().x; blockI++)
				for (site_t blockJ = 0; blockJ < geometry.GetBlockDimensions().y; blockJ++)
					for (site_t blockK = 0; blockK < geometry.GetBlockDimensions().z; blockK++)
					{
						const site_t blockNumber = geometry.GetBlockIdFromBlockCoordinates(
								blockI,
								blockJ,
								blockK);

						// Does block contain fluid sites?
						if (procForEachBlockFiltered.find(blockNumber) == procForEachBlockFiltered.end())
							continue;
						// Is block owned by this rank?
						if (procForEachBlockFiltered.at(blockNumber) != computeComms.Rank())
							continue;

						const BlockReadResult& blockReadResult = geometry.Blocks.at(blockNumber);

						//site_t blockIJData = blockNumber / geometry.GetBlockDimensions().z;
						//blockCoords.x = blockIJData / geometry.GetBlockDimensions().y;
						//blockCoords.y = blockIJData % geometry.GetBlockDimensions().y;
						//blockCoords.z = blockNumber % geometry.GetBlockDimensions().z;

						site_t m = 0;
						// Iterate over sites within the block.
						for (site_t localSiteI = 0; localSiteI < blockSize; localSiteI++)
							for (site_t localSiteJ = 0; localSiteJ < blockSize; localSiteJ++)
								for (site_t localSiteK = 0; localSiteK < blockSize; localSiteK++)
								{
									//site_t localSiteID =
									//	geometry.GetSiteIdFromSiteCoordinates(localSiteI, localSiteJ, localSiteK);
									//site_t siteIJData = localSiteID / block_size;
									//localSiteCoords.x = siteIJData  / block_size;
									//localSiteCoords.y = siteIJData  % block_size;
									//localSiteCoords.z = localSiteID % block_size;
									//util::Vector3D<site_t> globalSiteCoords = blockCoords * blockSize + localSiteCoords;

								  //if (blockReadResult.Sites[m].isFluid)
									if (blockReadResult.Sites[m].targetProcessor != SITE_OR_BLOCK_SOLID)
										myfile <<
											blockReadResult.Sites[m].isFluid << " " <<
											blockReadResult.Sites[m].targetProcessor << " " <<
											blockNumber << " " <<
											blockI * blockSize + localSiteI << " " <<
											blockJ * blockSize + localSiteJ << " " <<
											blockK * blockSize + localSiteK << std::endl; m++;
								}
					} myfile.close();
		}

		// The header section of the config file contains a number of records.
		site_t GeometryReader::GetHeaderLength(site_t blockCount) const
		{
			return io::formats::geometry::HeaderRecordLength * blockCount;
		}

		void GeometryReader::RereadBlocks(Geometry& geometry, const std::vector<idx_t>& movesPerProc,
				const std::vector<idx_t>& movesList,
				std::unordered_map<site_t, proc_t>& procForEachBlock,
				std::unordered_map<site_t, proc_t>& procForEachBlockFiltered)
		{
			// Initialise the map that will store changes (to be reverted).
			std::unordered_map<site_t, proc_t> origProcForEachBlock;
			std::unordered_map<site_t, proc_t>::const_iterator got;

			// Set the proc for each block to be the current proc whenever a site on that block is
			// going to be moved to the current proc.
			idx_t moveIndex = 0;
			for (proc_t fromProc = 0; fromProc < computeComms.Size(); ++fromProc)
			{
				for (idx_t moveNumber = 0; moveNumber < movesPerProc[fromProc]; ++moveNumber)
				{
					idx_t block  = movesList[3 * moveIndex];
					idx_t toProc = movesList[3 * moveIndex + 2];

					++moveIndex;

					if (toProc == (idx_t) computeComms.Rank())
					{
						got = procForEachBlock.find(block);
						origProcForEachBlock.insert({got->first,got->second});
							procForEachBlock.at(block) = computeComms.Rank();
					}
				}
			}

			// Reread the blocks into the GlobalLatticeData now.
			ReadInBlocksWithHalo(geometry, procForEachBlock, procForEachBlockFiltered, computeComms.Rank());

			// Revert changes to procForEachBlock.
			for (std::unordered_map<site_t, proc_t>::const_iterator iter = origProcForEachBlock.begin(); iter != origProcForEachBlock.end(); ++iter)
			{
				procForEachBlock.at(iter->first) = iter->second;
			}
		}

		void GeometryReader::ImplementMoves(Geometry& geometry,
				const std::unordered_map<site_t, proc_t>& procForEachBlock,
				const std::unordered_map<site_t, proc_t>& procForEachBlockFiltered,
				const std::vector<idx_t>& movesFromEachProc,
				const std::vector<idx_t>& movesList) const
		{
			log::Logger::Log<log::Debug, log::OnePerCore>("----> ImplementMoves(): procForEachBlockFiltered.size() == %i", procForEachBlockFiltered.size());
			// First all, set the proc rank for each site to what it originally was before
			// domain decomposition optimisation. Go through each block.
			for (site_t block = 0; block < geometry.GetBlockCount(); ++block)
			{
				// If this proc has owned a fluid site on this block either before or after optimisation,
				// the following will be non-null.
				if (geometry.Blocks.find(block) != geometry.Blocks.end())
				{
#ifdef HEMELB_USE_PARMETIS
					// Get the original proc for that block.
					proc_t originalProc = procForEachBlock.at(block);
#else
					// Get the original proc for that block.
					proc_t originalProc = procForEachBlockFiltered.at(block);
#endif
					// For each site on that block...
					for (site_t siteIndex = 0; siteIndex < geometry.GetSitesPerBlock(); ++siteIndex)
					{
						// if the site is non-solid...
						if (geometry.Blocks.at(block).Sites[siteIndex].targetProcessor != SITE_OR_BLOCK_SOLID)
						{
							// set its rank to be the rank it had before optimisation.
							geometry.Blocks.at(block).Sites[siteIndex].targetProcessor
								= ConvertTopologyRankToGlobalRank(originalProc);
						}
					}
				}
			}

//#ifdef HEMELB_USE_ZOLTAN
#ifdef HEMELB_USE_PARMETIS
			// Now implement the moves suggested by ParMETIS.
			idx_t moveIndex = 0;

			// For each source proc, go through as many moves as it had.
			for (proc_t fromProc = 0; fromProc < computeComms.Size(); ++fromProc)
			{
				for (idx_t moveNumber = 0; moveNumber < movesFromEachProc[fromProc]; ++moveNumber)
				{
					// For each move, get the block, site and destination proc.
					idx_t block  = movesList[3 * moveIndex];
					idx_t site   = movesList[3 * moveIndex + 1];
					idx_t toProc = movesList[3 * moveIndex + 2];

					// Only implement the move if we have read that block's data.
					if (geometry.Blocks.find(block) != geometry.Blocks.end())
					{
						// Implement the move.
						geometry.Blocks.at(block).Sites[site].targetProcessor
							= ConvertTopologyRankToGlobalRank((proc_t) toProc);
					} ++moveIndex;
				}
			}
#endif
		}

		proc_t GeometryReader::ConvertTopologyRankToGlobalRank(proc_t topologyRankIn) const
		{
			// If the global rank is not equal to the topology rank, we are not using rank 0 for
			// LBM.
			return (hemeLbComms.Rank() == computeComms.Rank())
				? topologyRankIn
				: (topologyRankIn + 1);
		}

		bool GeometryReader::ShouldValidate() const
		{
#ifdef HEMELB_VALIDATE_GEOMETRY
			return true;
#else
			return false;
#endif
		}
	}
}
