
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <chrono>
#include <random>

#include "colloids/ParticleSet.h"
#include "colloids/BodyForces.h"
#include "colloids/BoundaryConditions.h"
#include "log/Logger.h"
#include "io/writers/xdr/XdrMemWriter.h"
#include "io/formats/formats.h"
#include "io/formats/colloids.h"

namespace hemelb
{
	namespace colloids
	{
		struct ParticleSorter
		{
			int localRank;
			ParticleSorter(int rank) :
				localRank(rank)
			{
			}

			bool operator()(const Particle& a, const Particle& b)
			{
				// ORDER BY isLocal, ownerRank, particleId
				if (a.ownerRank == b.ownerRank)
					return (a.particleId < b.particleId);
				else if (a.ownerRank == localRank)
					return true;
				else if (b.ownerRank == localRank)
					return false;
				else
					return (a.ownerRank < b.ownerRank);
			}
		};

		ParticleSet::ParticleSet(geometry::LatticeData& latDatLBM,
				io::xml::Element& particlesElem,
				lb::MacroscopicPropertyCache& propertyCache,
				const configuration::SimConfig* simConfig,
				std::vector<proc_t>& neighbourProcessors,
				const net::IOCommunicator& ioComms_,
				const std::string& outputPath) :
			ioComms(ioComms_), localRank(ioComms.Rank()), latDatLBM(latDatLBM), propertyCache(propertyCache), path(outputPath), net(ioComms)
		{
			// open the file, unless it already exists, for writing only, creating it if it doesn't exist
			file = net::MpiFile::Open(ioComms, path, MPI_MODE_EXCL | MPI_MODE_WRONLY | MPI_MODE_CREATE);

			if (ioComms.OnIORank())
			{
				// write out the header information from core 0
				buffer.resize(io::formats::colloids::MagicLength);
				io::writers::xdr::XdrMemWriter writer(&buffer[0], io::formats::colloids::MagicLength);
				writer << (uint32_t) io::formats::HemeLbMagicNumber;
				writer << (uint32_t) io::formats::colloids::MagicNumber;
				writer << (uint32_t) io::formats::colloids::VersionNumber;
				file.Write(buffer);
			}

			HEMELB_MPI_CALL(MPI_File_seek_shared, (file, 0, MPI_SEEK_END));

			// add an element into scanMap for each neighbour rank with zero for both counts
			// sorting the list of neighbours allows the position in the map to be predicted
			// & giving the correct position in the map makes insertion significantly faster
			std::sort(neighbourProcessors.begin(), neighbourProcessors.end());
			for (std::vector<proc_t>::const_iterator iter = neighbourProcessors.begin();
					iter != neighbourProcessors.end(); iter++)
				scanMap.insert(scanMap.end(), scanMapContentType(*iter,
							scanMapElementType(0, 0)));

			//if (localRank == 1) {
			//	for (unsigned long rank = 2; rank < ioComms.Size(); rank++)
			//		scanMap.insert(scanMap.end(), scanMapContentType(rank,
			//					scanMapElementType(0, 0)));
			//} else {
			//	for (std::vector<proc_t>::const_iterator iter = neighbourProcessors.begin();
			//			iter != neighbourProcessors.end(); iter++)
			//		scanMap.insert(scanMap.end(), scanMapContentType(*iter,
			//					scanMapElementType(0, 0)));
			//}

			//if (localRank > 1 && (scanMap.find(1) == scanMap.end()))
			//	scanMap.insert(scanMap.begin(), scanMapContentType(1,
			//				scanMapElementType(0, 0)));

			// the local rank is added last, because its position cannot be easily predicted
			scanMap.insert(scanMapContentType(localRank, scanMapElementType(0, 0)));

			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank            = iterMap->first;
				const unsigned int& numberOfParticles  = iterMap->second.first;
				const unsigned int& numberOfVelocities = iterMap->second.second;
				log::Logger::Log<log::Debug, log::OnePerCore>("scanMap[%i] = {%i, %i}",
						neighbourRank,
						numberOfParticles,
						numberOfVelocities);
			}

			// assume we are at the <particles> node
			bool first = true; particlesFromFile = 0;
			for (io::xml::Element particleElem = particlesElem.GetChildOrNull("subgridParticle");
					particleElem != io::xml::Element::Missing();
					particleElem  = particleElem.NextSiblingOrNull("subgridParticle"))
			{
				// create the particle object from the settings in the config file
				//Particle nextParticle(latDatLBM, lbmParams, simConfig, particleElem);
				Particle nextParticle(latDatLBM, simConfig, particleElem);

				if (first)
				{
					propertyCache.velocityCache.SetRefreshFlag();
					// the first particle (particleId = 0) is the master particle
					// it is used to create replica particles
					memcpy(&masterParticle, &nextParticle, sizeof(nextParticle));
					first = false; continue;
				}

				// count the total number of particles from file
				// particleId of particles due to insertion is offset by this amount
				particlesFromFile++;

				// check the particle is valid, i.e. in fluid, and is locally owned
				if (nextParticle.IsValid() && nextParticle.GetOwnerRank() == localRank)
				{
					// add the particle to the list of known particles ...
					particles.push_back(nextParticle);
					// ... and keep the count of local particles up-to-date
					scanMap[localRank].first++;
				}
			}

			io::xml::Element radiusElem = particlesElem.GetChildOrThrow("sphereRadius");
			configuration::GetDimensionalValue(radiusElem, "lattice", sphereRadius);

			io::xml::Element insertElem = particlesElem.GetChildOrThrow("emissionCount");
			configuration::GetDimensionalValue(insertElem, "dimensionless", nEmitted);

			io::xml::Element centreElem = particlesElem.GetChildOrThrow("sphereCentre");
			centreElem.GetAttributeOrThrow("x", sphereCentre.x);
			centreElem.GetAttributeOrThrow("y", sphereCentre.y);
			centreElem.GetAttributeOrThrow("z", sphereCentre.z);
		}

		ParticleSet::~ParticleSet()
		{
			particles.clear();
		}

		const void ParticleSet::OutputInformation(const LatticeTimeStep currentTimestep)
		{
			// ensure the buffer is large enough
			const unsigned int maxSize = io::formats::colloids::RecordLength * particles.size()
				+ io::formats::colloids::HeaderLength;
			if (buffer.size() < maxSize)
			{
				buffer.resize(maxSize);
			}

			// create an XDR writer and write all the particles for this processor
			io::writers::xdr::XdrMemWriter writer(&buffer.front(), maxSize);

			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				if ((particle.GetOwnerRank() == localRank) && (particle.GetCreationTimestep() < currentTimestep+2))
				{
					particle.WriteToStream(currentTimestep, * ((io::writers::Writer*) &writer));
				}
			}

			// get the number of bytes written
			const unsigned int count = writer.getCurrentStreamPosition();

			// find how far we currently are into the file
			MPI_Offset positionBeforeWriting;
			HEMELB_MPI_CALL(MPI_File_get_position_shared, (file, &positionBeforeWriting));

			log::Logger::Log<log::Debug, log::OnePerCore>("from offsetEOF: %i\n", positionBeforeWriting);

			// go past the header (which we'll write at the end)
			unsigned int sizeOfHeader = io::formats::colloids::HeaderLength;
			HEMELB_MPI_CALL(MPI_File_seek_shared, (file, sizeOfHeader, MPI_SEEK_END));

			// collective write: the effect is as though all writes are done
			// in serialised order, i.e. as if rank 0 writes first, followed
			// by rank 1, and so on, until all ranks have written their data
			HEMELB_MPI_CALL(MPI_File_write_ordered,
					(file, &buffer.front(), count, MPI_CHAR, MPI_STATUS_IGNORE));

			// the collective ordered write modifies the shared file pointer
			// it should point to the byte following the highest rank's data
			// (should be true for all ranks but) we only need it for rank 0
			MPI_Offset positionAferWriting;
			HEMELB_MPI_CALL(MPI_File_get_position_shared, (file, &positionAferWriting));

			log::Logger::Log<log::Debug, log::OnePerCore>("new offsetEOF: %i\n", positionBeforeWriting);

			// now write the header section, only on rank 0
			if (ioComms.OnIORank())
			{
				writer << (uint32_t) io::formats::colloids::HeaderLength;
				writer << (uint32_t) io::formats::colloids::RecordLength;
				writer << (uint64_t) (positionAferWriting - positionBeforeWriting - io::formats::colloids::HeaderLength);
				writer << (uint64_t) currentTimestep;
				HEMELB_MPI_CALL(MPI_File_write_at,
						(file, positionBeforeWriting, &buffer[count], sizeOfHeader, MPI_CHAR, MPI_STATUS_IGNORE));
			}

			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank            = iterMap->first;
				const unsigned int& numberOfParticles  = iterMap->second.first;
				const unsigned int& numberOfVelocities = iterMap->second.second;
				log::Logger::Log<log::Debug, log::OnePerCore>("scanMap[%i] = {%i, %i}",
						neighbourRank,
						numberOfParticles,
						numberOfVelocities);
			}
		}

		const void ParticleSet::CalculateBodyForces()
		{
			for (std::vector<Particle>::iterator iter = particles.begin();
					iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				if (particle.GetOwnerRank() == localRank)
					particle.CalculateBodyForces();
			}
		}

		const void ParticleSet::InterpolateFluidVelocity(const LatticeTimeStep currentTimestep)
		{
			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				particle.InterpolateFluidVelocity(latDatLBM, propertyCache);
			}
			propertyCache.velocityCache.SetRefreshFlag();
		}

		const void ParticleSet::UpdateVelocities()
		{
			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				if (particle.GetOwnerRank() == localRank)
					particle.UpdateVelocity();
			}
		}

		const void ParticleSet::ResetVelocities()
		{
			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				if (particle.GetOwnerRank() == localRank)
					particle.ResetVelocity();
			}
		}

		const void ParticleSet::EmitParticles(const LatticeTimeStep currentTimestep)
		{
			// clear list of particles to broadcast
			Particles.clear();
			Particles.resize(nEmitted);

			// generate particles on IO rank only
			if (ioComms.OnIORank())
			{
				// construct a trivial random generator engine from a time-based seed
				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
				std::default_random_engine generator(seed);
				std::uniform_real_distribution<> distribution(-0.5,0.5);

				// the new particle (based on masterParticle)
				Particle newParticle;

				LatticePosition globalPosition;
				for (unsigned long iEmitted = 0; iEmitted < nEmitted; iEmitted++)
				{
					// replicate masterParticle to create a new particle
					memcpy(&newParticle, &masterParticle, sizeof(masterParticle));

					// generate globalPosition
					globalPosition.x = distribution(generator);
					globalPosition.y = distribution(generator);
					globalPosition.z = distribution(generator);
					LatticePosition normalise = globalPosition.GetNormalised();
					globalPosition = sphereCentre + normalise*sphereRadius;

					// set globalPosition
					newParticle.SetGlobalPosition(globalPosition);

					// set the particle's particleId
					newParticle.SetParticleId(particlesFromFile+tEmitted);

					// set the particle's timestep of creation
					newParticle.SetCreationTimestep(currentTimestep);

					// add the particle to list of known particles
					//particlesList.push_front(newParticle);
					Particles[iEmitted] = newParticle;

					// update count of emitted particles
					// for the purpose of generating a particleId
					tEmitted++;
				}
			}

			// broadcast particlesList
			std::vector<Particle>::iterator iterBegin = Particles.begin();
			for (proc_t neighbourRank = 0; neighbourRank < ioComms.Size(); neighbourRank++)
			{
				if (neighbourRank != ioComms.GetIORank()) {

					const unsigned int& numberOfParticlesToRecv_ = nEmitted;
					const unsigned int& numberOfParticlesToSend_ = nEmitted;

					if (ioComms.OnIORank())
						net.RequestSend(
								& ((PersistedParticle&) *iterBegin),
								numberOfParticlesToSend_, neighbourRank);
					else
						net.RequestReceive(
								& ((PersistedParticle&) *iterBegin),
								numberOfParticlesToRecv_, ioComms.GetIORank());
				}
			}
			net.Dispatch();

			// create list from vector
			std::list<Particle> particlesList(particles.begin(), particles.end());

			// iterate through new particles ...
			// and add particle to particlesList if I own particle
			for (std::vector<Particle>::iterator iter = Particles.begin(); iter != Particles.end(); iter++)
			{
				Particle& particle = *iter;
				particle.UpdatePosition(latDatLBM);

				//log::Logger::Log<log::Info, log::OnePerCore>("%f", particle.GetGlobalPosition().x);

				// check the particle is valid, i.e. in fluid, and is locally owned
				if (particle.IsValid() && particle.GetOwnerRank() == localRank)
				{
					// add the particle to the list of known particles ...
					particlesList.push_front(particle);
					// and keep the count of local particles up-to-date
					scanMap[localRank].first++;
				}
			}

			// count number of particles emitted on this rank
			int nEmittedHere = particlesList.size()-particles.size();
			if (nEmittedHere != 0)
				log::Logger::Log<log::Info, log::OnePerCore>("particles emitted (target): %03lu, particles emitted: %03i",
						nEmitted, nEmittedHere);

			// resize vector to recieve new particles
			particles.resize(particlesList.size());

			// list (particlesList) to vector (particles)
			std::copy(particlesList.begin(), particlesList.end(), particles.begin());

			//globalPosition.x = sphereCentre.x + (-0.5+iEmitted/(double)nEmitted)*sphereRadius;
			//globalPosition.y = sphereCentre.y;
			//globalPosition.z = sphereCentre.z;
		}

		const void ParticleSet::ApplyBoundaryConditions(const LatticeTimeStep currentTimestep)
		{
			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				if (particle.GetOwnerRank() == localRank)
				{
					BoundaryConditions::DoSomeThingsToParticle(currentTimestep, particle);
					if (particle.IsReadyToBeDeleted())
						log::Logger::Log<log::Debug, log::OnePerCore>("In ParticleSet::ApplyBoundaryConditions - timestep: %lu, particleId: %lu, IsReadyToBeDeleted: %s, markedForDeletion: %lu, lastCheckpoint: %lu\n",
								currentTimestep,
								particle.GetParticleId(),
								particle.IsReadyToBeDeleted() ?
								"YES" :
								"NO",
								particle.GetDeletionMarker(),
								particle.GetLastCheckpointTimestep());
				}
			}

			// shuffle (or partition) the particles in our vector containing all particles
			// so the first partition contains all the local particles that should be kept
			// and the other contains all the deletable-local plus the non-local particles
			std::vector<Particle>::iterator bound =
				std::partition(particles.begin(),
						particles.begin() + scanMap[localRank].first,
						std::not1(std::mem_fun_ref(&Particle::IsReadyToBeDeleted)));

			if (scanMap[localRank].first > (bound - particles.begin()))
				log::Logger::Log<log::Debug, log::OnePerCore>("In ParticleSet::ApplyBoundaryConditions - timestep: %lu, scanMap[localRank].first: %lu, bound-particles.begin(): %lu\n",
						currentTimestep,
						scanMap[localRank].first,
						bound - particles.begin());

			// the partitioning above may invalidate the scanMap used by the communication
			// the next communication function called is CommunicatePositions, which needs
			// the number of local particles to be correct - i.e. scanMap[localRank].first
			// - the rest of scanMap will be re-built by the CommunicatePositions function
			scanMap[localRank].first = bound - particles.begin();
		}

		const void ParticleSet::UpdatePositions(const LatticeTimeStep currentTimestep)
		{
			if (log::Logger::ShouldDisplay<log::Debug>())
				log::Logger::Log<log::Debug, log::OnePerCore>("In colloids::ParticleSet::UpdatePositions #particles == %i ...\n",
						localRank,
						particles.size());

			// only update the position for particles that are locally owned because
			// only the owner has velocity contributions from all neighbouring ranks
			ghostOnProc.clear(); unsigned int particleIdx = 0;
			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle = *iter;
				if (particle.GetOwnerRank() == localRank)
					particle.UpdatePosition(latDatLBM, currentTimestep, ghostOnProc, particleIdx);
				particleIdx++;
			}
		}

		const void ParticleSet::CommunicateParticlePositions(const LatticeTimeStep currentTimestep)
		{
			/** CommunicateParticlePositions
			 *
			 *  The global position of each particle is updated by the ownerRank process.
			 *  The ownerRank for each particle is verified when its position is updated.
			 *  Some (previously locally owned) particles may no longer be locally owned.
			 *
			 */

			// delete duplicates from migration list
			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank = iterMap->first;
				if (ghostOnProc.find(neighbourRank) != ghostOnProc.end())
					ghostOnProc.at(neighbourRank).erase(
							unique(
								ghostOnProc.at(neighbourRank).begin(),
								ghostOnProc.at(neighbourRank).end()),
							ghostOnProc.at(neighbourRank).end());
			}

			const unsigned int  particlesOffset  = scanMap[localRank].first;
			const unsigned int& particlesOffset_ = scanMap[localRank].first;

			if (scanMap.size() < 2)
			{
				return;
			}

			// reset the scanMap
			for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
				iterMap->second = scanMapElementType(0, 0);

			// the number of particles owned by this rank
			scanMap[localRank].first = particlesOffset;

			for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank = iterMap->first;
				if (neighbourRank != localRank) {

					if (ghostOnProc.find(neighbourRank) != ghostOnProc.end())
						iterMap->second.second = ghostOnProc.at(neighbourRank).size();

					unsigned int& numberOfParticlesToRecv_ = iterMap->second.first;
					unsigned int& numberOfParticlesToSend_ = iterMap->second.second;

					net.RequestSendR(
							numberOfParticlesToSend_, neighbourRank);
					net.RequestReceiveR(
							numberOfParticlesToRecv_, neighbourRank);
				}
			}
			net.Dispatch();

			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
				scanMap[localRank].second += iterMap->second.first;
			particles.resize(scanMap[localRank].second);

			// prepare vector of particles for transmission
			// populate from migration list
			Particles.clear();
			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank = iterMap->first;
				if (neighbourRank != localRank) {
					for (unsigned long i = 0; i < iterMap->second.second; i++)
						Particles.push_back(particles[ghostOnProc.at(neighbourRank)[i]]);
				}
			}

			std::vector<Particle>::iterator iterSendBegin = Particles.begin();
			std::vector<Particle>::iterator iterRecvBegin = particles.begin() + particlesOffset_;

			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank = iterMap->first;
				if (neighbourRank != localRank) {

					const unsigned int& numberOfParticlesToRecv_ = iterMap->second.first;
					const unsigned int& numberOfParticlesToSend_ = iterMap->second.second;

					net.RequestSend(
							& ((PersistedParticle&) *iterSendBegin),
							numberOfParticlesToSend_, neighbourRank);
					net.RequestReceive(
							& ((PersistedParticle&) *iterRecvBegin),
							numberOfParticlesToRecv_, neighbourRank);

					iterRecvBegin += numberOfParticlesToRecv_;
					iterSendBegin += numberOfParticlesToSend_;
				}
			}
			net.Dispatch();

			int numParts = particles.size();

			// remove particles owned by unknown ranks
			std::vector<Particle>::iterator newEndOfParticles =
				std::partition(
						particles.begin(),
						particles.end(),
						std::bind2nd(std::mem_fun_ref(&Particle::IsOwnerRankKnown), scanMap));
			particles.erase(newEndOfParticles, particles.end());

			if (numParts != particles.size()) log::Logger::Log<log::Info, log::OnePerCore>("particles with unknown ranks have been deleted");

			// sort the particles - local first, then in order of increasing owner rank
			std::sort(particles.begin(), particles.end(), ParticleSorter(latDatLBM.GetLocalRank()));

			// re-build the scanMap (to know how many particles this rank owns)
			for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
				iterMap->second = scanMapElementType(0, 0);
			for (std::vector<Particle>::const_iterator iterParticles = particles.begin(); iterParticles != particles.end(); iterParticles++)
				scanMap[iterParticles->GetOwnerRank()].first++;
		}

		const void ParticleSet::CommunicateFluidVelocities(const LatticeTimeStep currentTimestep)
		{
			/** CommunicateFluidVelocities
			 */

			const unsigned int  particlesOffset  = scanMap[localRank].first;
			const unsigned int& particlesOffset_ = scanMap[localRank].first;

			if (scanMap.size() < 2)
			{
				return;
			}

			// exchange counts
			for (scanMapIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank = iterMap->first;
				if (neighbourRank != localRank)
				{
					unsigned int& numberOfVelocitiesToRecv_ = iterMap->second.second;
					unsigned int& numberOfVelocitiesToSend_ = iterMap->second.first;

					net.RequestSendR(
							numberOfVelocitiesToSend_, neighbourRank);
					net.RequestReceiveR(
							numberOfVelocitiesToRecv_, neighbourRank);
				}
			}
			net.Dispatch();

			// sum counts
			unsigned int numberOfIncomingVelocities = 0;
			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
				numberOfIncomingVelocities += iterMap->second.second;
			velocityBuffer.resize(numberOfIncomingVelocities);

			std::vector<std::pair<unsigned long, util::Vector3D<double> > >::iterator
				iterRecvBegin = velocityBuffer.begin();
			std::vector<Particle>::iterator
				iterSendBegin = particles.begin() + particlesOffset_;

			for (scanMapConstIterType iterMap = scanMap.begin(); iterMap != scanMap.end(); iterMap++)
			{
				const proc_t& neighbourRank = iterMap->first;
				if (neighbourRank != localRank) {

					const unsigned int& numberOfVelocitiesToRecv_ = iterMap->second.second;
					const unsigned int& numberOfVelocitiesToSend_ = iterMap->second.first;

					net.RequestSend(
							& ((Particle&) *iterSendBegin),
							numberOfVelocitiesToSend_, neighbourRank);
					net.RequestReceive(
							& *iterRecvBegin,
							numberOfVelocitiesToRecv_, neighbourRank);

					iterRecvBegin += numberOfVelocitiesToRecv_;
					iterSendBegin += numberOfVelocitiesToSend_;
				}
			}
			net.Dispatch();

			// sum velocities
			velocityMap.clear();
			for (std::vector<std::pair<unsigned long, util::Vector3D<double> > >::const_iterator iterVelocityBuffer = velocityBuffer.begin(); iterVelocityBuffer != velocityBuffer.end(); iterVelocityBuffer++)
			{
				const unsigned long& particleId = iterVelocityBuffer->first;
				const util::Vector3D<double>& partialVelocity = iterVelocityBuffer->second;
				velocityMap[particleId] += partialVelocity;
			}

			// update local particles
			for (std::vector<Particle>::iterator iter = particles.begin(); iter != particles.end(); iter++)
			{
				Particle& particle  = *iter;
				uint64_t particleId = particle.GetParticleId();
				if (particle.GetOwnerRank() == localRank && velocityMap.find(particleId) != velocityMap.end())
					particle.AccumulateVelocity(velocityMap.at(particleId));
			}
		}
	}
}
