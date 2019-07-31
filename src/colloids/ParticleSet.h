
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_PARTICLESET_H
#define HEMELB_COLLOIDS_PARTICLESET_H

#include <vector>
#include <list>

#include "geometry/LatticeData.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "net/mpi.h"
#include "colloids/Particle.h"
#include "net/IOCommunicator.h"
#include "units.h"

namespace hemelb
{
	namespace colloids
	{
		/** represents the set of all particles known to the local process */
		class ParticleSet
		{
			public:
				/** constructor - gets local particle information from xml config file */
				ParticleSet(geometry::LatticeData& latDatLBM,
						io::xml::Element& xml,
						lb::MacroscopicPropertyCache& propertyCache,
						//const hemelb::lb::LbmParameters *lbmParams,
						const configuration::SimConfig* simConfig,
						std::vector<proc_t>& neighbourProcessors,
						const net::IOCommunicator& ioComms_,
						const std::string& outputPath);

				/** destructor - de-allocates all Particle objects created by this Set */
				~ParticleSet();

				/** updates the position of each particle */
				const void UpdatePositions(const LatticeTimeStep currentTimestep);

				/** updates the velocity of each particle using body forces */
				const void UpdateVelocities();

				/** resets the velocity of each particle */
				const void ResetVelocities();

				/** calculates the effect of all body forces on each particle */
				const void CalculateBodyForces();

				/** applies boundary conditions to all particles **/
				const void ApplyBoundaryConditions(const LatticeTimeStep currentTimestep);

				/** inject some delicious particles **/
				const void EmitParticles(const LatticeTimeStep currentTimestep);

				/** interpolates the fluid velocity to the location of each particle */
				const void InterpolateFluidVelocity(const LatticeTimeStep currentTimestep);

				/** communicates the positions of all particles to&from all neighbours */
				const void CommunicateParticlePositions(const LatticeTimeStep currentTimestep);

				/** communicates the partial fluid interpolations to&from all neighbours */
				const void CommunicateFluidVelocities(const LatticeTimeStep currentTimestep);

				const void OutputInformation(const LatticeTimeStep currentTimestep);

			private:
				const net::IOCommunicator& ioComms;

				/** cached copy of local rank (obtained from topology) */
				const proc_t localRank;

				/** voxel size */
				double voxelSize;

				/** time step length */
				double timeStepLength;

				/**
				 * contains all particles known to this process
				 * they are sorted using the less than operator
				 */
				std::vector<Particle> particles;

				/** contains particles to be sent to neighbouring rank */
				std::vector<Particle> Particles;

				/** particle to be replicated */
				Particle masterParticle;

				/** total number of particles read from file */
				unsigned long particlesFromFile;

				/** total number of particles emitted */
				unsigned long tEmitted = 0;

				/** number of particles to be inserted during step */
				unsigned long nEmitted = 0;

				/** defining sphere where replicated particles are introduced */
				LatticePosition sphereCentre;
				LatticeDistance sphereRadius;

				/** map neighbourRank -> {numberOfParticlesFromThere, numberOfVelocitiesFromThere} */
				typedef std::pair<unsigned int, unsigned int> scanMapElementType;
				typedef std::map <proc_t, scanMapElementType>::const_iterator scanMapConstIterType;
				typedef std::map <proc_t, scanMapElementType>::iterator scanMapIterType;
				typedef std::pair<proc_t, scanMapElementType> scanMapContentType;
				std::map<proc_t, scanMapElementType> scanMap;

				/** contiguous buffer into which MPI can write all the velocities from neighbours */
				std::vector<std::pair<unsigned long, util::Vector3D<double> > > velocityBuffer;

				/** map particleId -> sumOfvelocityContributionsFromNeighbours */
				std::map<unsigned long, util::Vector3D<double> > velocityMap;

				/** contains useful geometry manipulation functions */
				//const geometry::LatticeData& latDatLBM;
				geometry::LatticeData& latDatLBM;

				/**
				 * primary mechanism for interacting with the LB simulation
				 * - the velocity cache  : is used for velocity interpolation
				 * - the bodyForce cache : stores the colloid feedback forces
				 */
				lb::MacroscopicPropertyCache& propertyCache;

				/**
				 * abstracts communication via MPI.
				 */
				net::Net net;

				/**
				 * reusable output buffer.
				 */
				std::vector<char> buffer;

				/**
				 * path to write to.
				 */
				const std::string& path;

				/**
				 * MPI File handle to write with
				 */
				net::MpiFile file;

				/**
				 * particle migration list
				 */
				std::map<proc_t, std::vector<site_t> > ghostOnProc;
		};
	}
}

#endif /* HEMELB_COLLOIDS_PARTICLESET_H */
