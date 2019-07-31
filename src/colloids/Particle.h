
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_PARTICLE_H
#define HEMELB_COLLOIDS_PARTICLE_H

#include "net/mpi.h"
#include "colloids/PersistedParticle.h"
#include "geometry/LatticeData.h"
#include "io/xml/XmlAbstractionLayer.h"
#include "lb/MacroscopicPropertyCache.h"
#include "util/Vector3D.h"
#include "io/writers/Writer.h"

#include <set>
#include <time.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace hemelb
{
	namespace colloids
	{
		struct ParticleSorter;

		/**
		 * represents a single simulated biocolloid particle
		 *
		 * all persisted properties, i.e. those that are read in from a config file,
		 * are inherited from the PersistedParticle class (which handles the I/O)
		 */
		class Particle : PersistedParticle
		{
			public:
				/** constructor - gets initial values from an xml configuration file */
				Particle(const geometry::LatticeData& latDatLBM,
						const configuration::SimConfig *simConfig,
						io::xml::Element& xml);

				/** constructor - gets an invalid particle for making MPI data types */
				Particle() {};

				/** property getter for particleId */
				const unsigned long GetParticleId() const
				{
					return particleId;
				}

				/** property setter for particleId */
				const void SetParticleId(unsigned long id)
				{
					particleId = id;
				}

				/** property getter for ownerRank */
				const proc_t GetOwnerRank() const
				{
					return ownerRank;
				}

				/** property setter for ownerRank */
				const void SetOwnerRank(proc_t rank)
				{
					ownerRank = rank;
				}

				/** property getter for ownerRank */
				const proc_t GetCreationTimestep() const
				{
					return creationTimestep;
				}

				/** sets the timestep of creation */
				const void SetCreationTimestep(LatticeTimeStep timestep)
				{
					creationTimestep = timestep;
				}

				/** property getter for velocity */
				const LatticeVelocity GetVelocity() const
				{
					return velocity;
				}

				/** reset velocity */
				const void ResetVelocity()
				{
					velocity *= 0.0;
				}

				/** update velocity */
				const void UpdateVelocity()
				{
					velocity_ = velocity;
					velocity += bodyForce/GetDragCoefficient(false)+noise;
				}

				/** update velocity after applying boundary conditions */
				const void UpdateVelocity(const LatticeVelocity& x)
				{
					velocity += x;
				}

				/** property getter for bodyForce */
				const LatticeForceVector GetBodyForce() const
				{
					return bodyForce;
				}

				/** property getter for globalPosition */
				const LatticePosition GetGlobalPosition() const
				{
					return globalPosition;
				}

				/** property setter for globalPosition */
				const void SetGlobalPosition(LatticePosition position)
				{
					globalPosition = position;
				}

				/** property getter for smallRadius_a0 */
				const Dimensionless GetRadius(const bool units) const
				{
					if (units)	return (smallRadius_a0*voxelSize);
					else		return (smallRadius_a0);
				}

				const Dimensionless GetDragCoefficient(const bool units) const
				{
					// calculate the drag due to the fluid for this particle
					if (units)	return (6.0*PI*BLOOD_VISCOSITY_Pa_s
							*smallRadius_a0*voxelSize);
					else		return (6.0*PI*(hemelb::tau_-0.5)*Cs2
							*smallRadius_a0);
				}

				/** property getter for lastCheckpointTimestep */
				const LatticeTimeStep& GetLastCheckpointTimestep() const
				{
					return lastCheckpointTimestep;
				}

				/** property getter for markedForDeletionTimestep */
				const LatticeTimeStep& GetDeletionMarker() const
				{
					return markedForDeletionTimestep;
				}

				/** unsets the deletion marker - the particle will not be deleted */
				const void SetDeletionMarker()
				{
					markedForDeletionTimestep = SITE_OR_BLOCK_SOLID;
				}

				/** sets the deletion marker to the current timestep
				 *  the particle will be deleted after the next checkpoint
				 */
				const void SetDeletionMarker(LatticeTimeStep timestep)
				{
					if (timestep < markedForDeletionTimestep)
						markedForDeletionTimestep = timestep;
				}

				/** property getter for isValid */
				const bool IsValid() const
				{
					return isValid;
				}

				/** property getter for voxelSize */
				const PhysicalDistance GetVoxelSize() const
				{
					return voxelSize;
				}

				/**
				 * less than operator for comparing particle objects
				 *
				 * when used to sort a container of particle objects
				 * the ordering produced by this operator is:
				 * - increasing particleId
				 * - grouped by owner rank
				 * - with local rank first
				 */
				//const bool operator<(const Particle& other) const;

				/** determines if the owner rank of this particle is an existing key in map */
				const bool IsOwnerRankKnown(std::map<proc_t, std::pair<unsigned int, unsigned int> > map) const;

				const bool IsReadyToBeDeleted() const;

				/** for debug purposes only - outputs all properties to info log */
				const void OutputInformation() const;

				/** for serialisation into output file */
				const void WriteToStream(const LatticeTimeStep currentTimestep,
						io::writers::Writer& writer);

				/** calculates the drag coefficient */
				//const Dimensionless GetDragCoefficient(const bool units) const;

				/** updates the position of this particle using body forces and fluid velocity */
				const void UpdatePosition(const geometry::LatticeData& latDatLBM,
						const LatticeTimeStep currentTimestep,
						std::map<proc_t, std::vector<site_t> >& ghostOnProc,
						const unsigned int particleIdx);

				/** updates the position of this particle (first step) */
				const void UpdatePosition(const geometry::LatticeData& latDatLBM);

				/** calculates the effects of all body forces on this particle */
				const void CalculateBodyForces();

				/** interpolates the fluid velocity to the location of each particle */
				const void InterpolateFluidVelocity(const geometry::LatticeData& latDatLBM,
						const lb::MacroscopicPropertyCache& propertyCache);

				/** accumulate contributions to velocity from remote processes */
				const void AccumulateVelocity(util::Vector3D<double>& contribution)
				{
					velocity += contribution;
				}

				/** creates a derived MPI datatype that represents a single particle object
				 *  the fields included are all those from the PersistedParticle base class
				 *  note - this data type uses displacements rather than absolute addresses
				 *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
				 *  when you no longer need this type, remember to call MPI_Type_free
				 */
				const MPI_Datatype CreateMpiDatatypeWithPosition() const;

				/** creates a derived MPI datatype that represents a single particle object
				 *  the fields included in this type are: particleId and velocity(xyz) only
				 *  note - this data type uses displacements rather than absolute addresses
				 *  refer to Example 4.17 on pp114-117 of the MPI specification version 2.2
				 *  when you no longer need this type, remember to call MPI_Type_free
				 */
				const MPI_Datatype CreateMpiDatatypeWithVelocity() const;

			private:
				template<class T>
					double gen_norm(T &generator)
					{
						return generator();
					}

				/** rank owning this particle */
				proc_t ownerRank;

				/** velocity of this particle */
				LatticeVelocity velocity;
				LatticeVelocity velocity_;

				/** scaling of particle output */
				PhysicalDistance voxelSize;
				PhysicalTime timeStepLength;
				PhysicalPosition geometryOrigin;

				/** the effect of random thermal velocity contributions */
				LatticeVelocity noise;

				/** all body forces on this particle */
				LatticeForceVector bodyForce;

				/** is particle owned by known rank? */
				bool isValid;

				/** allow the sorter class to see our private members */
				friend struct ParticleSorter;

				/** convert velocity to lattice units */
				LatticeVelocity ConvertVelocityToLatticeUnits(const PhysicalVelocity& x) const
				{
					return x*timeStepLength/voxelSize;
				}

				/** convert velocity to physical units */
				PhysicalVelocity ConvertVelocityToPhysicalUnits(const LatticeVelocity& x) const
				{
					return x*voxelSize/timeStepLength;
				}

				/** convert force to lattice units */
				LatticeForceVector ConvertForceToLatticeUnits(const util::Vector3D<PhysicalForce>& x) const
				{
					return x*BLOOD_DENSITY_Kg_per_m3*pow(voxelSize,4)
						/(timeStepLength*timeStepLength);
				}
		};
	}

	namespace net
	{
		template<>
			MPI_Datatype MpiDataTypeTraits<colloids::Particle>::RegisterMpiDataType();
		template<>
			MPI_Datatype MpiDataTypeTraits<std::pair<unsigned long, util::Vector3D<double> > >::RegisterMpiDataType();
	}
}

#endif /* HEMELB_COLLOIDS_PARTICLE_H */
