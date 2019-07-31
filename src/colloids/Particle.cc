
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/Particle.h"
#include "colloids/BodyForces.h"
#include "geometry/LatticeData.h"
#include "log/Logger.h"
#include "units.h"

namespace hemelb
{
	namespace colloids
	{
		Particle::Particle(const geometry::LatticeData& latDatLBM,
				const configuration::SimConfig *simConfig,
				io::xml::Element& xml) :
			PersistedParticle(xml)
		{
			// updating position with zero velocity and zero body force is necessary
			// because of the side-effect that sets owner rank from the new position
			ownerRank  = SITE_OR_BLOCK_SOLID;
			velocity  *= 0.0;
			bodyForce *= 0.0;

			voxelSize      = simConfig->GetVoxelSize();
			timeStepLength = simConfig->GetTimeStepLength();
			geometryOrigin = simConfig->GetGeometryOrigin();

			UpdatePosition(latDatLBM);
		  //OutputInformation();
		}

		const bool Particle::IsOwnerRankKnown(
				std::map<proc_t, std::pair<unsigned int, unsigned int> > map) const
		{
			return map.count(ownerRank) > 0;
		}

		const bool Particle::IsReadyToBeDeleted() const
		{
			return markedForDeletionTimestep < lastCheckpointTimestep;
		}

		const void Particle::OutputInformation() const
		{
			log::Logger::Log<log::Info, log::OnePerCore>("In Colloids::Particle::OutputInformation, ID: %0.6i, owner: %0.7i, position: {%+4.2e,%+4.2e,%+4.2e}, velocity: {%+4.2e,%+4.2e,%+4.2e}",//bodyforce: {%g,%g,%g}",
					particleId, ownerRank,
					globalPosition.x, globalPosition.y, globalPosition.z,
					velocity.x, velocity.y, velocity.z);
					//bodyForce.x, bodyForce.y, bodyForce.z);
		}

		// this is the exact size that the xdr data produced for this particle will occupy
		// 5 fields * 8 bytes-per-field = 40 bytes
		// if more fields are included, change io::formats::colloids::RecordLength accordingly

		const void Particle::WriteToStream(const LatticeTimeStep currentTimestep,
				io::writers::Writer& writer)
		{
			lastCheckpointTimestep = currentTimestep;

			writer << (uint64_t) particleId;
			writer << (uint64_t) ownerRank;
			writer << globalPosition.x << globalPosition.y << globalPosition.z;
			writer << velocity.x << velocity.y << velocity.z;
			writer << bodyForce.x << bodyForce.y << bodyForce.z;
		}

		const void Particle::UpdatePosition(const geometry::LatticeData& latDatLBM)
		{
			// round the global position of the particle to the nearest site coordinates
			const util::Vector3D<site_t> siteGlobalPosition(
					(site_t)(0.5 + globalPosition.x),
					(site_t)(0.5 + globalPosition.y),
					(site_t)(0.5 + globalPosition.z));

			// convert the site coordinates into a local site index and find owner rank
			proc_t procId = latDatLBM.GetProcIdFromGlobalCoords(siteGlobalPosition);
			isValid = (procId != SITE_OR_BLOCK_SOLID);

			if (isValid && (ownerRank != procId))
			{
				log::Logger::Log<log::Trace, log::OnePerCore>("Changing owner of particle %i from %i to %i at step %i - %s",
						particleId, ownerRank, procId, -1, isValid ? "valid" : "INVALID");
				ownerRank = procId;
			}

			if (log::Logger::ShouldDisplay<log::Trace>())
				log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::UpdatePosition, id: %i, position is now: {%g,%g,%g}\n",
						particleId, globalPosition.x, globalPosition.y, globalPosition.z);
		}

		const void Particle::UpdatePosition(const geometry::LatticeData& latDatLBM,
				const LatticeTimeStep currentTimestep,
				std::map<proc_t, std::vector<site_t> >& ghostOnProc,
				const unsigned int particleIdx)
		{
			// first, update the position: newPosition = oldPosition + velocity + bodyForce/drag
			// then,  update the owner rank for the particle based on its new position
			if (currentTimestep > creationTimestep+2)
				globalPosition += velocity;

			// round the global position of the particle to the nearest site coordinates
			const util::Vector3D<site_t> siteGlobalPosition(
					(site_t)(0.5 + globalPosition.x),
					(site_t)(0.5 + globalPosition.y),
					(site_t)(0.5 + globalPosition.z));

			// convert the site coordinates into a local site index and find owner rank
			proc_t procId = latDatLBM.GetProcIdFromGlobalCoords(siteGlobalPosition);
			isValid = (procId != SITE_OR_BLOCK_SOLID);

			// this is only necessary if particles are stuck or leave the geometry
			// a deficiency caused by the lubrication boundary
#ifndef HEMELB_TRACER_PARTICLES
			// if site is solid, move particle only due to flow
			if (!isValid) {
				// revert position...
				globalPosition -= velocity;
				// ...and update only due to flow
				globalPosition += velocity_;

				// round the global position of the particle to the nearest site coordinates
				const util::Vector3D<site_t> siteGlobalPosition(
						(site_t)(0.5 + globalPosition.x),
						(site_t)(0.5 + globalPosition.y),
						(site_t)(0.5 + globalPosition.z));

				// convert the site coordinates into a local site index and find owner rank
				proc_t procId = latDatLBM.GetProcIdFromGlobalCoords(siteGlobalPosition);
				isValid = (procId != SITE_OR_BLOCK_SOLID);

				// still not valid!
				if (!isValid)
					log::Logger::Log<log::Trace, log::OnePerCore>("Particle %i from %i to %i at step %i - still not valid",
							particleId, ownerRank, procId, currentTimestep);
			}
#endif

			if (isValid && (ownerRank != procId))
			{
				log::Logger::Log<log::Trace, log::OnePerCore>("Changing owner of particle %i from %i to %i at step %i - %s",
						particleId, ownerRank, procId, currentTimestep, isValid ? "valid" : "INVALID");
				ownerRank = procId;
				ghostOnProc[procId].push_back(particleIdx);
			}

			if (log::Logger::ShouldDisplay<log::Trace>())
				log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::UpdatePosition, id: %i, position is now: {%g,%g,%g}\n",
						particleId, globalPosition.x, globalPosition.y, globalPosition.z);

			// check if particle is close to rank boundary (and needs to be communicated)
			std::set<proc_t> procsFound;
			// determine the global coordinates of the next neighbour site
			for (site_t x = ((site_t)globalPosition.x)-1;
					x < ((site_t)globalPosition.x)+3; x++)
				for (site_t y = ((site_t)globalPosition.y)-1;
						y < ((site_t)globalPosition.y)+3; y++)
					for (site_t z = ((site_t)globalPosition.z)-1;
							z < ((site_t)globalPosition.z)+3; z++)
					{
						const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

						// convert the global coordinates of the site into a local site index
						//site_t siteId;
						//bool isSiteValid = latDatLBM.GetContiguousSiteId(siteGlobalPosition, procId, siteId);
						//bool isSiteLocal = (procId == latDatLBM.GetLocalRank());

						proc_t procId = latDatLBM.GetProcIdFromGlobalCoords(siteGlobalPosition);
						bool isSiteLocal = (procId == latDatLBM.GetLocalRank());
						if (!isSiteLocal) procsFound.insert(procId);
					}

			for (std::set<proc_t>::const_iterator iter = procsFound.begin(); iter != procsFound.end(); iter++)
				ghostOnProc[*iter].push_back(particleIdx);
		}

		const void Particle::CalculateBodyForces()
		{
			log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateBodyForces, id: %i, position: {%g,%g,%g}\n",
					particleId, globalPosition.x, globalPosition.y, globalPosition.z);

			// delegate the calculation of body force to the BodyForces class
			bodyForce = ConvertForceToLatticeUnits(
					BodyForces::GetBodyForcesForParticle(*this));

			log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::CalculateBodyForces, id: %i, position: {%g,%g,%g}, bodyForce: {%g,%g,%g}",
					particleId, globalPosition.x, globalPosition.y, globalPosition.z,
					bodyForce.x, bodyForce.y, bodyForce.z);
		}

		/** modified dirac delta function according to Peskin */
		const Dimensionless diracOperation(const LatticePosition& relativePosition)
		{
			Dimensionless delta = 1.0;
			for (int xyz = 0; xyz < 3; xyz++)
			{
				const LatticeDistance rmod = fabs(relativePosition[xyz]);

				if (rmod <= 1.0)
					delta *= 0.125 * (3.0 - 2.0 * rmod + sqrt( 1.0 +  4.0 * rmod - 4.0 * rmod * rmod));
				else if (rmod <= 2.0)
					delta *= 0.125 * (5.0 - 2.0 * rmod - sqrt(-7.0 + 12.0 * rmod - 4.0 * rmod * rmod));
				else
					delta = 0.0;
			}
			return delta;
		}

		const void Particle::InterpolateFluidVelocity(const geometry::LatticeData& latDatLBM,
				const lb::MacroscopicPropertyCache& propertyCache)
		{
			/** InterpolateFluidVelocity
			 *    For each local neighbour lattice site
			 *    - get velocity for the site from the macroscopic cache object
			 *    - calculate site's contribution to the velocity interpolation
			 *    - increment particle velocity property with this contribution
			 *    - will require communication to transmit remote contributions
			 */

			log::Logger::Log<log::Debug, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, id: %i, position: {%g,%g,%g}\n",
					particleId, globalPosition.x, globalPosition.y, globalPosition.z);

			velocity *= 0.0;
			// determine the global coordinates of the next neighbour site:
			// nested loop - x, y, z directions semi-open interval [-2, +2)
			for (site_t x = ((site_t)globalPosition.x)-1;
					x < ((site_t)globalPosition.x)+3; x++)
				for (site_t y = ((site_t)globalPosition.y)-1;
						y < ((site_t)globalPosition.y)+3; y++)
					for (site_t z = ((site_t)globalPosition.z)-1;
							z < ((site_t)globalPosition.z)+3; z++)
					{
						const util::Vector3D<site_t> siteGlobalPosition(x, y, z);

						// convert the global coordinates of the site into a local site index
						proc_t procId;
						site_t siteId;
						bool isSiteValid = latDatLBM.GetContiguousSiteId(siteGlobalPosition, procId, siteId);
						bool isSiteLocal = (procId == latDatLBM.GetLocalRank());

						/** TODO: implement boundary conditions for invalid/solid sites */
						if (!isSiteValid || !isSiteLocal)
							continue;

						// read value of velocity for site index from macroscopic cache
						LatticeVelocity siteFluidVelocity = propertyCache.velocityCache.Get(siteId);

						// calculate term of the interpolation sum
						LatticePosition relativePosition(siteGlobalPosition);
						relativePosition -= globalPosition;
						LatticeVelocity partialInterpolation = siteFluidVelocity * diracOperation(relativePosition);

						// accumulate each term of the interpolation
						velocity += partialInterpolation;

						log::Logger::Log<log::Trace, log::OnePerCore>("In colloids::Particle::InterpolateFluidVelocity, particleId: %i, siteIndex: %i, fluidVelocity: {%g,%g,%g}, partialInterpolation: {%g,%g,%g}, velocitySoFar: {%g,%g,%g}\n",
								particleId, siteId, siteFluidVelocity.x, siteFluidVelocity.y, siteFluidVelocity.z,
								partialInterpolation.x, partialInterpolation.y, partialInterpolation.z,
								velocity.x, velocity.y, velocity.z);
					}
		}
	}
}
