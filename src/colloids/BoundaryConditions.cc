
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "colloids/BoundaryConditions.h"
#include "colloids/BoundaryConditionDeletion.h"
#include "colloids/BoundaryConditionLubrication.h"
#include "colloids/BoundaryConditionSpherical.h"
#include "geometry/Site.h"
#include "geometry/SiteData.h"

namespace hemelb
{
	namespace colloids
	{
		std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsWall;
		std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsIlet;
		std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsOlet;
		std::vector<BoundaryCondition*> BoundaryConditions::boundaryConditionsSphr;

		const geometry::LatticeData* BoundaryConditions::latticeData;

		const void BoundaryConditions::InitBoundaryConditions(
				const geometry::LatticeData* const latticeData,
				io::xml::Document& xml)
		{
			BoundaryConditions::latticeData = latticeData;

			std::map<std::string, BoundaryConditionFactory_Create> mapBCGenerators;
#ifndef HEMELB_TRACER_PARTICLES
			mapBCGenerators["lubrication"] = &(
					BoundaryConditionFactoryLubrication::Create);
#endif
			mapBCGenerators[   "deletion"] = &(
					BoundaryConditionFactoryDeletion::Create);
			mapBCGenerators[  "spherical"] = &(
					BoundaryConditionFactorySpherical::Create);

			io::xml::Element colloidsBC = xml.GetRoot().GetChildOrThrow("colloids").GetChildOrThrow("boundaryConditions");

			for (std::map<std::string, BoundaryConditionFactory_Create>::const_iterator
					iter  = mapBCGenerators.begin();
					iter != mapBCGenerators.end();
					iter++)
			{
				const std::string boundaryConditionClass = iter->first;
				const BoundaryConditionFactory_Create createFunction = iter->second;
				for(	// there must be at least one element for each type
						io::xml::Element bcNode = colloidsBC.GetChildOrThrow(boundaryConditionClass);
						!(bcNode == io::xml::Element::Missing());
						bcNode = bcNode.NextSiblingOrNull(boundaryConditionClass))
				{
					const std::string& appliesTo = bcNode.GetAttributeOrThrow("appliesTo");
					BoundaryCondition* nextBC = createFunction(bcNode);
					if		(appliesTo == "Wall")
						BoundaryConditions::boundaryConditionsWall.push_back(nextBC);
					else if (appliesTo == "Ilet")
						BoundaryConditions::boundaryConditionsIlet.push_back(nextBC);
					else if (appliesTo == "Olet")
						BoundaryConditions::boundaryConditionsOlet.push_back(nextBC);
					else if (appliesTo == "Sphr" && boundaryConditionClass == "spherical")
						BoundaryConditions::boundaryConditionsSphr.push_back(nextBC);
				}
			}
		}

		const bool BoundaryConditions::DoSomeThingsToParticle(
				const LatticeTimeStep currentTimestep,
				Particle& particle)
		{
			bool keep = true;

			// detect collision(s)
			proc_t procId;
			site_t siteId;
			const util::Vector3D<site_t> siteGlobalPosition(
					(site_t)(0.5+particle.GetGlobalPosition().x),
					(site_t)(0.5+particle.GetGlobalPosition().y),
					(site_t)(0.5+particle.GetGlobalPosition().z));
			const bool isValid = latticeData->GetContiguousSiteId(
					siteGlobalPosition, procId, siteId);
			if (!isValid) {return keep;}

			const lb::lattices::LatticeInfo latticeInfo = BoundaryConditions::latticeData->GetLatticeInfo();
			const geometry::Site<const geometry::LatticeData> site = latticeData->GetSite(siteId);
			const geometry::SiteData siteData = site.GetSiteData();
			const geometry::SiteType siteType = siteData.GetSiteType();
			const distribn_t* siteWallDistances = site.GetWallDistances();

			const bool isNearWall = siteData.IsWall();
			const bool isNearIlet = (siteType == geometry::INLET_TYPE);
			const bool isNearOlet = (siteType == geometry::OUTLET_TYPE);

			// if the particle is not near a boundary then simply keep it
			if (!isNearWall && !isNearIlet && !isNearOlet) {
				particle.SetDeletionMarker();
				return keep;
			}

			log::Logger::Log<log::Trace, log::OnePerCore>(
					"*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu, isNearWall: %s, isNearIlet: %s, isNearOlet: %s ***\n",
					particle.GetParticleId(),
					isNearWall ? "TRUE" : "FALSE",
					isNearIlet ? "TRUE" : "FALSE",
					isNearOlet ? "TRUE" : "FALSE");

			std::vector<LatticePosition> particleToWallVectors;
			// only use lattice vectors 1 to 6 (the face-of-a-cube vectors)
			for (Direction direction = 1; direction <= 6; ++direction)
			{
				// in general, this "distance" is a fraction of a non-unit lattice vector
				// however, we treat this fractional magnitude as a real lattice distance
				// because all of the face-of-a-cube lattice vectors will be unit vectors
				double thisDistance = siteWallDistances[direction - 1];

				// a negative distance to the wall from a site in any direction means that
				// the wall is further away than the nearest fluid site in that direction
				if (thisDistance < 0.0) continue;

				// the particle cannot be allowed to go past halfway between this site and
				// the next lattice site in this direction, because the next site is solid
				// the wall is assumed to be no further away than half the distance to the
				// solid site so that the particle never becomes nearest to a solid site
				if (thisDistance > 0.5) thisDistance = 0.5;

				// conversion from LatticeCoordinate to LatticePosition is done
				// auto-magically by the multiplication & its arithmetic traits
				const LatticePosition siteToWall = latticeInfo.GetVector(direction) * thisDistance;
				const LatticePosition particleToSite = (LatticePosition)siteGlobalPosition
					- particle.GetGlobalPosition();

				// particleToWall = siteToWall + projection of particleToSite in the siteToWall direction
				const LatticePosition particleToWallVector = siteToWall +
					siteToWall.GetNormalised() * siteToWall.GetNormalised().Dot(particleToSite);

				log::Logger::Log<log::Trace, log::OnePerCore>(
						"*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu, siteToWall: {%g,%g,%g}, particleToSite: {%g,%g,%g}, particleToWall: {%g,%g,%g}\n",
						particle.GetParticleId(),
						siteToWall.x, siteToWall.y, siteToWall.z,
						particleToSite.x, particleToSite.y, particleToSite.z,
						particleToWallVector.x, particleToWallVector.y, particleToWallVector.z);

				particleToWallVectors.push_back(particleToWallVector);
			}

			if (isNearWall)
				for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsWall.begin();
						iter != boundaryConditionsWall.end(); iter++)
				{
					BoundaryCondition& boundaryCondition = **(iter);
					keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
				}

			if (isNearIlet)
				for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsIlet.begin();
						iter != boundaryConditionsIlet.end(); iter++)
				{
					BoundaryCondition& boundaryCondition = **(iter);
					keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
				}

			if (isNearOlet)
				for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsOlet.begin();
						iter != boundaryConditionsOlet.end(); iter++)
				{
					BoundaryCondition& boundaryCondition = **(iter);
					keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
				}

			for (std::vector<BoundaryCondition*>::iterator iter = boundaryConditionsSphr.begin();
					iter != boundaryConditionsSphr.end(); iter++)
			{
				BoundaryCondition& boundaryCondition = **(iter);
				keep &= boundaryCondition.DoSomethingToParticle(particle, particleToWallVectors);
			}

			if (keep) particle.SetDeletionMarker();
			else {
				particle.SetDeletionMarker(currentTimestep);
				log::Logger::Log<log::Trace, log::OnePerCore>(
						"*** In BoundaryConditions::DoSomeThingsToParticle for id: %lu - attempting to set markedForDeletion to %lu (value actually becomes: %lu)\n",
						particle.GetParticleId(),
						currentTimestep,
						particle.GetDeletionMarker());
			} return keep;
		}
	}
}
