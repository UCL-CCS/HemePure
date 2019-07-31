
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_COLLOIDS_GRAVITICBODYFORCE_H
#define HEMELB_COLLOIDS_GRAVITICBODYFORCE_H

#include "colloids/BodyForces.h"

namespace hemelb
{
	namespace colloids
	{
		class GraviticBodyForce : public BodyForce
		{
			public:
				static BodyForce* ReadFromXml(io::xml::Element& xml)
				{
					LatticeForceVector field;
					io::xml::Element fieldElem = xml.GetChildOrThrow("field");
					fieldElem.GetAttributeOrThrow("x", field.x);
					fieldElem.GetAttributeOrThrow("y", field.y);
					fieldElem.GetAttributeOrThrow("z", field.z);

					return new GraviticBodyForce(field);
				};

				virtual const LatticeForceVector GetForceForParticle(const Particle& particle) const
				{
				  //return graviticAcceleration*particle.GetMass();
					return (graviticAcceleration*prefactor
						*particle.GetRadius(true)*particle.GetRadius(true)*particle.GetRadius(true)
						*(BLOOD_DENSITY_Kg_per_m3-MAGNE_DENSITY_Kg_per_m3));
				};

			protected:
				GraviticBodyForce(const LatticeForceVector constantForce) :
					graviticAcceleration(constantForce) {};

				const LatticeForceVector graviticAcceleration;
		};

		class GraviticBodyForceFactory : public BodyForceFactory<GraviticBodyForce> { };
	}
}
#endif /* HEMELB_COLLOIDS_GRAVITICBODYFORCE_H */
