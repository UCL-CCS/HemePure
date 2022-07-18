
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LB_STREAMERS_GUOZHENGSHIELASTICWALLDELEGATE_H
#define HEMELB_LB_STREAMERS_GUOZHENGSHIELASTICWALLDELEGATE_H

#include "lb/streamers/BaseStreamerDelegate.h"
#include "util/Bessel.h"
#include <complex>

namespace hemelb
{
  namespace lb
  {
    namespace streamers
    {

      /**
       * This class implements the boundary condition described by Guo, Zheng and Shi
       * in 'An Extrapolation Method for Boundary Conditions in Lattice-Boltzmann method'
       * Physics of Fluids, 14/6, June 2002, pp 2007-2010.
       */
      template<typename CollisionImpl>
      class GuoZhengShiElasticWallDelegate : public BaseStreamerDelegate<CollisionImpl>
      {
        public:
          typedef CollisionImpl CollisionType;
          typedef typename CollisionType::CKernel::LatticeType LatticeType;

          GuoZhengShiElasticWallDelegate(CollisionType& delegatorCollider, kernels::InitParams& initParams) :
            collider(delegatorCollider),bValues(initParams.boundaryObject)
          {

          }

 
          inline void StreamLink(const LbmParameters* lbmParams,
                                 geometry::LatticeData* const latDat,
                                 const geometry::Site<geometry::LatticeData>& site,
                                 kernels::HydroVars<typename CollisionType::CKernel>& hydroVars,
                                 const Direction& iPrime)
          {
            Direction i = LatticeType::INVERSEDIRECTIONS[iPrime];
            // Get the distance to the boundary.
            double wallDistance = 1.0; //For elastic walls we are always assuming wall distance is 1

            // Set up for GZS - do the extrapolation from this site - u_ghost

            // Now we work out the hypothetical velocity of the ghost site on the other side
            // of the boundary.
            // Assume that the ghost velocity is linearly decreasing from the local fluid velocity 
	    // over a distance of dr - the amount of wall stretch at that location.
            
            // Hence velocityWall = velocityFluid * (dr-1)/dr
            distribn_t fWall[LatticeType::NUMVECTORS];
            kernels::HydroVars<typename CollisionType::CKernel> hydroVarsWall(fWall);
	   
	    distribn_t gamma = lbmParams->ElasticWallStiffness; //Wall stiffness 	    
	    
	    //LBM pressure version
	    //double dr = std::max((hydroVars.density - 1.0)/(3.0*gamma),0.0) + 1.0;
	    double dr = (hydroVars.density - 1.0)/(3.0*gamma) + 1.0;
	    
	    //Dynamic pressure version
	   // double dr = std::max(0.5*hydroVars.density*hydroVars.momentum.GetMagnitudeSquared()/gamma,0.0) + 1.0;

	    hydroVars.wallExtension = dr-1.0;

            hydroVarsWall.density = hydroVars.density;
            //hydroVarsWall.momentum = hydroVars.momentum * (dr - 1.)/dr;
            //hydroVarsWall.momentum = hydroVars.momentum * 0.43;
            
	    //double F = 0.525;
	    //double F = 0.85; //cylinder results
	    //double F = 0.025; //forearem arteries, dx=5e-5
	    //double F = 0.5; //forearem arteries, dx=2.1e-4
	    distribn_t F = lbmParams->BoundaryVelocityRatio; //Ratio of velocities at edge of simulation domain 	    
	    
	    if (dr > (1.0-F)){
		    hydroVarsWall.momentum = hydroVars.momentum * (F + dr - 1.)/dr; //Rejig from paper writeup July 2021
	    } else {
		    
		    hydroVarsWall.momentum = hydroVars.momentum * 0.0;
	    }

            // Find the non-equilibrium distribution in the unstreamed direction.
            std::copy(hydroVars.GetFNeqPtr(),
                      hydroVars.GetFNeqPtr() + LatticeType::NUMVECTORS,
                      hydroVarsWall.GetFNeqPtr());

            // Finally, we want to collide and stream, using the chosen collision kernel.
            //
            // hydroVarsWall contains the correct density and momentum, but
            // we must also set f, f_eq and f_neq to ensure that the collision
            // will work properly for more complex operators.

            // Calculate equilibrium values
            LatticeType::CalculateFeq(hydroVarsWall.density,
                                      hydroVarsWall.momentum.x,
                                      hydroVarsWall.momentum.y,
                                      hydroVarsWall.momentum.z,
                                      hydroVarsWall.GetFEqPtr());

            // For the wall site, construct f_old  = f_eq + f_neq
            for (unsigned j = 0; j < LatticeType::NUMVECTORS; ++j)
            {
              fWall[j] = hydroVarsWall.GetFEq()[j] + hydroVarsWall.GetFNeq()[j];
            }
            // Perform collision
            collider.Collide(lbmParams, hydroVarsWall);
            // stream
            distribn_t* fNew = latDat->GetFNew(site.GetIndex() * LatticeType::NUMVECTORS);
            fNew[i] = hydroVarsWall.GetFPostCollision()[i];

          }

        private:
	  typedef std::complex<double> Complex;          
          // the collision
          CollisionType collider;
	  iolets::BoundaryValues* bValues;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H */
