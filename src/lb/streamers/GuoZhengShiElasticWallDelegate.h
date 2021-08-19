
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

	    //Womersley Pressure version
	    //Data to import
	    //do similar to gamma
	    //Ensure in LU
	    double womersleyNumber = 3.544907701811032;
	    double period = 6000.0;
	    double cO = 1.178511301977579; //sqrt(thick*E/2*rho)
	    double poisson = 0.5;
	    double pressureGrad = -7.406795112344784e-07;
	  	    
	    omega = 2.0 * PI / period;
	    
	    Complex lambda = pow(Complex(0.0,1.0),1.5) * womersleyNumber;
            Complex g = 2.0 * util::BesselJ1ComplexArgument(lambda) / ( lambda * util::BesselJ0ComplexArgument(lambda));

	    //Build terms of frequency equation, we will ASSUME that the wall thickness is 0.1*Radius and rho_wall == rho_fluid
	    //Quadratic equation of form Ax^2 + Bx + C = 0
	    double thick = 0.1;
	    Complex A = (g - 1.0) * (poisson * poisson - 1);
       	    Complex B = thick * (g - 1.0) + (2.0 * poisson - 0.5) * g - 2.0;
	    Complex C = 2.0*thick + g;
    
	    Complex rt1 = (-B + sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
	    Complex rt2 = (-B - sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
    
	    Complex nu;
            if ((rt1.real() > rt2.real()) || ((rt1.real() == rt2.real()) && (rt1.imag() > rt2.imag())))
	    {
		    nu = rt1;
	    }
	    else
	    {
		    nu = rt2;
	    }
    
	    Complex M = (2.0 + nu * (2.0 * poisson - 1.0)) / (nu * (2.0 * poisson - g));

	    Complex c = cO*sqrt(2.0/((1.0-poisson)*(1.0-poisson)*nu));
	    Complex H = abs(Complex(0.0,1.0)*c*pressureGrad/omega);

	    double mod1 = abs(H)*cos(std::arg(H)) - abs(H*(1.0-M)/c)*cos(std::arg(H*(1.0-M)/c));
	    double mod2 = abs(H)*sin(std::arg(H)) - abs(H*(1.0-M)/c)*sin(std::arg(H*(1.0-M)/c));
	    Pm = sqrt(mod1*mod1 + mod2*mod2);
	    Pt = atan(mod2/mod1);
	    if (mod1<0){
		Pt += PI;
	    }
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
	   
	    double gamma = lbmParams->ElasticWallStiffness; //Wall stiffness 	    
	    
	    //LBM pressure version
	    //double dr = std::max((hydroVars.density - 1.0)/(3.0*gamma),0.0) + 1.0;
	    double dr = (hydroVars.density - 1.0)/(3.0*gamma) + 1.0;
	    
	    //Dynamic pressure version
	   // double dr = std::max(0.5*hydroVars.density*hydroVars.momentum.GetMagnitudeSquared()/gamma,0.0) + 1.0;
	    


	    const distribn_t* neighbourFOld;
	    util::Vector3D<float> ioletNormal = LatticeVector(1,0,0);
	    
	    // Find the neighbour's global location and which proc it's on.
	    LatticeVector neighbourGlobalLocation = site.GetGlobalSiteCoords() + LatticeVector(lrint(ioletNormal.x), lrint(ioletNormal.y), lrint(ioletNormal.z));
	    //LatticeVector neighbourGlobalLocation = site.GetGlobalSiteCoords() + LatticeVector(0,0,-1);
	    //std::cout<< "neighbour at " << neighbourGlobalLocation.x << "," << neighbourGlobalLocation.y << "," << neighbourGlobalLocation.z << std::endl;
	    proc_t neighbourProcessor = latDat->GetProcIdFromGlobalCoords(neighbourGlobalLocation);
	    //std::cout << "On Proc " << neighbourProcessor << std::endl;
	    
	    LatticeVelocity neighbourMomentum;
	    if (neighbourProcessor == latDat->GetLocalRank())
	    {
		      // If it's local, get a Site object for it.
		      geometry::Site<geometry::LatticeData> nextSiteOut =
			  latDat->GetSite(latDat->GetContiguousSiteId(neighbourGlobalLocation));
		      neighbourFOld = nextSiteOut.GetFOld<LatticeType> ();

		    
		    LatticeVelocity neighbourVelocity;
		    distribn_t neighbourFEq[LatticeType::NUMVECTORS];
			  // Go ahead and calculate the density, momentum and eqm distribution.
		    distribn_t neighbourDensity;
		    //LatticeVelocity neighbourMomentum;
			    // Note that nextNodeOutVelocity is passed as the momentum argument, this
			    // is because it is immediately divided by density when the function returns.
		    LatticeType::CalculateDensityMomentumFEq(neighbourFOld,
							     neighbourDensity,
							     neighbourMomentum.x,
							     neighbourMomentum.y,
							     neighbourMomentum.z,
							     neighbourVelocity.x,
							     neighbourVelocity.y,
							     neighbourVelocity.z,
							     neighbourFEq);
	    }

	    hydroVars.coverageFactor = dr-1.0;

            hydroVarsWall.density = hydroVars.density;
            //hydroVarsWall.momentum = hydroVars.momentum * (dr - 1.)/dr;
            //hydroVarsWall.momentum = hydroVars.momentum * 0.43;
            
	    //double F = 0.525;
	    //double F = 0.85; //cylinder results
	    //double F = 0.025; //forearem arteries, dx=5e-5
	    double F = 0.5; //forearem arteries, dx=2.1e-4
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

	double Pt, Pm, omega;
      };

    }
  }
}

#endif /* HEMELB_LB_STREAMERS_GUOZHENGSHIDELEGATE_H */
