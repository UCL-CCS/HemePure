
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.
#ifndef HEMELB_LB_IOLETS_INOUTLETREADWRITEVELOCITY_H
#define HEMELB_LB_IOLETS_INOUTLETREADWRITEVELOCITY_H

#include "lb/iolets/InOutLetVelocity.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      class InOutLetReadWriteVelocity : public InOutLetVelocity
      {
        public:
          InOutLetReadWriteVelocity();

          InOutLet* Clone() const;

          bool IsCommsRequired() const
          {
            return true;
          }

          void DoComms(const BoundaryCommunicator& boundaryComm, const LatticeTimeStep timeStep);

          void SetWeightsFilePath(const std::string& path)
          {
            weightsFilePath = path;
          }

          void SetArea(const distribn_t& a)
          {
            area = a;
          }

          const LatticeTimeStep& GetCouplingFrequency()
          {
            return couplingFrequency;
          }
          void SetCouplingFrequency(const LatticeTimeStep& frequency)
          {
            couplingFrequency = frequency;
          }

          const std::string& GetFlowRateFilePath()
          {
            return flowRateFilePath;
          }
          void SetFlowRateFilePath(const std::string& path)
          {
            flowRateFilePath = path;
          }

          const std::string& GetPressureFilePath()
          {
            return pressureFilePath;
          }
          void SetPressureFilePath(const std::string& path)
          {
            pressureFilePath = path;
          }

          const double GetFlowRateConversionFactor()
          {
            return flowRateConversionFactor;
          }
          void SetFlowRateConversionFactor(const double& factor)
          {
            flowRateConversionFactor = factor;
          }

          const double GetPressureConversionFactor()
          {
            return pressureConversionFactor;
          }
          void SetPressureConversionFactor(const double& factor)
          {
            pressureConversionFactor = factor;
          }

          PhysicalSpeed ConvertFlowRateToVelocity(const distribn_t flowRate)
          {
            return flowRate / weights_sum;
          }

          void SetSmoothingFactor(const double& factor)
          {
            smoothingFactor = factor;
          }

          void SetWarmup(const LatticeTimeStep& warmup)
          {
            warmUpLength = warmup;
          }

          LatticeVelocity GetVelocity(const LatticePosition& x, const LatticeTimeStep t) const;

          void Initialise(const util::UnitConverter* unitConverter);

          bool useWeightsFromFile;

          void DoPreStreamCoupling(const site_t& siteID,
                                   const LatticeTimeStep& timeStep,
                                   const LatticeVector& sitePos,
                                   const LatticeDensity& density,
                                   const LatticeVelocity& velocity);

          void DoPostStreamCoupling(const site_t& siteID,
                                    const LatticeTimeStep& timeStep,
                                    const LatticeVector& sitePos);

        private:
          std::string flowRateFilePath, pressureFilePath;
          std::string weightsFilePath;
          std::map<std::vector<int>, double> weights_table;
          const util::UnitConverter* units;
          distribn_t area, densitySum, densityAvg;
          site_t siteCount;
          PhysicalTime startTime;
          LatticeTimeStep warmUpLength, couplingTimeStep, couplingFrequency;
          double weights_sum, flowRateConversionFactor, pressureConversionFactor, smoothingFactor;
          LatticeSpeed maxVelocity, maxVelocityNew;
      };

    }
  }
}
#endif /* HEMELB_LB_IOLETS_INOUTLETREADWRITEVELOCITY_H */
