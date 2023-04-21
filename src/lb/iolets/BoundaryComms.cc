
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/iolets/BoundaryComms.h"
#include "lb/iolets/BoundaryValues.h"
#include "net/IOCommunicator.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace lb
  {
    namespace iolets
    {

      BoundaryComms::BoundaryComms(SimulationState* iSimState, int centreRank, const BoundaryCommunicator& boundaryComm) :
          nProcs(boundaryComm.Size()), bcComm(boundaryComm)
      {
        // Let centreRank be the BoundaryControlling(BC) proc
        bcComm.SetBCProcRank(centreRank);

        // Only BC proc sends
        if (bcComm.IsCurrentProcTheBCProc())
        {
          sendRequest = new MPI_Request[nProcs];
          sendStatus = new MPI_Status[nProcs];
        }
        else
        {
          sendRequest = NULL;
          sendStatus = NULL;
        }
      }

      BoundaryComms::~BoundaryComms()
      {

        if (bcComm.IsCurrentProcTheBCProc())
        {
          delete[] sendRequest;
          delete[] sendStatus;
        }
      }

      void BoundaryComms::Send(distribn_t* density)
      {
        if (bcComm.IsCurrentProcTheBCProc())
        {
          for (int proc = 0; proc < nProcs; proc++)
          {
            HEMELB_MPI_CALL(
                MPI_Isend, (
                    density,
                    1,
                    net::MpiDataType(*density),
                    proc,
                    100,
                    bcComm,
                    &sendRequest[proc]
                ));
          }
        }
      }

      void BoundaryComms::Receive(distribn_t* density)
      {
        HEMELB_MPI_CALL(
            MPI_Irecv, (
                density,
                1,
                net::MpiDataType(*density),
                bcComm.GetBCProcRank(),
                100,
                bcComm,
                &receiveRequest
            ));
      }

      void BoundaryComms::FinishSend()
      {
        // Don't move on to next step with BC proc until all messages have been sent
        // Precautionary measure to make sure proc doesn't overwrite, before message is sent
        if (bcComm.IsCurrentProcTheBCProc())
        {
          HEMELB_MPI_CALL(
              MPI_Waitall, (nProcs, sendRequest, sendStatus)
          );
        }
      }

      void BoundaryComms::Wait()
      {
        HEMELB_MPI_CALL(
            MPI_Wait, (&receiveRequest, &receiveStatus)
        );
      }

      void BoundaryComms::WaitAllComms()
      {
        FinishSend();
        Wait();
      }

    }
  }
}
