
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "lb/kernels/momentBasis/DHumieresD3Q15MRTBasis.h"

namespace hemelb
{
  namespace lb
  {
    namespace kernels
    {
      namespace momentBasis
      {

        /*
         *  Kinetic moments defined in d'Humieres et al. 2002. To get the matrix below, columns 8 and 9 are respectively permuted
         *  with columns 14 and 11 to match HemeLB's lattice velocity ordering.
         *
         *  See publication for  meaning of e, epsilon, etc.
         */
        const distribn_t DHumieresD3Q15MRTBasis::REDUCED_MOMENT_BASIS[DHumieresD3Q15MRTBasis::NUM_KINETIC_MOMENTS][Lattice::NUMVECTORS] =
            { { -2, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, 1 }, // e
              { 16, -4, -4, -4, -4, -4, -4, 1, 1, 1, 1, 1, 1, 1, 1 }, // epsilon
              { 0, -4, 4, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1 }, // q_x
              { 0, 0, 0, -4, 4, 0, 0, 1, -1, 1, -1, -1, 1, -1, 1 }, // q_y
              { 0, 0, 0, 0, 0, -4, 4, 1, -1, -1, 1, 1, -1, -1, 1 }, // q_z
              { 0, 2, 2, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 }, // 3*p_xx
              { 0, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0 }, // p_ww = p_yy - p_zz
              { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, -1, -1, -1, -1 }, // p_xy
              { 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, -1, -1, 1, 1 }, // p_yz
              { 0, 0, 0, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1 }, // p_zx
              { 0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, -1, 1, 1, -1 } // m_xyz
            };

        const distribn_t DHumieresD3Q15MRTBasis::BASIS_TIMES_BASIS_TRANSPOSED[NUM_KINETIC_MOMENTS] =
            { 18., 360., 40., 40., 40., 12., 4., 8., 8., 8., 8. };

        void DHumieresD3Q15MRTBasis::SetUpCollisionMatrix(std::array<distribn_t, NUM_KINETIC_MOMENTS>& collisionMatrix,
                                                          distribn_t relaxationRate)
        {
          // Relaxation values taken from d'Humieres 2002.
          collisionMatrix.at(0) = 1.6; // e (s1)
          collisionMatrix.at(1) = 1.2; // epsilon (s2)
          collisionMatrix.at(2) = 1.6; // q_x (s4)
          collisionMatrix.at(3) = 1.6; // q_y (s4)
          collisionMatrix.at(4) = 1.6; // q_z (s4)
          collisionMatrix.at(5) = relaxationRate; // 3p_xx (s9)
          collisionMatrix.at(6) = relaxationRate; // p_ww (s9)
          collisionMatrix.at(7) = relaxationRate; // p_xy (s11 = s9)
          collisionMatrix.at(8) = relaxationRate; // p_yz (s11 = s9)
          collisionMatrix.at(9) = relaxationRate; // p_zx (s11 = s9)
          collisionMatrix.at(10) = 1.2; // m_xyz (s14)
        }

      }
    }
  }
}
