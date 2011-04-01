/*
 * Copyright (C) 2011 Dale Lukas Peterson
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * Ensures that eigenvalues of linearized dynamic equations match the values
 * published in Meijaard et al. 2007.
 * Every eigenvalue they published in Table 2 of that paper is hand coded as
 * written in that paper and compared with the values calculated by the model
 * which makes use of the gyrostat formulation of the equations of motion.
 * Most eigenvalues match to within 1e-14, the larger eigenvalues only match to
 * about 1e-12 due to the nature of floating point comparisons.
 * Each eigenvalue in their table has 14 digits presented after the decimal
 * place.  This isn't ideal since the number of significant digits at most can
 * be not quite 16 for IEEE double precision numbers.  For example some of
 * their eigenvalues have two digits to the left of the decimal place, while
 * others have none -- so they are effectively reporting results with different
 * amounts of numerical precision.  Ideally, one would use standard floating
 * point notation where one digit is to the left of the decimal place, and 15
 * are after it, along with the exponent.
 */
#include "whipple.h"
#include "whippleutils.h"
#include "gtest/gtest.h"

TEST(Whipple, Eigenvalues) {
  Whipple * bb = new Whipple();
  int is_complex[4] = {0, 0, 0, 0};
  int N_real = 0, N_complex = 0;

  double evals[] = {-3.13164324790656, 0.0,
                        -3.13423125066578, 0.0,
                        -3.07158645641514, 0.0,
                        -2.63366137253667, 0.0,
                        -1.42944427361326, 0.0,
                        -0.32286642900409, 0.0,
                        -0.00406690076970, 0.0,
                        0.10268170574766, 0.0,
                        0.14327879765713, 0.0,
                        0.15790184030917, 0.0,
                        0.16105338653172, 0.0,
                        3.13164324790656, 0.0,
                        3.52696170990070, -0.8077402751993,
                        2.68234517512745, -1.68066296590675,
                        1.70675605663975, -2.31582447384325,
                        0.41325331521125, -3.07910818603206,
                        -0.77534188219585, -4.46486771378823,
                        -1.52644486584142, -5.87673060598709,
                        -2.13875644258362, -7.19525913329805,
                        -2.69348683581097, -8.46037971396931,
                        -3.21675402252485, -9.69377351531791,
                        -3.72016840437287, -10.90681139476287,
                        -5.53094371765393, 0.0,
                        3.52696170990070, +0.8077402751993,
                        2.68234517512745, +1.68066296590675,
                        1.70675605663975, +2.31582447384325,
                        0.41325331521125, +3.07910818603206,
                        -0.77534188219585, +4.46486771378823,
                        -1.52644486584142, +5.87673060598709,
                        -2.13875644258362, +7.19525913329805,
                        -2.69348683581097, +8.46037971396931,
                        -3.21675402252485, +9.69377351531791,
                        -3.72016840437287, +10.90681139476287,
                        5.53094371765393, 0.0,
                        -7.11008014637442, 0.0,
                        -8.67387984831735, 0.0,
                        -10.35101467245920, 0.0,
                        -12.15861426576450, 0.0,
                        -14.07838969279820, 0.0,
                        -16.08537123098030, 0.0,
                        -18.15788466125262, 0.0,
                        -20.27940894394569, 0.0,
                        -22.43788559040858, 0.0,
                        -24.62459635017404, 0.0};

  gsl_matrix_complex_view m = gsl_matrix_complex_view_array(evals, 4, 11);

  bb->u1 = 0.0;
  bb->u3 = 0.0;
  bb->evalConstants();

  // Values from Table 2 part (a)
  // bifurcation speed, weave speed, capsize speed
  double speeds[] = {0.68428307889246, // bifurcation speed
                     4.29238253634111, // weave speed
                     6.02426201538837};// capsize speed
  for (int j = 0; j < 3; ++j) {
    bb->u5 = -speeds[j]/(bb->rf+bb->rft);
    bb->calcEvals();
    if (j == 0) {
      ASSERT_NEAR(
          GSL_REAL(gsl_vector_complex_get(bb->evals, 1)),
          3.78290405129320,
          1.0e-14);
      ASSERT_NEAR(
            GSL_REAL(gsl_vector_complex_get(bb->evals, 2)),
            3.78290405129320,
            1.0e14);
      ASSERT_NEAR(
            GSL_IMAG(gsl_vector_complex_get(bb->evals, 1)),
            0.0,
            1e-7); // not exactly zero imaginary part
      ASSERT_NEAR(
            GSL_IMAG(gsl_vector_complex_get(bb->evals, 2)),
            0.0,
            1.0e-7); // not exactly real
    }
    if (j == 1) {
      ASSERT_TRUE((fabs(
              GSL_IMAG(gsl_vector_complex_get(bb->evals, 1))
              - 3.43503384866144) <= 1e-14)
             ^
             (fabs(
              GSL_IMAG(gsl_vector_complex_get(bb->evals, 1))
              + 3.43503384866144) <= 1e-14));
      ASSERT_TRUE((fabs(
              GSL_IMAG(gsl_vector_complex_get(bb->evals, 2))
              - 3.43503384866144) <= 1e-14)
             ^
             (fabs(
              GSL_IMAG(gsl_vector_complex_get(bb->evals, 2))
              + 3.43503384866144) <= 1e-14));
      ASSERT_TRUE(fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 1)))
             <= 1e-14); // not purely imag
      ASSERT_TRUE(fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 2)))
             <= 1e-14); // not purely imag
    }
    if (j == 2) {
      ASSERT_TRUE(fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 0)))
          <= 1e-14);
      ASSERT_TRUE(fabs(GSL_IMAG(gsl_vector_complex_get(bb->evals, 0)))
          <= 1e-14);
    }
  } // for j

  // Verifying values in Table 2 parts (b) and (c)
  for (int j = 0; j < 11; ++j) {
    bb->u5 = -j/(bb->rf+bb->rft);
    bb->calcEvals();
    if (j == 0) {
        // assert that all are real
        for (int i = 0; i < 4; ++i) {
          ASSERT_TRUE(GSL_IMAG(gsl_vector_complex_get(bb->evals, j)) == 0.0);
        } // for i
        // Now match the magnitudes in Meijaard2007
        ASSERT_TRUE((fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 0)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 0, j)))
                <= 1e-14)
               ^
                (fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 0)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 1, j)))
                <= 1e-14)
               );
        ASSERT_TRUE((fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 1)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 0, j)))
                <= 1e-14)
               ^
                (fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 1)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 1, j)))
                <= 1e-14)
               );
        ASSERT_TRUE((fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 2)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 2, j)))
                <= 1e-14)
               ^
                (fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 2)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 3, j)))
                <= 1e-14)
               );
        ASSERT_TRUE((fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 3)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 2, j)))
                <= 1e-14)
               ^
                (fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 3)) -
                GSL_REAL(gsl_matrix_complex_get(&(m.matrix), 3, j)))
                <= 1e-14)
               );
        // Make sure they are of opposite sign
        ASSERT_NE(signbit(GSL_REAL(gsl_vector_complex_get(bb->evals, 0))),
                  signbit(GSL_REAL(gsl_vector_complex_get(bb->evals, 1))));
        ASSERT_NE(signbit(GSL_REAL(gsl_vector_complex_get(bb->evals, 2))),
                  signbit(GSL_REAL(gsl_vector_complex_get(bb->evals, 3))));
    } else {
      is_complex[0] = is_complex[1] = is_complex[2] = is_complex[3] = 0;
      N_real = N_complex = 0;
      for (int i = 0; i < 4; ++i) {
        if (GSL_IMAG(gsl_vector_complex_get(bb->evals, i)) != 0.0) {
          is_complex[i] = 1;
          N_complex += 1;
        } else
          N_real += 1;
      } // for i
      // Assert that we have two real and two complex
      ASSERT_EQ(N_real, 2);
      ASSERT_EQ(N_complex, 2);

      // Eigenvalues are sorted by complex magnitude from smallest to
      // largest, so if we detect a complex eig, we know that the next one is
      // its complex conjugate, but we need to skip the next eigenvalue
      // Unfortunately, the order of the complex eigenvalues can have either
      // the negative imaginary one first or second.
      for (int i = 0; i < 4; ++i) {
        if (i == 0) {
          ASSERT_NEAR(GSL_REAL(gsl_vector_complex_get(bb->evals, i)),
                      GSL_REAL(gsl_matrix_complex_get(&(m.matrix), i, j)),
                      1.0e-13); // all passed with 1e-13, not 1e-14
        }
        if (i == 1) {
          // double check that they are complex conjugates, probably a redudant
          // test
          // all passed
          ASSERT_EQ(GSL_REAL(gsl_vector_complex_get(bb->evals, i)),
                    GSL_REAL(gsl_vector_complex_get(bb->evals, i+1)));
          ASSERT_EQ(GSL_IMAG(gsl_vector_complex_get(bb->evals, i)),
                    -GSL_IMAG(gsl_vector_complex_get(bb->evals, i+1)));

          // Checking that the real parts are the same
          // all pass
          ASSERT_NEAR(
                 GSL_REAL(gsl_vector_complex_get(bb->evals, i)),
                 GSL_REAL(gsl_matrix_complex_get(&(m.matrix), i, j)),
                 1e-13);  // all pass with 1e-13, not 1e-14
          ASSERT_NEAR(
                 GSL_REAL(gsl_vector_complex_get(bb->evals, i)),
                 GSL_REAL(gsl_matrix_complex_get(&(m.matrix), i+1, j)),
                 1e-13);  // all pass with 1e-13, not 1e-14
          ASSERT_NEAR(
                 GSL_REAL(gsl_vector_complex_get(bb->evals, i+1)),
                 GSL_REAL(gsl_matrix_complex_get(&(m.matrix), i, j)),
                 1e-13);  // all pass with 1e-13, not 1e-14
          ASSERT_NEAR(
                 GSL_REAL(gsl_vector_complex_get(bb->evals, i+1)),
                 GSL_REAL(gsl_matrix_complex_get(&(m.matrix), i+1, j)),
                 1e-13);  // all pass with 1e-13, not 1e-14

          // Checking that the imaginary parts are the same
          ASSERT_TRUE((fabs(
                 GSL_IMAG(gsl_vector_complex_get(bb->evals, i)) -
                 GSL_IMAG(gsl_matrix_complex_get(&(m.matrix), i, j)))
                 <= 1e-13)
                 ^
                 (fabs(
                 GSL_IMAG(gsl_vector_complex_get(bb->evals, i)) -
                 GSL_IMAG(gsl_matrix_complex_get(&(m.matrix), i+1, j)))
                 <= 1e-13));  // all pass with 1e-13, not 1e-14
          ASSERT_TRUE((fabs(
                 GSL_IMAG(gsl_vector_complex_get(bb->evals, i+1)) -
                 GSL_IMAG(gsl_matrix_complex_get(&(m.matrix), i, j)))
                 <= 1e-13)
                 ^
                 (fabs(
                 GSL_IMAG(gsl_vector_complex_get(bb->evals, i+1)) -
                 GSL_IMAG(gsl_matrix_complex_get(&(m.matrix), i+1, j)))
                 <= 1e-13));  // all pass with 1e-13, not 1e-14
          i += 1; // skip to i=3
        } // i == 1, 2
        if (i == 3) {
          ASSERT_NEAR(GSL_REAL(gsl_vector_complex_get(bb->evals, i)),
                      GSL_REAL(gsl_matrix_complex_get(&(m.matrix), i, j)),
                      1e-12); // all passed with 1e-12, not 1e-13
        } // i == 3
      } // for i
    } // else
  } // for j

  // Close files and free memory
  delete bb;
}
