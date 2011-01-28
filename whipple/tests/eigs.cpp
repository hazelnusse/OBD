/* whippleeig.cpp
 * 
 * Copyright (C) 2010 Dale Lukas Peterson
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
#include <assert.h>
#include "whipple.h"
#include "whippleutils.h"

/*
meijaard_evals = [5.53094371765393 + 0.00000000000000i,
                  3.52696170990070 + 0.807740275199300i,
                  2.68234517512745 + 1.68066296590675i,
                  1.70675605663975 + 2.31582447384325i,
                  0.413253315211250 + 3.07910818603206i,-0.775341882195850 + 4.46486771378823i,-1.52644486584142 + 5.87673060598709i,-2.13875644258362 + 7.19525913329805i,-2.69348683581097 + 8.46037971396931i,-3.21675402252485 + 9.69377351531791i,-3.72016840437287 + 10.9068113947629i;
3.13164324790656 + 0.00000000000000i,3.52696170990070 - 0.807740275199300i,2.68234517512745 - 1.68066296590675i,1.70675605663975 - 2.31582447384325i,0.413253315211250 - 3.07910818603206i,-0.775341882195850 - 4.46486771378823i,-1.52644486584142 - 5.87673060598709i,-2.13875644258362 - 7.19525913329805i,-2.69348683581097 - 8.46037971396931i,-3.21675402252485 - 9.69377351531791i,-3.72016840437287 - 10.9068113947629i;
-3.13164324790656 + 0.00000000000000i,-3.13423125066578 + 0.00000000000000i,-3.07158645641514 + 0.00000000000000i,-2.63366137253667 + 0.00000000000000i,-1.42944427361326 + 0.00000000000000i,-0.322866429004090 + 0.00000000000000i,-0.00406690076970000 + 0.00000000000000i,0.102681705747660 + 0.00000000000000i,0.143278797657130 + 0.00000000000000i,0.157901840309170 + 0.00000000000000i,0.161053386531720 + 0.00000000000000i;
-5.53094371765393 + 0.00000000000000i,-7.11008014637442 + 0.00000000000000i,-8.67387984831735 + 0.00000000000000i,-10.3510146724592 + 0.00000000000000i,-12.1586142657645 + 0.00000000000000i,-14.0783896927982 + 0.00000000000000i,-16.0853712309803 + 0.00000000000000i,-18.1578846612526 + 0.00000000000000i,-20.2794089439457 + 0.00000000000000i,-22.4378855904086 + 0.00000000000000i,-24.6245963501740 + 0.00000000000000i;];
*/

int main(int argc, char ** argv)
{
  Whipple * bb = new Whipple();

  bb->u1 = 0.0;
  bb->u3 = 0.0;
  bb->evalConstants();

  for (int i = 0; i < 11; ++i) {
    bb->u5 = -i/(bb->rf+bb->rft);
    bb->calcEvals();
    cout << "v = " << i << "\n";
    bb->printEvals();
    switch (i) {
      case 0:
        // assert that all are real
        for (int j = 0; j < 4; ++j) {
          assert(GSL_IMAG(gsl_vector_complex_get(bb->evals, j)) == 0.0);
        } // for j
        // Now match the magnitudes in Meijaard2007
        assert(fabs(fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 0))) -
               3.13164324790656) <= 1e-14);
        assert(fabs(fabs(GSL_REAL(gsl_vector_complex_get(bb->evals, 1))) -
               3.13164324790656) <= 1e-14);
        // Make sure they are of opposite sign
        assert(signbit(GSL_REAL(gsl_vector_complex_get(bb->evals, 0))) !=
               signbit(GSL_REAL(gsl_vector_complex_get(bb->evals, 1))));

        break;
      case 1:
        // assert that all are complex
        for (int j = 0; j < 4; ++j) {
          assert(GSL_REAL(gsl_vector_complex_get(bb->evals, j)) != 0.0);
          assert(GSL_IMAG(gsl_vector_complex_get(bb->evals, j)) != 0.0);
        } // for j
        break;

      default:
        break;
    } // switch()
  } // for i

  // Close files and free memory
  delete bb;
  return 0;
} // main

