/* whipplesteadycalcs.cpp
 * 
 * Copyright (C) 2010 Dale Lukas Peterson
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 3 of the License, or (at your option)
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
#include "whipplesteadycalcs.h"

void Whipple::writeBndryRecord_dt(const char * filename) const
{
  ofstream fp(filename, ios::out);
  if (fp.is_open()) {
    fp << "import numpy as np\n";
    fp << "boundary_dt = np.dtype([('q3', np.float64),\n"
          "                        ('q1_z', np.float64),\n"
          "                        ('q2_z', np.float64),\n"
          "                        ('q1_i', np.float64),\n"
          "                        ('q2_i', np.float64),\n"
          "                        ('q1_min', np.float64),\n"
          "                        ('q2_min', np.float64),\n"
          "                        ('q1_max', np.float64),\n"
          "                        ('q2_max', np.float64)])\n";
    fp.close();
  } else {
    cerr << "Unable to open " << filename << "for writing.\n";
    cerr << "Aborting.\n";
    exit(0);
  }
} // writeBndryRecord_dt()

// Process command line specified options to generate:
// 1) Boundary curves -- steady_boundary_curves.dat
// 2) Steady turning quantities within feasible region -- steady_feasible.dat
// 3) Specified velocity level curves -- steady_velocity.dat
// 4) Specified steer torque level curves -- steady_torque.dat
// 5) Specified friction coefficient curves -- steady_mew.dat
//
// Item 1 is done always, because these boundaries are used in the remaining
// calculations.  Items 2-5 are only done if whipplesteady is invoked with the
// correct command line options (see whipplesteady --help).
void Whipple::steadyCalcs(steadyOpts_t * opt)
{
  double delta = M_PI / double (opt->N - 1);
  int ig_index;
  string filename;
  gsl_matrix * M = gsl_matrix_alloc(opt->N, 9);
  gsl_vector * steer = &(gsl_matrix_column(M, 0).vector);
  for (int i = 0; i < opt->N; ++i)
    gsl_vector_set(steer, i, delta * i);
  gsl_vector * lean_zero = &(gsl_matrix_column(M, 1).vector);
  gsl_vector * pitch_zero = &(gsl_matrix_column(M, 2).vector);
  gsl_vector * lean_inf = &(gsl_matrix_column(M, 3).vector);
  gsl_vector * pitch_inf = &(gsl_matrix_column(M, 4).vector);
  gsl_vector * lean_min = &(gsl_matrix_column(M, 5).vector);
  gsl_vector * pitch_min = &(gsl_matrix_column(M, 6).vector);
  gsl_vector * lean_max = &(gsl_matrix_column(M, 7).vector);
  gsl_vector * pitch_max = &(gsl_matrix_column(M, 8).vector);
  
  // Find static equilibrium values of lean and pitch
  ig_index = staticEq(lean_zero, pitch_zero, steer, this);
  
  // Find infinite speed equilibrium values of lean and pitch
  infspeed(lean_inf, pitch_inf,
           gsl_vector_get(lean_zero, ig_index),
           gsl_vector_get(pitch_zero, ig_index),
           ig_index, steer, this);

  // Find configurational limits
  cfglim(lean_max, pitch_max, lean_min, pitch_min, steer, this);

  //if (opt->all)         // generate all quantities within feasible region
  //  steadyMesh(lean_zero, pitch_zero, lean_inf, pitch_inf,  // used for bounds
  //             steer, ig_index, opt);
  //if (opt->iso_v)       // generate iso velocity curves
  //  bb->steadyVel(opt);
  //if (opt->iso_t)       // generate iso torque curves
  //  bb->steadyTorque(opt);
  //if (opt->iso_mew)     // generate iso mew curves
  //  bb->steadyMew(opt);

  // Write boundary data to file
  filename = opt->outfolder; filename += "boundary.dat";
  FILE * fp = fopen(filename.c_str(), "wb");
  gsl_matrix_fwrite(fp, M);
  fclose(fp);
  // Write Numpy data type information to file
  filename = opt->outfolder; filename += "boundary_record.py";
  writeBndryRecord_dt(filename.c_str());

  gsl_matrix_free(M);
} // steadyCalcs()

static int static_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_vector_set(f, 0, bike->F[0]);      // Kinematic constraint equation
  gsl_vector_set(f, 1, bike->F[9]);      // Coefficient of g in steady lean eom
  return GSL_SUCCESS;
} // static_f()

static int static_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_matrix_set(J, 0, 0, bike->dF[0]);
  gsl_matrix_set(J, 0, 1, bike->dF[1]);
  gsl_matrix_set(J, 1, 0, bike->dF[27]);
  gsl_matrix_set(J, 1, 1, bike->dF[28]);
  return GSL_SUCCESS;
} // static_df()

static int static_fdf(const gsl_vector * x, void * params,
                      gsl_vector * f, gsl_matrix* J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_vector_set(f, 0, bike->F[0]);      // Kinematic constraint equation
  gsl_vector_set(f, 1, bike->F[9]);      // Coefficient of g in steady lean eom
  gsl_matrix_set(J, 0, 0, bike->dF[0]);
  gsl_matrix_set(J, 0, 1, bike->dF[1]);
  gsl_matrix_set(J, 1, 0, bike->dF[27]);      // Kinematic constraint equation
  gsl_matrix_set(J, 1, 1, bike->dF[28]);      // Coefficient of g in steady lean eom
  return GSL_SUCCESS;
} // static_fdf()

// Given a vector of steer values, calculate the lean values associated with
// static equilibrium.  Also, return the indices of the steer/lean vectors
// which most nearly cause the u5^2 coefficient to go to zero.
static int staticEq(gsl_vector * lean, gsl_vector * pitch,
             const gsl_vector * steer, Whipple * bike)
{
  int i, iter, max_iter = 50, status;
  double ftol = 1e-15;
  gsl_vector * x = gsl_vector_alloc(2);         // vector to store the solution
  gsl_vector * u5s_coefs = zeros(steer->size);
  const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, 2);
  gsl_multiroot_function_fdf f = {&static_f, &static_df, &static_fdf, 2, bike};
  gsl_vector_set(x, 0, bike->q1);
  gsl_vector_set(x, 1, bike->q2);
  gsl_multiroot_fdfsolver_set(s, &f, x);

  // for loop to loop over all values of steer
  for (i = 0; i < lean->size; ++i) {
    bike->q3 = gsl_vector_get(steer, i);  // steer as a parameter
    iter = 0;
    do
    {
      ++iter;
      status = gsl_multiroot_fdfsolver_iterate(s);
      if (status)
        iterateError(status, "staticEq()");
      status = gsl_multiroot_test_residual(s->f, ftol);
    } while (status == GSL_CONTINUE && iter < max_iter);
    // cout << "f = (" << s->f->data[0] << ", " << s->f->data[1] << ")\n";

    // Increase the tolerance by an order of magnitude to improve convergence
    //if (iter == max_iter) {
    //  increaseftol(&ftol, &i, "StaticEq()", bike->q3);
    //  continue;
    //} // if

    // Store the lean into the lean vector
    gsl_vector_set(lean, i, gsl_vector_get(s->x, 0));
    gsl_vector_set(pitch, i, gsl_vector_get(s->x, 1));

    // Store the square of the coefficient of the u5^2 term;
    gsl_vector_set(u5s_coefs, i, bike->F[10] * bike->F[10]);
    //ftol = 1e-15;  // reset the error tolerance
  } // for

  // Assign a large value to the u5s_coefs vector near steer = 0 and steer = PI
  // This ensure the minimum will be near PI/2 where the two boudary curves
  // cross
  for (i = 0; i < 5; ++i) {
    gsl_vector_set(u5s_coefs, i, 10000.0);
    gsl_vector_set(u5s_coefs, u5s_coefs->size - 1 - i, 10000.0);
  }

  // Free dynamically allocated variables
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
  i = gsl_vector_min_index(u5s_coefs);
  gsl_vector_free(u5s_coefs);
  return i;
} // staticEq()

static int cfglim_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_vector_set(f, 0, bike->F[0]);      // Kinematic constraint equation
  gsl_vector_set(f, 1, bike->dF[1]);     // Partial derivative of KC w.r.t. pitch
  return GSL_SUCCESS;
} // cfglim_f()

static int cfglim_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_matrix_set(J, 0, 0, bike->dF[0]);
  gsl_matrix_set(J, 0, 1, bike->dF[1]);
  gsl_matrix_set(J, 1, 0, bike->dF[33]);
  gsl_matrix_set(J, 1, 1, bike->dF[34]);
  return GSL_SUCCESS;
} // cfglim_df()

static int cfglim_fdf(const gsl_vector * x, void * params, gsl_vector * f,
                   gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_vector_set(f, 0, bike->F[0]);      // Kinematic constraint equation
  gsl_vector_set(f, 1, bike->dF[1]);     // Partial derivative of KC w.r.t. pitch
  gsl_matrix_set(J, 0, 0, bike->dF[0]);
  gsl_matrix_set(J, 0, 1, bike->dF[1]);
  gsl_matrix_set(J, 1, 0, bike->dF[33]);
  gsl_matrix_set(J, 1, 1, bike->dF[34]);
  return GSL_SUCCESS;
} // cfglim_fdf()

static void cfglim(gsl_vector * lean_max, gsl_vector * pitch_max,
                gsl_vector * lean_min, gsl_vector * pitch_min,
                const gsl_vector * steer, Whipple * bike)
{
  int i, N = steer->size, iter = 0, iter_max = 50, status;
  double ftol = 1e-14;
  gsl_vector * x = gsl_vector_alloc(2);         // vector to store the solution
  gsl_vector * lean, * pitch;
  const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, 2);
  gsl_multiroot_function_fdf f = {&cfglim_f, &cfglim_df,
                                  &cfglim_fdf, 2, bike};

  // Maximum lean initial guess
  gsl_vector_set(x, 0, M_PI/3.0);
  gsl_vector_set(x, 1, M_PI/2.0);
  lean = lean_max;            // set lean to point at max lean vector
  pitch = pitch_max;          // set pitch to point at max pitch vector
  for (int c = 0; c < 2; gsl_vector_set(x, 0, -M_PI/3.0), // min lean i.g.
                         gsl_vector_set(x, 1, M_PI/2.0),
                         lean = lean_min,     // point at min lean vector
                         pitch = pitch_min,   // point at min pitch vector
                         ++c) {
    gsl_multiroot_fdfsolver_set(s, &f, x);
    for (i = N / 2; i < N - 1; ++i) {
      bike->q3 = gsl_vector_get(steer, i);    // steer as a parameter
      iter = 0;
      do
      {
        ++iter;
        status = gsl_multiroot_fdfsolver_iterate(s);
        if (status)
          break;//iterateError(status, "cfglim()");
        status = gsl_multiroot_test_residual(s->f, ftol);
      } while(status == GSL_CONTINUE && iter < iter_max);
      //cout << "f = (" << s->f->data[0] << ", " << s->f->data[1] << ")\n";
      
      // Increase the tolerance by an order of magnitude to improve convergence
      //if (iter == iter_max) {
      //  increaseftol(&ftol, &i, "cfglim()", bike->q3);
      //  continue;
      //} // if

      // Store the lean and pitch
      gsl_vector_set(lean, i, gsl_vector_get(s->x, 0));
      gsl_vector_set(pitch, i, gsl_vector_get(s->x, 1));
      // Reset ftol in case it had been increased due to convergence issues
      //ftol = 1e-15;
    } // for i (steer from PI/2 to PI)
    gsl_vector_set(lean, i, gsl_vector_get(lean, i-1));
    gsl_vector_set(pitch, i, gsl_vector_get(pitch, i-1));

    gsl_vector_set(x, 0, gsl_vector_get(lean, N/2));
    gsl_vector_set(x, 1, gsl_vector_get(pitch, N/2));
    gsl_multiroot_fdfsolver_set(s, &f, x);
    for (i = N / 2 - 1; i > 0; --i) {
      bike->q3 = gsl_vector_get(steer, i);  // steer as a parameter
      iter = 0;
      do
      {
        ++iter;
        status = gsl_multiroot_fdfsolver_iterate(s);
        if (status)
          break;//iterateError(status, "cfglim()");
        status = gsl_multiroot_test_residual(s->f, ftol);
      } while(status == GSL_CONTINUE && iter < iter_max);
      //cout << "f = (" << s->f->data[0] << ", " << s->f->data[1] << ")\n";
      
      // Increase the tolerance by an order of magnitude to improve convergence
      //if (iter == iter_max) {
      //  increaseftol(&ftol, &i, "cfglim()", bike->q3);
      //  continue;
      //} // if

      // Store the lean and pitch
      gsl_vector_set(lean, i, gsl_vector_get(s->x, 0));
      gsl_vector_set(pitch, i, gsl_vector_get(s->x, 1));
      // Reset ftol in case it had been increased due to convergence issues
      //ftol = 1e-15;
    } // for i (steer from PI/2 to O)
    gsl_vector_set(lean, 0, gsl_vector_get(lean, 1));
    gsl_vector_set(pitch, 0, gsl_vector_get(pitch, 1));
  } // for c

  // Free dynamically allocated variables
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
} // cfglim()

static int inf_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_vector_set(f, 0, bike->F[0]);      // Kinematic constraint equation
  gsl_vector_set(f, 1, bike->F[10]);    // Coefficient of u5^2 in steady eom
  return GSL_SUCCESS;
} // inf_f()

static int inf_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_matrix_set(J, 0, 0, bike->dF[0]);
  gsl_matrix_set(J, 0, 1, bike->dF[1]);
  gsl_matrix_set(J, 1, 0, bike->dF[30]);
  gsl_matrix_set(J, 1, 1, bike->dF[31]);
  return GSL_SUCCESS;
} // inf_df()

static int inf_fdf(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = gsl_vector_get(x, 0);       // Lean
  bike->q2 = gsl_vector_get(x, 1);       // Pitch
  bike->steadyEqns();
  gsl_vector_set(f, 0, bike->F[0]);      // Kinematic constraint equation
  gsl_vector_set(f, 1, bike->F[10]);     // Coefficient of u5^2 in steady eom
  gsl_matrix_set(J, 0, 0, bike->dF[0]);
  gsl_matrix_set(J, 0, 1, bike->dF[1]);
  gsl_matrix_set(J, 1, 0, bike->dF[30]);
  gsl_matrix_set(J, 1, 1, bike->dF[31]);
  return GSL_SUCCESS;
} // inf_fdf()

static void infspeed(gsl_vector * lean, gsl_vector * pitch,
                    double lean_ig, double pitch_ig, int ig_index,
                    const gsl_vector * steer, Whipple * bike)
{
  // We need 
  int i, iter, status, iter_max = 100;
  double ftol = 1e-14;

  gsl_vector * x = gsl_vector_alloc(2);         // vector to store the solution
  const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, 2);
  gsl_multiroot_function_fdf f = {&inf_f, &inf_df, &inf_fdf, 2, bike};

  // Setup the initial conditions
  bike->q1 = lean_ig;
  bike->q2 = pitch_ig;
  bike->q3 = gsl_vector_get(steer, ig_index);
  bike->calcPitch();
  gsl_vector_set(x, 0, lean_ig);
  gsl_vector_set(x, 1, bike->q2);
  gsl_multiroot_fdfsolver_set(s, &f, x);

  for (i = ig_index; i > 0; --i) {
    bike->q3 = gsl_vector_get(steer, i);
    iter = 0;

    do {
      status = gsl_multiroot_fdfsolver_iterate(s);
      if (status)
        break;
      status = gsl_multiroot_test_residual(s->f, ftol);
    } while(status == GSL_CONTINUE && ++iter < iter_max);
    // Increase the tolerance by an order of magnitude to improve convergence
    //if (iter == iter_max) {
    //  gsl_vector_set(x, 0, gsl_vector_get(lean, i-1));
    //  gsl_vector_set(x, 1, gsl_vector_get(pitch, i-1));
    //  gsl_multiroot_fdfsolver_set(s, &f, x);
    //  increaseftol(&ftol, &i, "infspeed()", bike->q3);
    //  continue;
    //} // if
    
    // Store the lean into the lean vector
    gsl_vector_set(lean, i, gsl_vector_get(s->x, 0));
    gsl_vector_set(pitch, i, gsl_vector_get(s->x, 1));
  } // for
  gsl_vector_set(lean, i, gsl_vector_get(lean, 1));
  gsl_vector_set(pitch, i, gsl_vector_get(pitch, 1));
  
  // Setup the initial conditions
  bike->q1 = lean_ig;
  bike->q2 = pitch_ig;
  bike->q3 = gsl_vector_get(steer, ig_index);
  bike->calcPitch();
  gsl_vector_set(x, 0, lean_ig);
  gsl_vector_set(x, 1, bike->q2);
  gsl_multiroot_fdfsolver_set(s, &f, x);

  for (i = ig_index + 1; i < steer->size - 1; ++i) {
    bike->q3 = gsl_vector_get(steer, i);
    
    iter = 0;
    do {
      status = gsl_multiroot_fdfsolver_iterate(s);

      if (status)
        break;

      status = gsl_multiroot_test_residual(s->f, ftol);
    } while (status == GSL_CONTINUE && ++iter < iter_max);
    // Increase the tolerance by an order of magnitude to improve convergence
    //if (iter == iter_max) {
    //  gsl_vector_set(x, 0, gsl_vector_get(lean, i-1));
    //  gsl_vector_set(x, 1, gsl_vector_get(pitch, i-1));
    //  gsl_multiroot_fdfsolver_set(s, &f, x);
    //  increaseftol(&ftol, &i, "infspeed()", bike->q3);
    //  continue;
    //} // if
    
    // Store the lean into the lean vector
    gsl_vector_set(lean, i, gsl_vector_get(s->x, 0));
    gsl_vector_set(pitch, i, gsl_vector_get(s->x, 1));
  } // for
  gsl_vector_set(lean, i, gsl_vector_get(lean, i - 1));
  gsl_vector_set(pitch, i, gsl_vector_get(pitch, i - 1));
  
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
}

static void iterateError(int status, const char * routine)
{
  if (status == GSL_EBADFUNC)
    cerr << "GSL_EBADFUN";
  else if (status == GSL_ENOPROG)
    cerr << "GSL_ENOPROG";
  cerr << " encountered in " << routine << ".  Aborting.\n";
  exit(0);
}

static void increaseftol(double * ftol, int * i,
                         const char * routine, double steer)
{
  (*ftol) = (*ftol) * 10.0;
  (*i) = (*i) - 1;
  cerr << "Max iterations encountered in " << routine << " at steer = " << steer
       << "\nIncreasing tolerance to ftol = " << *ftol << '\n';
} // increaseftol()
