#include "whipplesteadyboundaries.h"

void Whipple::steadyBoundaries(steadyOpts_t * options)
{
  size_t N = 1001;   // number of steer angles to look at
  int crossing_ig_index;
  gsl_vector * steer = linspaceN(0.0, M_PI, N); // Steer angles
  gsl_vector * lean_zero = zeros(N);            // Zero speed lean angles
  gsl_vector * pitch_zero = zeros(N);           // Zero speed lean angles
  gsl_vector * lean_min = zeros(N);             // Min lean config lean
  gsl_vector * pitch_min = zeros(N);            // Min lean config pitch
  gsl_vector * lean_max = zeros(N);             // Max lean config lean
  gsl_vector * pitch_max = zeros(N);            // Max lean config pitch
  
  // Find static equilibrium values of lean and pitch
  // This is the zero-speed boundary of the steady turning region
  crossing_ig_index = staticEq(lean_zero, pitch_zero, steer, this);
  // write_results("./plots/zero_speed_steady.dat", lean_zero, steer);
  
  // configuration limit boundary curves
  cfglim(lean_max, pitch_max, lean_min, pitch_min, steer, this);
  //write_results("./plots/min_lean.dat", lean_min, steer);
  //write_results("./plots/max_lean.dat", lean_max, steer);
  

  gsl_vector_free(steer);
  gsl_vector_free(lean_zero);
  gsl_vector_free(pitch_zero);
  gsl_vector_free(lean_min);
  gsl_vector_free(pitch_min);
  gsl_vector_free(lean_max);
  gsl_vector_free(pitch_max);
} // steadyBoundaries()

static int static_f(const gsl_vector * x, void * params, gsl_vector * f)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = x->data[0];       // Lean
  bike->q2 = x->data[1];       // Pitch
  bike->steadyEqns();
  f->data[0] = bike->F[0];      // Kinematic constraint equation
  f->data[1] = bike->F[9];      // Coefficient of g in steady lean eom
  return GSL_SUCCESS;
} // static_f()

static int static_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = x->data[0];       // Lean
  bike->q2 = x->data[1];       // Pitch
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
  bike->q1 = x->data[0];       // Lean
  bike->q2 = x->data[1];       // Pitch
  bike->steadyEqns();
  f->data[0] = bike->F[0];      // Kinematic constraint equation
  f->data[1] = bike->F[9];      // Coefficient of g in steady lean eom
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
  x->data[0] = bike->q1;
  x->data[1] = bike->q2;
  gsl_multiroot_fdfsolver_set(s, &f, x);

  // for loop to loop over all values of steer
  for (i = 0; i < lean->size; ++i) {
    bike->q3 = steer->data[i];  // steer as a parameter
    iter = 0;
    do
    {
      status = gsl_multiroot_fdfsolver_iterate(s);
      if (status)
        iterateError(status, "staticEq()");
      status = gsl_multiroot_test_residual(s->f, ftol);
    } while (status == GSL_CONTINUE && ++iter < max_iter);

    // Increase the tolerance by an order of magnitude to improve convergence
    if (iter == max_iter) {
      increaseftol(&ftol, &i, "StaticEq()", bike->q3);
      continue;
    } // if

    // Store the lean into the lean vector
    lean->data[i] = s->x->data[0];
    pitch->data[i] = s->x->data[1];

    // Store the square of the coefficient of the u5^2 term;
    u5s_coefs->data[i] = bike->F[10] * bike->F[10];
    ftol = 1e-15;  // reset the error tolerance
  } // for

  // Assign a large value to the u5s_coefs vector near steer = 0 and steer = PI
  // This ensure the minimum will be near PI/2 where the two boudary curves
  // cross
  for (i = 0; i < 5; ++i)
    u5s_coefs->data[i] = u5s_coefs->data[lean->size - 1 - i] = 10000.0;

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
  bike->q1 = x->data[0];       // Lean
  bike->q2 = x->data[1];       // Pitch
  bike->steadyEqns();
  f->data[0] = bike->F[0];      // Kinematic constraint
  f->data[1] = bike->dF[1];     // Partial derivative of KC w.r.t. pitch
  return GSL_SUCCESS;
} // cfglim_f()

static int cfglim_df(const gsl_vector * x, void * params, gsl_matrix * J)
{
  Whipple * bike = (Whipple *) params;
  bike->q1 = x->data[0];       // Lean
  bike->q2 = x->data[1];       // Pitch
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
  bike->q1 = x->data[0];       // Lean
  bike->q2 = x->data[1];       // Pitch
  bike->steadyEqns();
  f->data[0] = bike->F[0];     // Kinematic constraint
  f->data[1] = bike->dF[1];    // Partial derivative of KC w.r.t. pitch
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
  double ftol = 1e-15;
  gsl_vector * x = gsl_vector_alloc(2);         // vector to store the solution
  gsl_vector * lean, * pitch;
  const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_newton;
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(T, 2);
  gsl_multiroot_function_fdf f = {&cfglim_f, &cfglim_df,
                                  &cfglim_fdf, 2, bike};

  // Maximum lean initial guess
  x->data[0] = M_PI / 3.0;    // Lean
  x->data[1] = M_PI / 4.0;    // Pitch
  lean = lean_max;            // set lean to point at max lean vector
  pitch = pitch_max;          // set pitch to point at max pitch vector
  for (int c = 0; c < 2; x->data[0] = - M_PI / 3.0, // min lean i.g.
                         x->data[1] = M_PI / 4.0,
                         lean = lean_min,     // point at min lean vector
                         pitch = pitch_min,   // point at min pitch vector
                         ++c) {
    gsl_multiroot_fdfsolver_set(s, &f, x);
    for (i = N / 2; i < N - 1; ++i) {
      bike->q3 = steer->data[i];  // steer as a parameter
      iter = 0;
      do
      {
        status = gsl_multiroot_fdfsolver_iterate(s);
        if (status)
          iterateError(status, "cfglim()");
        status = gsl_multiroot_test_residual(s->f, ftol);
      } while(status == GSL_CONTINUE && ++iter < iter_max);
      
      // Increase the tolerance by an order of magnitude to improve convergence
      if (iter == iter_max) {
        increaseftol(&ftol, &i, "cfglim()", bike->q3);
        continue;
      } // if

      // Store the lean into the lean vector
      lean->data[i] = s->x->data[0];
      pitch->data[i] = s->x->data[1];
    } // for i (steer from PI/2 to PI)
    lean->data[i] = lean_max->data[i-1];    // can't solve at steer = 0
    pitch->data[i] = pitch_max->data[i-1];  // can't solve at steer = 0

    x->data[0] = lean->data[N / 2];
    x->data[1] = pitch->data[N / 2];
    gsl_multiroot_fdfsolver_set(s, &f, x);
    for (i = N / 2 - 1; i > 0; --i) {
      bike->q3 = steer->data[i];  // steer as a parameter
      iter = 0;
      do
      {
        status = gsl_multiroot_fdfsolver_iterate(s);
        if (status)
          iterateError(status, "cfglim()");
        status = gsl_multiroot_test_residual(s->f, ftol);
      } while(status == GSL_CONTINUE && ++iter < iter_max);
      
      // Increase the tolerance by an order of magnitude to improve convergence
      if (iter == iter_max) {
        increaseftol(&ftol, &i, "cfglim()", bike->q3);
        continue;
      } // if

      // Store the lean into the lean vector
      lean->data[i] = s->x->data[0];
      pitch->data[i] = s->x->data[1];
    } // for i (steer from PI/2 to O)
    lean->data[0] = lean->data[i+1];    // can't solve at steer = 0
    pitch->data[0] = pitch->data[i+1];  // can't solve at steer = 0
  } // for c

  // Free dynamically allocated variables
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
} // cfglim()

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
  (*ftol) *= 10;
  --(*i);
  cerr << "Max iterations encountered in " << routine << " at steer = " << steer
       << "\nIncreasing tolerance to ftol = " << *ftol << '\n';
} // increaseftol()
