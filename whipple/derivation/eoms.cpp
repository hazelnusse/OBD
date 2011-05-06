#include <iostream>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

void writelatexHeader(void);
ex dot(const matrix & a, const matrix & b);
matrix cross(const matrix &a, const matrix & b);

int main(int argc, char ** argv)
{
  // Symbols for all real quantities in derivation
  realsymbol lr("lr", "l_r"), lf("lf", "l_f"),
             xr("xr", "x_r"), xf("xf", "x_f"),
             zr("zr", "z_r"), zf("zf", "z_f"),
             ls("ls", "l_s"),
             x("x", "x"), xp("xp", "\\dot{x}"), xpp("xpp", "\\ddot{x}"),
             y("y", "y"), yp("yp", "\\dot{y}"), ypp("ypp", "\\ddot{y}"),
             psi("psi", "\\psi"),
             psip("psip", "\\dot{\\psi}"),
             psipp("psipp", "\\ddot{\\psi}"),
             phi("phi", "\\phi"),
             phip("phip", "\\dot{\\phi}"),
             phipp("phipp", "\\ddot{\\phi}"),
             theta("theta", "\\theta"),
             thetap("thetap", "\\dot{\\theta}"),
             thetapp("thetapp", "\\ddot{\\theta}"),
             delta("delta", "\\delta"),
             deltap("deltap", "\\dot{\\delta}"),
             deltapp("deltapp", "\\ddot{\\delta}"),
             thetar("thetar", "\\theta_r"),
             thetarp("thetarp", "\\dot{\\theta}_r"),
             thetarpp("thetarpp", "\\ddot{\\theta}_r"),
             thetaf("thetaf", "\\theta_f"),
             thetafp("thetafp", "\\dot{\\theta}_f"),
             thetafpp("thetafp", "\\ddot{\\theta}_f"),
             wr("wr", "\\omega_r"), wrp("wrp", "\\dot{\\omega}_r"),
             wf("wf", "\\omega_f"), wfp("wfp", "\\dot{\\omega}_f"),
             wx("wx", "\\omega_x"), wxp("wxp", "\\dot{\\omega}_x"),
             wy("wy", "\\omega_y"), wyp("wyp", "\\dot{\\omega}_y"),
             wz("wz", "\\omega_z"), wzp("wzp", "\\dot{\\omega}_z"),
             ws("ws", "\\omega_s"), wsp("wsp", "\\dot{\\omega}_s");

  // Symbols for all positive real quantities
  possymbol rf("rf", "r_f"), rft("rft", "r_{ft}"),
            rr("rr", "r_r"), rrt("rrt", "r_{rt}"),
            Irxx("Irxx", "I_{rxx}"), Iryy("Iryy", "I_{ryy}"),
            Irzz("Irzz", "I_{rzz}"), Irxz("Irxz", "I_{rxz}"),
            Ifxx("Ifxx", "I_{fxx}"), Ifyy("Ifyy", "I_{fyy}"),
            Ifzz("Ifzz", "I_{fzz}"), Ifxz("Ifxz", "I_{fxz}"),
            Irs("Irs", "I_{rs}"), Ifs("Ifs", "I_{fs}"),
            mr("mr", "m_{r}"), mf("mf", "m_{f}");


  // List of generalized speeds
  symbol speeds[6] = {wx, wy, wz, ws, wr, wf};

  // Rotation matrices
  matrix N_A(3,3), A_B(3,3), B_R(3,3), R_F(3,3),
         A_N(3,3), B_A(3,3), R_B, F_R(3,3);
  N_A =    cos(psi),   -sin(psi),           0,
           sin(psi),    cos(psi),           0,
                  0,           0,           1;
  A_B =           1,           0,           0,
                  0,    cos(phi),   -sin(phi),
                  0,    sin(phi),    cos(phi);
  B_R =  cos(theta),           0,  sin(theta),
                  0,           1,           0,
        -sin(theta),           0,  cos(theta);
  R_F =  cos(delta), -sin(delta),           0,
         sin(delta),  cos(delta),           0,
                  0,           0,           1;
  A_N = N_A.transpose();
  B_A = A_B.transpose();
  R_B = B_R.transpose();
  F_R = R_F.transpose();

  // A few basis vectors
  matrix e_x(3, 1, lst(1, 0, 0)),
         e_y(3, 1, lst(0, 1, 0)),
         e_z(3, 1, lst(0, 0, 1));

  // Gyrostat inertia matrices
  matrix * IRCG = new matrix(3, 3), * IFCG = new matrix(3, 3);
  *IRCG = Irxx,    0, Irxz,
             0, Iryy,    0,
          Irxz,    0, Irzz;
  *IFCG = Ifxx,    0, Ifxz,
             0, Ifyy,    0,
          Ifxz,    0, Ifzz;

  // Front Gyrostat Calculations
  // All front gyrostat vectors and matrix quantities are expressed in the fork
  // coordinate system
  //
  // Position of front wheel center relative to front wheel contact
  //
  // Step 1:  Form unit vector pointing from front wheel center to center of
  // front tire casing directly above front wheel contact
  // gz = unitvec(n3> - dot(n3>, f2>)*f2>)
  matrix * gz = new matrix(3, 1);
  *gz = F_R.mul(R_B.mul(B_A.mul(e_z))).sub(e_y.mul_scalar(dot(e_y, F_R.mul(R_B.mul(B_A.mul(e_z))))));
  // TODO:  Autolev generates a more compact representation of the term in the
  // dot(gz, gz)... verify that they are equivalent
  *gz = (*gz).mul_scalar(1.0/sqrt(1.0 - pow(sin(phi)*cos(delta) + sin(delta)*sin(theta)*cos(phi), 2)));
  //*gz = (*gz).mul_scalar(1.0/sqrt(dot(*gz, *gz)));

  // Step 2:  Form position from front wheel contact to front wheel center
  matrix * r_fn_fbo = new matrix(3, 1);
  *r_fn_fbo = ((*gz).mul_scalar(-rf)).add(F_R.mul(R_B.mul(B_A.mul(e_z))).mul_scalar(-rft));

  (*r_fn_fbo)(0, 0) =  (rf/(sqrt(1.0 - pow(sin(phi)*cos(delta) + sin(delta)*sin(theta)*cos(phi), 2))) + rft)*(sin(theta)*cos(delta)*cos(phi)-sin(phi)*sin(delta));
  (*r_fn_fbo)(2, 0) = -(rf/(sqrt(1.0 - pow(sin(phi)*cos(delta) + sin(delta)*sin(theta)*cos(phi), 2))) + rft)*cos(theta)*cos(phi);

  // Position of front gyrostat mass center relative to front wheel center,
  matrix * r_fbo_fgo = new matrix(3, 1, lst(xf, 0, zf));

  // Front gyrostat carrier angular velocity
  matrix * w_fa_n = new matrix(3, 1, lst(wx, wy, wz));

  // Front gyrostat carrier partial angular velocities
  // the r-th column represents the r-th partial angular velocity
  matrix * p_w_fa_n = new matrix(3, 6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*p_w_fa_n)(i, j) = diff((*w_fa_n)(i, 0), speeds[j]);
    }
  }

  // Front gyrostat rotor angular velocity
  matrix * w_fb_n = new matrix(3, 1, lst(wx, wf, wz));

  // First vector term in Eqn 23 of Mitiguy and Reckdahl
  matrix * w_fa_n_x_wf_fy = new matrix(3, 1);
  * w_fa_n_x_wf_fy = cross(*w_fa_n, matrix(3, 1, lst(0, wf, 0)));

  // Second vector term in Eqn 23 of Mitiguy and Reckdahl
  matrix * alf_f_spin = new matrix(3, 1, lst(0, wfp, 0));

  // Front gyrostat rotor partial angular velocities
  // the r-th column represents the r-th partial angular velocity
  matrix * p_w_fb_n = new matrix(3, 6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*p_w_fb_n)(i, j) = diff((*w_fb_n)(i, 0), speeds[j]);
    }
  }

  // Front gyrostat mass center velocity
  matrix * v_fgo_n = new matrix(3, 1);
  *v_fgo_n = (cross(*w_fb_n, *r_fn_fbo)).add(cross(*w_fa_n, *r_fbo_fgo));
  (*v_fgo_n)(0,0) = collect((*v_fgo_n)(0,0), lst(wx, wy, wz, wf));
  (*v_fgo_n)(1,0) = collect((*v_fgo_n)(1,0), lst(wx, wy, wz, wf));
  (*v_fgo_n)(2,0) = collect((*v_fgo_n)(2,0), lst(wx, wy, wz, wf));

  // Front gyrostat mass center partial velocities
  // the r-th column represents the r-th partial velocity
  matrix * p_v_fgo_n = new matrix(3, 6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*p_v_fgo_n)(i, j) = diff((*v_fgo_n)(i, 0), speeds[j]);
    }
  }

  // Front gyrostat carrier angular acceleration
  matrix * alf_fa_n = new matrix(3, 1, lst(wxp, wyp, wzp));

  // Front gyrostat rotor angular acceleration
  matrix * alf_fb_n = new matrix(3, 1);
  for (int i = 0; i < 3; ++i) {
    (*alf_fb_n)(i, 0) = diff((*w_fb_n)(i, 0), wx)*wxp
                      + diff((*w_fb_n)(i, 0), wf)*wfp
                      + diff((*w_fb_n)(i, 0), wz)*wzp;
  }
  *alf_fb_n = (*alf_fb_n).add(cross(*w_fa_n, *w_fb_n));

  // Front gyrostat mass center acceleration
  matrix * a_fgo_n = new matrix(3, 1);
  for (int i = 0; i < 3; ++i) {
    (*a_fgo_n)[i] = diff((*v_fgo_n)(i, 0), wx)*wxp
                  + diff((*v_fgo_n)(i, 0), wy)*wyp
                  + diff((*v_fgo_n)(i, 0), wz)*wzp
                  + diff((*v_fgo_n)(i, 0), ws)*wsp
                  + diff((*v_fgo_n)(i, 0), wr)*wrp
                  + diff((*v_fgo_n)(i, 0), wf)*wfp
                  + diff((*v_fgo_n)(i, 0), phi)*phip
                  + diff((*v_fgo_n)(i, 0), theta)*thetap
                  + diff((*v_fgo_n)(i, 0), delta)*deltap;
  }
  *a_fgo_n = (*a_fgo_n).add(cross(*w_fa_n, *v_fgo_n));

  // Front gyrostat steady turning mass center acceleration
  matrix * a_fgo_n_steady = new matrix(3, 1);
  *a_fgo_n_steady = cross(*w_fa_n, *v_fgo_n);

  // TFCG
  matrix * TFCG = new matrix(3,1);
  *TFCG = (*IFCG).mul(*alf_fa_n).add(cross(*w_fa_n, (*IFCG).mul(*w_fa_n)));

  // TFCG_steady
  matrix * TFCG_steady = new matrix(3,1);
  *TFCG_steady = cross(*w_fa_n, (*IFCG).mul(*w_fa_n));

  // Rear Gyrostat Calculations
  // All front gyrostat vectors and matrix quantities are expressed in the R
  // coordinate system
  //
  // Position of rear wheel center relative to rear wheel contact
  matrix * r_rn_rbo = new matrix(3, 1);
  *r_rn_rbo = (R_B.mul(e_z)).mul_scalar(-rr).add(R_B.mul(B_A.mul(e_z)).mul_scalar(-rrt));

  // Position of rear gyrostat mass center relative to rear wheel center
  matrix * r_rbo_rgo = new matrix(3, 1, lst(xr, 0, zr));

  // Rear gyrostat carrier angular velocity
  matrix * w_ra_n = new matrix(3, 1);
  *w_ra_n = R_F.mul((*w_fa_n).sub(e_z.mul_scalar(dot(*w_fa_n, e_z))).add(e_z.mul_scalar(ws)));

  // Rear gyrostat carrier partial angular velocities
  // the r-th column represents the r-th partial angular velocity
  matrix * p_w_ra_n = new matrix(3, 6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*p_w_ra_n)(i, j) = diff((*w_ra_n)(i, 0), speeds[j]);
    }
  }

  // Rear gyrostat rotor angular velocity
  matrix * w_rb_n = new matrix(3, 1,
                               lst((*w_ra_n)[0], wr, (*w_ra_n)[2]));

  // Rear gyrostat rotor partial angular velocities
  // the r-th column represents the r-th partial angular velocity
  matrix * p_w_rb_n = new matrix(3, 6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*p_w_rb_n)(i, j) = diff((*w_rb_n)(i, 0), speeds[j]);
    }
  }

  // First vector term in Eqn 23 of Mitiguy and Reckdahl
  matrix * w_ra_n_x_wr_ry = new matrix(3, 1);
  * w_ra_n_x_wr_ry = cross(*w_ra_n, matrix(3, 1, lst(0, wr, 0)));

  // Second vector term in Eqn 23 of Mitiguy and Reckdahl
  matrix * alf_r_spin = new matrix(3, 1, lst(0, wrp, 0));

  // Rear gyrostat angular accleration
  matrix * alf_ra_n = new matrix(3, 1);
  for (int i = 0; i < 3; ++i) {
    (*alf_ra_n)[i] = diff((*w_ra_n)[i], wx)*wxp
                   + diff((*w_ra_n)[i], wy)*wyp
                   + diff((*w_ra_n)[i], ws)*wsp
                   + diff((*w_ra_n)[i], delta)*deltap;
  }

  // Rear gyrostat rotor angular acceleration
  matrix * alf_rb_n = new matrix(3, 1);
  for (int i = 0; i < 3; ++i) {
    (*alf_rb_n)[i] = diff((*w_rb_n)[i], wx)*wxp
                   + diff((*w_rb_n)[i], wy)*wyp
                   + diff((*w_rb_n)[i], ws)*wsp
                   + diff((*w_rb_n)[i], wr)*wrp
                   + diff((*w_rb_n)[i], delta)*deltap;
  }
  *alf_rb_n = (*alf_rb_n).add(cross(*w_ra_n, *w_rb_n));

  // Rear gyrostat mass center velocity
  matrix * v_rgo_n = new matrix(3, 1);
  *v_rgo_n = cross(*w_rb_n, *r_rn_rbo).add(cross(*w_rb_n, *r_rbo_rgo));
  (*v_rgo_n)[0] = collect((*v_rgo_n)[0], lst(wx, wy, wz, ws, wr, wf));
  (*v_rgo_n)[1] = collect((*v_rgo_n)[1], lst(wx, wy, wz, ws, wr, wf));
  (*v_rgo_n)[2] = collect((*v_rgo_n)[2], lst(wx, wy, wz, ws, wr, wf));

  // Front gyrostat mass center partial velocities
  // the r-th column represents the r-th partial velocity
  matrix * p_v_rgo_n = new matrix(3, 6);
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*p_v_rgo_n)(i, j) = diff((*v_rgo_n)(i, 0), speeds[j]);
    }
  }

  // Rear gyrostat mass center acceleration
  matrix * a_rgo_n = new matrix(3, 1);
  for (int i = 0; i < 3; ++i) {
    (*a_rgo_n)[i] = diff((*v_rgo_n)(i, 0), wx)*wxp
                  + diff((*v_rgo_n)(i, 0), wy)*wyp
                  + diff((*v_rgo_n)(i, 0), wz)*wzp
                  + diff((*v_rgo_n)(i, 0), ws)*wsp
                  + diff((*v_rgo_n)(i, 0), wr)*wrp
                  + diff((*v_rgo_n)(i, 0), wf)*wfp
                  + diff((*v_rgo_n)(i, 0), phi)*phip
                  + diff((*v_rgo_n)(i, 0), theta)*thetap
                  + diff((*v_rgo_n)(i, 0), delta)*deltap;
  }
  *a_rgo_n = (*a_rgo_n).add(cross(*w_ra_n, *v_rgo_n));

  // Rear gyrostat steady mass center acceleration
  matrix * a_rgo_n_steady = new matrix(3, 1);
  *a_rgo_n_steady = cross(*w_ra_n, *v_rgo_n);

  // TRCG
  matrix * TRCG = new matrix(3,1);
  *TRCG = (*IRCG).mul(*alf_ra_n).add(cross(*w_ra_n, (*IRCG).mul(*w_ra_n)));

  // TRCG_steady
  matrix * TRCG_steady = new matrix(3,1);
  *TRCG_steady = cross(*w_ra_n, (*IRCG).mul(*w_ra_n));

  // Holonomic constraint
  // -rrt*nz> - rf*bz> + lr*rx> + ls*rz> + lf*fx> + rf*gz> + rft*nz>
  ex hc;
  hc = (rft - rrt)
     + collect(dot(A_B.mul(e_z).mul_scalar(-rr), e_z)
     + dot(A_B.mul(B_R.mul(e_x)).mul_scalar(lr), e_z)
     + dot(A_B.mul(B_R.mul(e_z)).mul_scalar(ls), e_z), cos(phi))
     + dot(((A_B.mul(B_R.mul(R_F))).mul(e_x)).mul_scalar(lf), e_z)
     + rf*sqrt(1.0 - pow(sin(phi)*cos(delta) + sin(delta)*sin(theta)*cos(phi), 2))
     ;
     //  add(A_B.mul(B_R.mul(e_x.mul_scalar(lr).add(e_z.mul_scalar(ls))))).
     //  add((A_B.mul(B_R.mul(R_F.mul((*gz))))).mul_scalar(rf)))(2, 0),
     //  lst(rft, rrt, rr, rf, lr, ls, lf));

  ex Frstar_R[6], Frstar_F[6], *tmpr, *tmpf;
  for (int r = 0; r < 6; ++r) {
    Frstar_R[r] = 0;
    Frstar_F[r] = 0;

    /*
    // First term in Eqn 22 of Mitiguy and Reckdahl
    // This one is the nastiest expression of all, it comes from dotting the
    // partial velocities with the accelerations
    tmpr = new ex;
    *tmpr = 0;
    tmpf = new ex;
    *tmpf = 0;
    for (int i = 0; i < 3; ++i) {
      *tmpr += (*p_v_rgo_n)(i, r)*((*a_rgo_n)(i, 0));
      *tmpf += (*p_v_fgo_n)(i, r)*((*a_fgo_n)(i, 0));
    }
    Frstar_R[r] -= mr*(*tmpr);  // multiply the by rear assembly mass and subtract
    Frstar_F[r] -= mf*(*tmpf);  // multiply the by front assembly mass and subtract
    delete tmpr;
    delete tmpf;
    */

    // Second term in Eqn 22 of Mitiguy and Reckdahl
    tmpr = new ex;
    *tmpr = 0;
    tmpf = new ex;
    *tmpf = 0;
    for (int i = 0; i < 3; ++i) {
      *tmpr += (*p_w_ra_n)(i, r)*((*TRCG)(i, 0));
      *tmpf += (*p_w_fa_n)(i, r)*((*TFCG)(i, 0));
    }
    Frstar_R[r] -= *tmpr;  // Minus so that we end up with Eqn 24
    Frstar_F[r] -= *tmpf;  // Minus so that we end up with Eqn 24
    delete tmpr;
    delete tmpf;

    // Second term in Eqn 23 of Mitiguy and Reckdahl
    tmpr = new ex;
    *tmpr = 0;
    tmpf = new ex;
    *tmpf = 0;
    for (int i = 0; i < 3; ++i) {
      *tmpr += (*p_w_ra_n)(i, r)*((*w_ra_n_x_wr_ry)(i, 0)); // First term in parenthesis in Eq (23)
      *tmpr += (*p_w_rb_n)(i, r)*((*alf_r_spin)(i,0));
      *tmpf += (*p_w_fa_n)(i, r)*((*w_fa_n_x_wf_fy)(i, 0)); // First term in parenthesis in Eq (23)
      *tmpf += (*p_w_fb_n)(i, r)*((*alf_f_spin)(i,0));
    }
    Frstar_R[r] -= Irs*(*tmpr);  // Minus so that we end up with Eqn 24
    Frstar_F[r] -= Ifs*(*tmpf);  // Minus so that we end up with Eqn 24
    delete tmpr;
    delete tmpf;
  }

  // Output LaTeX
  cout << latex;
  writelatexHeader();
  cout << "Rotation matrices: \n"
       << "\\begin{align}\n  {}^NR^{A} &= " << N_A << "\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^AR^{B} &= " << A_B << "\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^BR^{R} &= " << B_R << "\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^RR^{F} &= " << R_F << "\n" << "\\end{align}\n"
       << "The ordering of the generalized speeds, for purposes of calculating"
          " partial velocity and partial angular velocity matrices, is assumed"
          " to be $["
          << speeds[0] << ", "
          << speeds[1] << ", "
          << speeds[2] << ", "
          << speeds[3] << ", "
          << speeds[4] << ", "
          << speeds[5] << "]$\n\n"
       << "\\begin{align}\n  \\bs{g}_z&= "
       << *gz << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  \\bs{r}^{FBO/FN} &= "
       << *r_fn_fbo << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  \\bs{r}^{FGO/FBO} &= "
       << *r_fbo_fgo << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{FA} &= "
       << *w_fa_n << "_{F}\n" << "\\end{align}\n"
       << "Partial angular velocities of $FA$ in $N$ "
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{FA}_{partials} &= "
       << *p_w_fa_n << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{FB} &= "
       << *w_fb_n << "_{F}\n" << "\\end{align}\n"
       << "Partial angular velocities of $FB$ in $N$ "
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{FB}_{partials} &= "
       << *p_w_fb_n << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{FA} \\times  \\omega_f "
       << e_y << "_{F} &= "
       << *w_fa_n_x_wf_fy << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{v}^{FGO} &= "
       << *v_fgo_n << "_{F}\n" << "\\end{align}\n"
       << "Partial velocities of $FGO$ in $N$ "
       << "\\begin{align}\n  {}^N\\bs{v}^{FGO}_{partials} &= "
       << *p_v_fgo_n << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\alpha}^{FA} &= "
       << *alf_fa_n << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\alpha}^{FB} &= "
       << *alf_fb_n << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{a}^{FGO} &= "
       << *a_fgo_n << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{a}^{FGO}_{steady} &= "
       << *a_fgo_n_steady << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{T}^{FCG} &= "
       << *TFCG << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{T}^{FCG}_{steady} &= "
       << *TFCG_steady << "_{F}\n" << "\\end{align}\n"
       << "\\begin{align}\n  F_x^{F*} &= " << Frstar_F[0] << "\\\\\n"
       << "  F_y^{F*} &= " << Frstar_F[1] << "\\\\\n"
       << "  F_z^{F*} &= " << Frstar_F[2] << "\\\\\n"
       << "  F_s^{F*} &= " << Frstar_F[3] << "\\\\\n"
       << "  F_r^{F*} &= " << Frstar_F[4] << "\\\\\n"
       << "  F_f^{F*} &= " << Frstar_F[5] << "\\end{align}\n"
       << "\\begin{align}\n  \\bs{r}^{RBO/RN} &= "
       << *r_rn_rbo << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  \\bs{r}^{RGO/RBO} &= "
       << *r_rbo_rgo << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{RA} &= "
       << *w_ra_n << "_{R}\n" << "\\end{align}\n"
       << "Partial angular velocities of $RA$ in $N$ "
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{RA}_{partials} &= "
       << *p_w_ra_n << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\alpha}^{RA} &= "
       << *alf_ra_n << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{RB} &= "
       << *w_rb_n << "_{R}\n" << "\\end{align}\n"
       << "Partial angular velocities of $RB$ in $N$ "
       << "\\begin{align}\n  {}^N\\bs{\\omega}^{RB}_{partials} &= "
       << *p_w_rb_n << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{\\alpha}^{RB} &= "
       << *alf_rb_n << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{v}^{RGO} &= "
       << *v_rgo_n << "_{R}\n" << "\\end{align}\n"
       << "Partial velocities of $RGO$ in $N$ "
       << "\\begin{align}\n  {}^N\\bs{v}^{RGO}_{partials} &= "
       << *p_v_rgo_n << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{a}^{RGO} &= "
       << *a_rgo_n << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{a}^{RGO}_{steady} &= "
       << *a_rgo_n_steady << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{T}^{RCG} &= "
       << *TRCG << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  {}^N\\bs{T}^{RCG}_{steady} &= "
       << *TRCG_steady << "_{R}\n" << "\\end{align}\n"
       << "\\begin{align}\n  F_x^{R*} &= " << Frstar_R[0] << "\\\\\n"
       << "  F_y^{R*} &= " << Frstar_R[1] << "\\\\\n"
       << "  F_z^{R*} &= " << Frstar_R[2] << "\\\\\n"
       << "  F_s^{R*} &= " << Frstar_R[3] << "\\\\\n"
       << "  F_r^{R*} &= " << Frstar_R[4] << "\\\\\n"
       << "  F_f^{R*} &= " << Frstar_R[5] << "\\end{align}\n"
       << "Holonomic constraint\n"
       << "\\begin{align}\n  f_{hc}(\\phi, \\theta, \\delta) &= "
       << hc << "\n" << "\\end{align}\n"
       << "\\end{document}\n";


  // Free dynamically allocated memory
  // Front gyrostat quantities
  delete w_fa_n;
  delete p_w_fa_n;
  delete w_fb_n;
  delete w_fa_n_x_wf_fy;
  delete alf_f_spin;
  delete p_w_fb_n;
  delete alf_fa_n;
  delete alf_fb_n;
  delete r_fbo_fgo;
  delete r_fn_fbo;
  delete gz;
  delete v_fgo_n;
  delete p_v_fgo_n;
  delete a_fgo_n;
  delete a_fgo_n_steady;
  delete TFCG;
  delete TFCG_steady;

  // Rear gyrostat quantities
  delete r_rn_rbo;
  delete r_rbo_rgo;
  delete w_ra_n;
  delete p_w_ra_n;
  delete w_rb_n;
  delete w_ra_n_x_wr_ry;
  delete alf_r_spin;
  delete p_w_rb_n;
  delete alf_ra_n;
  delete alf_rb_n;
  delete v_rgo_n;
  delete p_v_rgo_n;
  delete a_rgo_n;
  delete a_rgo_n_steady;
  delete TRCG;
  delete TRCG_steady;

  // Leave main
  return 0;
}
/*
  // Steady turning conditions
  exmap m;
  m[phip] = 0;
  m[thetap] = 0;
  m[deltap] = 0;
  m[wxp] = 0;
  m[wyp] = 0;
  m[wzp] = 0;
  m[wsp] = 0;
  m[wfp] = 0;
  m[wrp] = 0;


  ex v_gfo[3][6];
  symbol speeds[6] = {wx, wy, wz, ws, wr, wf};
  cout << "[\n";
  for (int i = 0; i < 3; ++i) {
    cout << "[";
    for (int j = 0; j < 6; ++j)
      cout << (v_gfo[i][j] = diff(v_gfo_n[i], speeds[j])) << ", ";
    cout << "]\n";
  }
  cout << "]\n";
  ex v_gfo_grad[3][3][6];
  symbol coords[3] = {phi, theta, delta};
  for (int i = 0; i < 3; ++i) {
    cout << "[\n";
    for (int j = 0; j < 3; ++ j) {
      cout << "[";
      for (int k = 0; k < 6; ++k) {
        cout << (v_gfo_grad[i][j][k] = diff(v_gfo[j][k], coords[j])) << ", ";
      }
      cout << "]\n";
    }
    cout << "]\n\n";
  }

*/

void writelatexHeader(void)
{
  cout <<
"\\documentclass[landscape,a0paper,11pt]{article}\n"
"\\usepackage[margin=1in,centering]{geometry}\n"
"\\usepackage{amsmath}\n"
"\\usepackage{amssymb}\n"
"\\renewcommand{\\b}[1]{ \\mathbf{ #1 } }\n"
"\\newcommand{\\bs}[1]{ \\boldsymbol{ #1 } }\n"
"\\begin{document}\n\n";
}

ex dot(const matrix & a, const matrix & b)
{
  return a(0, 0)*b(0, 0) + a(1, 0)*b(1, 0) + a(2, 0)*b(2, 0);
}

matrix cross(const matrix &a, const matrix & b)
{
  return matrix(3, 1, lst(a(1,0)*b(2,0) - a(2,0)*b(1,0),
                          a(2,0)*b(0,0) - a(0,0)*b(2,0),
                          a(0,0)*b(1,0) - a(1,0)*b(0,0)));
}

