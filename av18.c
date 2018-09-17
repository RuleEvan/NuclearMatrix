
double av18( ) {

double woods_saxon(double r) {
  double a = 0.2; // [fm]
  double r_0 = 0.5 // [fm]
  double w = 1.0 + exp((r - r_0)/a);
  w = 1.0/w;
  
  return w;
}

double t_mu(double r) {
  double x = MU*r;
  double c = 2.1; // [fm]^-2
  double t= (1.0 + 3.0/x + 3.0/pow(x, 2.0));
  t *= exp(-x)/x;
  t *= pow(1.0 - exp(-c*r*r), 2.0);
  
  return t;
}


double v_c_01_pp(double r) {
  double I = -11.27028;
  double P = 3346.6874;
  double Q = 1859.5627;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_01_np(double r) {
  double I = -10.66788;
  double P = 3126.5542;
  double Q = 1746.4298;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_01_nn(double r) {
  double I = -11.27028;
  double P = 3342.7664;
  double Q = 1857.4367;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_l2_01(double r) {
  double I = 0.12472;
  double P = 16.7780;
  double Q = 9.0972;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_00(double r) {
  double I = -2.09971;
  double P = 1204.4301;
  double Q = 511.9380;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_l2_00(double r) {
  double I = -0.31452;
  double P = 217.4559;
  double Q = 117.9063;
  double R = 0;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_11_pp(double r) {
  double I = -7.62701;
  double P = 1815.4920;
  double Q = 969.3863;
  double R = 1847.8059;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_11_np(double r) {
  double I = -7.62701;
  double P = 1813.5315;
  double Q = 966.2483;
  double R = 1847.8059;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_11_nn(double r) {
  double I = -7.62701;
  double P = 1811.5710;
  double Q = 967.2603;
  double R = 1847.8059;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_l2_11(double r) {
  double I = 0.06709;
  double P = 342.0669;
  double Q = 185.4713;
  double R = -615.2339;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_t_11(double r) {
  double I = 1.07985;
  double P = 0.0;
  double Q = -190.0949;
  double R = -811.2040;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_ls_11(double r) {
  double I = -0.62697;
  double P = -570.5571;
  double Q = -309.3605;
  double R = 819.1222;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_ls2_11(double r) {
  double I = 0.74129;
  double P = 9.3418;
  double Q = 5.0652;
  double R = -376.4384;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_c_10(double r) {
  double I = -8.62770;
  double P = 2605.2682;
  double Q = 1459.6345;
  double R = 411.9733;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_l2_10(double r) {
  double I = -0.13201;
  double P = 253.4350;
  double Q = 137.4144;
  double R = -1.0076;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_t_10(double r) {
  double I = 1.485601;
  double P = 0.0;
  double Q = -1126.8359;
  double R = 370.1324;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_ls_10(double r) {
  double I = 0.10180;
  double P = 86.0658;
  double Q = 46.6655;
  double R = -356.5175;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_ls2_10(double r) {
  double I = 0.07357;
  double P = -217.5791;
  double Q = -117.9731;
  double R = 18.3935;

  double x = MU*r;
  double v = I*pow(t_mu(r), 2.0);
  v += woods_saxon(r)*(P + x*Q + x*x*R);

  return v;
}

double v_ci_11(double r) {
  double v = 1.0/3.0*(v_c_11_pp(r) + v_c_11_nn(r) + v_c_11_np(r));
  
  return v;
}

double v_ci_01(double r) {
  double v = 1.0/3.0*(v_c_01_pp(r) + v_c_01_nn(r) + v_c_01_np(r));
  
  return v;
}
 
double v_c(double r) {
  double v = 9.0*v_ci_11(r) + 3.0*v_c_10_np(r) + 3.0*v_ci_01(r) + v_c_00_np(r);
  v *= 1.0/16.0;
  
  return v;
}

double v_tau(double r) {
  double v = 3.0*v_ci_11(r) - 3.0*v_c_10_no(r) + v_ci_01(r) - v_c_00_np(r);
  v *= 1.0/16.0;

  return v;
}

double v_sigma(double r) {
   double v = 3.0*v_ci_11(r) + v_c_10_no(r) -3.0*v_ci_01(r) - v_c_00_np(r);
  v *= 1.0/16.0;

  return v;
}

double v_sigma_tau(double r) {
   double v = v_ci_11(r) - v_c_10_no(r) - v_ci_01(r) + v_c_00_np(r);
  v *= 1.0/16.0;

  return v;
}

double v_t(double r) {
  double v = 3.0*v_t_11(r) + v_t_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_ls(double r) {
  double v = 3.0*v_ls_11(r) + v_ls_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_ls2(double r) {
  double v = 3.0*v_ls2_11(r) + v_ls2_10(r);
  v *= 1.0/4.0;
  
  return v;
}

double v_t_tau(double r) {
  double v = v_t_11(r) - v_t_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_ls_tau(double r) {
double v = v_ls_11(r) - v_ls_10(r);
v *= 1.0/4.0;

return v;
}

double v_ls2_tau(double r) {
  double v = v_ls2_11(r) - v_ls2_10(r);
  v *= 1.0/4.0;

  return v;
}

double v_cd_11(double r) {
  double v = 0.5*(v_c_11_pp(r) + v_c_11_nn(r)) - v_c_11_np(r);
  v *= 1.0/6.0;

  return v;
}

double v_cd_01(double r) {
  double v = 0.5*(v_c_01_pp(r) + v_c_01_nn(r)) - v_c_01_np(r);
  v *= 1.0/6.0;

  return v;
}

double v_cap_t(double r) {
  double v = 3.0*v_cd_11(r) + v_cd_01(r);
  v *= 1.0/4.0;

  return v;
}

double v_sigma_cap_t(double r) {
  double v = v_cd_11(r) - v_cd_01(r);
  v *= 1.0/4.0;

  return v;
}

double v_t_cap_t(double r) {
  double v = 0.5*(v_t_11_pp(r) + v_t_11_nn(r)) - v_t_11_np(r);
  v *= 1.0/6.0;

  return v;
}

double v_ca_11(double r) {
  double v = v_c_11_pp(r) - v_c_11_nn(r);
  v *= 1.0/4.0;

  return v;
}

double v_ca_01(double r) {
  double v = v_c_01_pp(r) - v_c_01_nn(r);
  v *= 1.0/4.0;

  return v;
}

double v_tau_z(double r) {
  double v = 3.0*v_ca_11(r) + v_ca_01(r);
  v *= 1.0/4.0;

  return v;
}

double v_sigma_tau_z(double r) {
  double v = v_ca_11(r) - v_ca_01(r);
  v *= 1.0/4.0;

  return v;
}
