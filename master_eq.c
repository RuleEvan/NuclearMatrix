#include "master_eq.h"

void master_equation() {
  double g_1_pp = 0.36;
  double g_2_pp = 2.0;
  double g_3_pp = -0.62;
  double g_4_pp = -1.9;
  double g_5_pp = -8.0;
  double g_t_pp = 1.0;
  double gtp = 1.0;
  double gtnn = 1.0;
  double gtpn = 1.0;
  double gnunn = 1.0;
  double gtpp = 1.0;

  double c_sl_6 = 1.0;
  double c_sr_6 = 1.0;
  double c_v_l_6 = 1.0;
  double c_v_r_6 = 1.0;
  double c_t_6 = 1.0;
  double c_pp_l_9 = 1.0;
  double c_pp_r_9 = 1.0;
  double c_pn_l_9 = 1.0;
  double c_pn_r_9 = 1.0;
  double c_nn_l_9 = 1.0;
  double c_nn_r_9 = 1.0;

  double g_ev_l = 1.0;
  double g_v_me_l = 1.0;
  double g_v_me_r = 1.0;
  double g_v_pn_l = 1.0;
  double g_v_pn_tilde = 1.0;
  double c_v_9 = 1.0;
  double c_v_tilde_9 = 1.0;
  double g_6_nn = 1.0;
  double g_7_nn = 1.0;
  double g_v_nn_l = 1.0;
  double g_v_pm_l = 1.0;
  double m_betabeta = 1.0;

  double c_vl_7 =1.0;
  double c_vr_7 = 1.0;
  double c_v_r_r = 1.0;
  double g_v_pn = 1.0;

  double m_f = compute_matrix_element_tau_plus(-1);
  double m_f_sd = compute_matrix_element_tau_plus(-2);
  double m_gt_aa = compute_matrix_element_sigma_tau_plus(-1);
  double m_gt_aa_sd = compute_matrix_element_sigma_tau_plus(-2);
  double m_t_aa = compute_matrix_element_TT(-1);
  double m_t_aa_sd = compute_matrix_element_TT(-2);
  double m_gt_ap = compute_matrix_element_sigma_tau_plus(-3);
  double m_gt_ap_sd = compute_matrix_element_sigma_tau_plus(-4);
  double m_t_ap = compute_matrix_element_TT(-3);
  double m_t_ap_sd = compute_matrix_element_TT(-4);
  double m_gt_pp = compute_matrix_element_sigma_tau_plus(-5);
  double m_gt_pp_sd = compute_matrix_element_sigma_tau_plus(-6);
  double m_t_pp = compute_matrix_element_TT(-5);
  double m_t_pp_sd = compute_matrix_element_TT(-6);
  double m_gt_mm = compute_matrix_element_sigma_tau_plus(-7);
  double m_gt_mm_sd = compute_matrix_element_sigma_tau_plus(-8);
  double m_t_mm = compute_matrix_element_TT(-7);
  double m_t_ap_mm = compute_matrix_element_TT(-8);

  double m_gt = m_gt_aa + m_gt_ap + m_gt_pp + m_gt_mm;
  double m_t = m_t_ap + m_t_pp + m_t_mm;
  double m_ps = 0.5*m_gt_ap + m_gt_pp + 0.5*m_t_ap + m_t_pp;
  double m_t6 = 2.0*(gtp - gtnn)/pow(G_AXIAL, 2.0)*pow(PION_MASS/M_NEUTRON, 2.0)*m_f_sd - 8.0*G_TENSOR/G_MAGNETIC*(m_gt_mm + m_t_mm);
  m_t6 += gtpn*pow(PION_MASS/(2.0*M_NEUTRON), 2.0)*(m_gt_ap_sd + m_t_ap_sd) + gtpp*pow(PION_MASS/(2.0*M_NEUTRON), 2.0)*(m_gt_pp_sd + m_t_pp_sd);
  
  double m_nu_3 = -V_UD*V_UD*(-1.0/pow(G_AXIAL, 2.0)*m_f + m_gt + m_t + 2.0*pow(PION_MASS/G_AXIAL, 2.0)*gnunn*m_f_sd);
  double m_nu_6 = V_UD*((2700.0/M_NEUTRON)*(c_sl_6 - c_sr_6)+pow(PION_MASS, 2.0)/(M_NEUTRON*HIGGS_VEV*1000.0)*(c_vl_7 - c_vr_7)*m_ps) + V_UD*c_t_6*m_t6;
  double m_nu_9 = -1.0/(2.0*pow(M_NEUTRON, 2.0))*c_pp_l_9*(0.5*m_gt_ap_sd + m_gt_pp_sd + 0.5*m_t_ap_sd + m_t_pp_sd);
  m_nu_9 += pow(PION_MASS/M_NEUTRON, 2.0)/2.0*c_pn_l_9*(m_gt_ap_sd + m_t_ap_sd) - 2.0/pow(G_AXIAL, 2.0)*pow(PION_MASS/M_NEUTRON, 2.0)*c_nn_l_9*m_f_sd;
  double m_r_9 = -1.0/(2.0*pow(M_NEUTRON, 2.0))*c_pp_r_9*(0.5*m_gt_ap_sd + m_gt_pp_sd + 0.5*m_t_ap_sd + m_t_pp_sd);
  m_r_9 += pow(PION_MASS/M_NEUTRON, 2.0)/2.0*c_pn_r_9*(m_gt_ap_sd + m_t_ap_sd) - 2.0/pow(G_AXIAL, 2.0)*pow(PION_MASS/M_NEUTRON, 2.0)*c_nn_r_9*m_f_sd;

  double m_e_l_6 = -V_UD*c_v_l_6/3.0*(pow(G_VECTOR/G_AXIAL, 2.0)*m_f + 1.0/3.0*(2.0*m_gt_aa + m_t_aa) + 6.0*g_ev_l/pow(G_AXIAL, 2.0)*m_f_sd);
  double m_e_r_6 = -V_UD*c_v_r_r/3.0*(pow(G_VECTOR/G_AXIAL, 2.0)*m_f - 1.0/3.0*(2.0*m_gt_aa + m_t_aa) + 6.0*g_ev_l/pow(G_AXIAL, 2.0)*m_f_sd);

  double m_me_l_6 = V_UD*c_v_l_6/6.0*(pow(G_VECTOR/G_AXIAL, 2.0)*m_f - 1.0/3.0*(m_gt_aa - 4.0*m_t_aa) - 3.0*(m_gt_ap + m_gt_pp + m_t_ap + m_t_pp) - 12.0*g_v_me_l/pow(G_AXIAL, 2.0)*m_f_sd);
  double m_me_r_6 = V_UD*c_v_r_6/6.0*(pow(G_VECTOR/G_AXIAL, 2.0)*m_f + 1.0/3.0*(m_gt_aa - 4.0*m_t_aa) + 3.0*(m_gt_ap + m_gt_pp + m_t_ap + m_t_pp) - 12.0*g_v_me_r/pow(G_AXIAL, 2.0)*m_f_sd);

  double m_m_6 = V_UD*c_v_l_6*(2.0*G_AXIAL/G_MAGNETIC*(m_gt_mm + m_t_mm) + pow(PION_MASS/M_NEUTRON, 2.0)*(-2.0/pow(G_AXIAL, 2.0)*g_v_nn_l*m_f_sd + 0.5*g_v_pn_l*(m_gt_ap_sd + m_t_ap_sd)));
  double m_m_9 = pow(PION_MASS/M_NEUTRON, 2.0)*(-2.0/pow(G_AXIAL, 2.0)*(g_6_nn*c_v_9 + g_7_nn*c_v_tilde_9)*m_f_sd + 0.5*(g_v_pn*c_v_9 + g_v_pn_tilde*c_v_tilde_9)*(m_gt_ap_sd + m_t_ap_sd));

  double a_nu = m_betabeta/M_ELECTRON*m_nu_3 + M_NEUTRON/M_ELECTRON*m_nu_6 + pow(M_NEUTRON, 2.0)/(M_ELECTRON*HIGGS_VEV*1000.0)*m_nu_9;
  double a_r = pow(M_NEUTRON, 2.0)/(M_ELECTRON*HIGGS_VEV*1000.0)*m_r_9;
  double a_e = m_e_l_6 + m_e_r_6;
  double a_me = m_me_l_6 + m_me_r_6;
  double a_m = M_NEUTRON/M_ELECTRON*m_m_6 + pow(M_NEUTRON, 2.0)/(M_ELECTRON*HIGGS_VEV*1000.0);

  double t_half = G_0(1)*(a_nu*a_nu + a_r*a_r) - 2.0*(G_0(1) - G_0(4))*a_nu*a_r + 4.0*G_0(2)*a_e*a_e;
  t_half += 2.0*G_0(4)*(a_me*a_me + a_me*(a_nu + a_r)) - 2.0*G_0(3)*((a_nu + a_r)*a_e + 2.0*a_me*a_e);
  t_half += G_0(9)*a_m*a_m + G_0(6)*(a_nu - a_r)*a_m;
  t_half *= pow(G_AXIAL, 4.0);
  printf("Half_life: %g\n", 1.0/t_half);

  return;
}  
