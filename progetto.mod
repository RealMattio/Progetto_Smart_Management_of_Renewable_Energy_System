# Sets
set J ordered;
set I ordered;

#Parameters
param Pnom_i;
param Pnom_hvac;
param alpha;
param beta;
param eta_c;
param eta_h;
param R;
param Pnom_b;
param Eb;
param eta_ch;
param eta_dsc;
param Tsp;
param Delta;
param SoCmin;
param SoCmax;
param Pnom_pev;
param Epev;
param eta_pev;
param Eabp{I};
param Tabp{I};
param Tabp_done{I};
param Eabp_done{I};
param d_abp_k1{I};
param s_abp_k1{I};
param t_abp_done{I};
param Dabp{I};
param Pmax_abp{I};
param Pmin_abp{I};
param Tk;
param SoCk;
param SoCpevk;
param SoCpev_obj;
param P_PV_forecast{J};
param Tex_forecast{J};
param c{J};
param tdSCPV;
param Dts;
param UR_hvac{J};
param UR_pev{J};
param UR_abp{J};

# Variables
var Pc{J} >= 0;
var Ph{J} >= 0;
var Pch{J} >= 0;
var Pdsc{J} >= 0;
var Ppev{J} >= 0;
var Pabp{J,I} >= 0;
var Pabp_tot{J} >= 0;
var Pcurt{J}>=0;
var Pi{J}>=0;
var T{J};
var SoC{J};
var eps{J}>=0;
var d_c{J} binary;
var d_h{J} binary;
var d_b{J} binary;
var d_abp{J,I} binary;
var s_abp{J,I} binary;
var t_abp{J,I} binary;


# Objective function
minimize total_cost : sum {j in J} (c[j]*Pi[j] + tdSCPV*Pcurt[j]+100*eps[j]);

#Constraints for HVAC
subject to con_hvac_1 {j in J: ord(j)=1}:       T[j] == alpha*Tk-beta*R*(eta_c*Pc[j]-eta_h*Ph[j])+beta*Tex_forecast[j];
subject to con_hvac_2 {j in J: ord(j)>1}:       T[j] == alpha*T[j-1]-beta*R*(eta_c*Pc[j]-eta_h*Ph[j])+beta*Tex_forecast[j];
subject to con_hvac_3 {j in J}:                 T[j]*UR_hvac[j] <= (Tsp+Delta)*UR_hvac[j]+eps[j];
subject to con_hvac_4 {j in J}:                 T[j]*UR_hvac[j] >= (Tsp-Delta)*UR_hvac[j]-eps[j];
subject to con_hvac_5 {j in J}:                 Pc[j]<= d_c[j]*Pnom_hvac;
subject to con_hvac_6 {j in J}:                 Ph[j]<= d_h[j]*Pnom_hvac;
subject to con_hvac_7 {j in J}:                 d_h[j] + d_c[j]<= UR_hvac[j];

#Constraints for battery
subject to con_battery_1  {j in J: ord(j)=1}:   SoC[j] == SoCk + Dts/Eb*(eta_ch*Pch[j]-1/eta_dsc*Pdsc[j]);
subject to con_battery_2  {j in J: ord(j)>1}:   SoC[j] == SoC[j-1] + Dts/Eb*(eta_ch*Pch[j]-1/eta_dsc*Pdsc[j]);
subject to con_battery_3  {j in J}:             SoC[j] <= SoCmax;
subject to con_battery_4  {j in J}:             SoC[j] >= SoCmin;
subject to con_battery_5  {j in J}:             Pch[j]<= d_b[j]*Pnom_b;
subject to con_battery_6  {j in J}:             Pdsc[j]<= (1-d_b[j])*Pnom_b;

#Contraints for PEV
subject to con_pev_1:                           sum {j in J} (Dts*eta_pev*Ppev[j]) == Epev*(SoCpev_obj-SoCpevk);
subject to con_pev_2 {j in J}:                  Ppev[j]<= UR_pev[j]*Pnom_pev;

#Constraints for ABP
subject to con_abp_1 {i in I}:                            sum {j in J} (Dts*Pabp[j,i]) = Eabp[i]-Eabp_done[i]; 
subject to con_abp_2 {i in I}:                            sum {j in J} (d_abp[j,i]) = Tabp[i]-Tabp_done[i];
subject to con_abp_3 {j in J, i in I}:                    Pabp[j,i]<= d_abp[j,i]*Pmax_abp[i];
subject to con_abp_4 {j in J, i in I}:                    Pabp[j,i]>= d_abp[j,i]*Pmin_abp[i];
subject to con_abp_5 {j in J, i in I}:                    d_abp[j,i]+s_abp[j,i]<=1;
subject to con_abp_6 {j in J, i in I: ord(j)=1}:          d_abp_k1[i]-d_abp[j,i]-s_abp[j,i]<=0;
subject to con_abp_7 {j in J, i in I: ord(j)>1}:          d_abp[j-1,i]-d_abp[j,i]-s_abp[j,i]<=0;
subject to con_abp_8 {j in J, i in I: ord(j)=1}:          s_abp_k1[i]-s_abp[j,i]<=0;
subject to con_abp_9 {j in J, i in I: ord(j)>1}:          s_abp[j-1,i]-s_abp[j,i]<=0;
subject to con_abp_10 {j in J, i in I: ord(i)>1}:          d_abp[j,i]-s_abp[j,i-1]<=0;
subject to con_abp_11 {j in J, i in I: ord(i)>1}:         t_abp[j,i] == s_abp[j,i-1]-d_abp[j,i]-s_abp[j,i];
subject to con_abp_12 {i in I: ord(i)>1}:                 sum {j in J} (t_abp[j,i]) <= Dabp[i]-t_abp_done[i]; 
subject to con_abp_13 {j in J, i in I}:                   d_abp[j,i] <= UR_abp[j];
subject to con_abp_14 {j in J}:                           Pabp_tot[j] == sum {i in I} (Pabp[j,i]); 

#Constraints for system
subject to con_system_1 {j in J}:                 Pi[j] == Pc[j]+Ph[j]+Pch[j]-Pdsc[j]-P_PV_forecast[j]+Pcurt[j]+Ppev[j]+Pabp_tot[j];
subject to con_system_2 {j in J}:                 Pcurt[j] <= P_PV_forecast[j];
subject to con_system_3 {j in J}:                 Pi[j] <= Pnom_i;