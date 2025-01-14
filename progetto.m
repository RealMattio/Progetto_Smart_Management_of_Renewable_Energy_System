%% Problem data

% Network import 
parms.Pnom_i = 250; % maximal imported power [kW]
parms.Pnom_e = 120; % maximal exported power [kW]

% Non-Renewable DG
parms.Pnom_g = 80; % maximal power DG [kW]
parms.eta_g = 0.6; % DGs efficiency
parms.cf = 0.45; % price of fuel [€/kWh]

% HVAC data
parms.R = 1.9; %walls thermal resistance [°C/kW]
parms.C = 3.925; % thermal capacitance of air [kWh/°C]
parms.eta_c = 1.9; % cooling system COP
parms.eta_h = 0.85; % hating system efficiency
parms.Pnom_hvac = 12; % nominal power of the cooling system [kW]

% PV generator data
PV_panel_Area = 1.6275; % Single panel area [m2]
npanels = 332; % number of panels
PV_Area =  npanels*PV_panel_Area; % Total PV area
eta_pv = 0.222; % PV system efficiency
parms.Pnom_PV = 120; % Nominal power at 1kW/m2 [kW]

% Battery data
parms.Pnom_b = 150; % Battery nominal power [kW]
parms.Eb = 150; % Battery capaciy [kWh]
parms.eta_ch = 0.95; % Battery charging efficiency [pu]
parms.eta_dsc = 0.95; % Battery discharging efficiency [pu]

% ABP data
parms.n_abp_phases = 4; % number of phases [-]
parms.Eabp = [11 22 31 14]'; % phases energy [kWh]
parms.Tabp = [2 2 2 2]'; % phase time execution [time steps]
parms.Dabp = [0 2 2 2]'; % maximal delay among phases [time steps]
parms.Pmax_abp = [11 22 31 14]'; % maximal power for each phase [kW]
parms.Pmin_abp = [0 0 0 0]'; % minimal power for each phase [kW]

% Selection of days
months_days = [31 28 31 30 31 30 31 31 30 31 30 31];
m1 = 6; %starting month
d1 = 15; %staring day
m2 = 6; %end month
d2 = 18; %end day
idxs = (sum(months_days(1:m1-1))+(d1-1))*24+1:(sum(months_days(1:m2-1))+d2)*24+1;

% Control specifications data
% ATTENTION: non è implementato il controllo sull'impostazione del setpoint
% della temperatura! Pertanto in questa fase è necessario che il periodo di
% simulazione corrisponda allo stesso mese
parms.Tsp_winter = 20; % temperature set-point from April to September [°C]
parms.Tsp_summer = 22; % temperature set-point from October to March
if m1 >= 4 & m1 <= 9 % if the starting month is a summer month Tsp will be setuped as Tsp summer
    parms.Tsp = parms.Tsp_summer;
else
    parms.TSP = parms.Tsp_winter;
end
parms.Delta  = 2; % temperature regulation tolerance [°C] 
parms.SoCmin = 0.1; % desired minimal battery SoC [pu]
parms.SoCmax = 0.9; % desired maximal battery SoC [pu]
parms.Dts = 1; % control sampling time [h]


% External temperature data 
load T_ex_rome_campus_bio_medico_2022.mat
% [month,day,hour,forecasted temperature (°C), actual temperature (°C)]
T_ex = T_ex(idxs,:);

% Solar irradiation data
load FALSA_previsione_irraggiamento.mat
% [day,forecasted Ir (°C), actual Ir (°C)]
Ir = Ir(idxs,:);

% VANNO AGGIUNTI GLI UFFICI!!

% Prices data 
load energy_prices_january_july_2021

priceF1 = 0.53276; % F1 price [eur/kWh]
priceF2 = 0.54858; % F2 price [eur/kWh]
priceF3 = 0.46868; % F3 price [eur/kWh]
% assuming working days
cday = [ones(6,1)*priceF3;
        priceF2;
        ones(11,1)*priceF1;
        ones(4,1)*priceF2;
        ones(2,1)*priceF3];  
c = [cday;cday;cday;cday]; %two days+1 day for forecasts
%parms.tdSCPV = 8.22/1000; % PV self-consumption tarif discount [€/kWh]


%% Control parameters
parms.T = 18; % control time horizon
battery_compensation = 1; % 1 to activate the error compensation with battery

%% Simulation parameters
Tf = 24*3;
T0 = parms.Tsp-1; % Initial temperature
SoC0 = 0.5; % Initial battery state of charge


%% Offline computation
parms.alpha = exp(-parms.Dts/(parms.R*parms.C));
parms.beta = 1 - parms.alpha;

%% Initialization
Pc = zeros(Tf,1);        % air cooling power [kW] 
Ph = zeros(Tf,1);        % air heating power [kW]
Pch = zeros(Tf,1);       % battery charging power [kW]
Pdsc = zeros(Tf,1);      % battery discharging power [kW]
Ppv = zeros(Tf,1);       % PV generation [kW]
Pi = zeros(Tf,1);        % Imported power [kW]
Pscpv = zeros(Tf,1);     % Self-consumed PV power [kW]
Pcurt = zeros(Tf,1);     % PV Curtailement [kW]

Pabp = zeros(Tf,1);      % Pabp Power Consumption [kW]
T = zeros(Tf+1,1);       % Internal air temperature [°C]
T(1) = T0;               % Initial Internal air temperature [°C]
SoC = zeros(Tf+1,1);     % Battery state of charge [°C]
SoC(1) = SoC0;           % Initial Battery state of charge [°C]

UR_hvac = zeros(Tf,1);   % User Requirements for HVAC

UR_abp = zeros(Tf,1);    % User Requirements for ABP 
%d1 = 0;                  % day counter for PEV (to define user requirements)
d2 = 0;                  % day counter for ABP (to define user requirements)

% intialize variables to manage ABP
abp_varsk.tabp_donek =  zeros(parms.n_abp_phases,1); 
abp_varsk.s_abpk = zeros(parms.n_abp_phases,1); 
abp_varsk.d_abpk = zeros(parms.n_abp_phases,1);
abp_varsk.Eabp_donek =  parms.Eabp; 
abp_varsk.Tabp_donek =  parms.Tabp; 

rng(11)

ur_hvac = [zeros(6,1);ones(12,1);zeros(6,1)];
ur_hvac = repmat(ur_hvac,4,1);

k_start_abp = randi(5)+18; % abp starting time
for k=1:Tf

    % get measurement
    Tk = T(k);
    SoCk = SoC(k);
    
    % get forecasts
    ck = c(k:k+parms.T-1);
    T_ex_forecast_k = [T_ex(k,2);
                       T_ex(k+1:k+parms.T-1,3)];
    P_PV_forecast_k = parms.Pnom_PV*Ir(k:k+parms.T-1,2);
    
    % get user requirements
    
    UR_hvac_k = ur_hvac(k:k+parms.T-1);
    

    UR_abp_k = zeros(parms.T,1);
    if k >= k_start_abp+d2*24
        UR_abp_k(1:k_start_abp+d2*24+12-k) = ones(k_start_abp+12+d2*24-k,1); % complete the cycle in 12 hours
        if k==k_start_abp+d2*24
            abp_varsk.Eabp_donek = zeros(parms.n_abp_phases,1);
            abp_varsk.Tabp_donek = zeros(parms.n_abp_phases,1); 
        end
        if k_start_abp+12+d2*24-k == 0
            d2 = d2+1;
            k_start_abp = randi(5)+18;
        end
    end
    if abp_varsk.Tabp_donek(end)==parms.Tabp(end) % abp cycle termined 
        abp_varsk.tabp_donek = zeros(parms.n_abp_phases,1); 
        abp_varsk.s_abpk = zeros(parms.n_abp_phases,1); 
        abp_varsk.d_abpk = zeros(parms.n_abp_phases,1); 
    end
    
    % compute control 
    [Pc(k),Ph(k),hat_Pchk,hat_Pdsck,hat_Pik,Pabp(k),abp_varsk] = compute_control_step(parms,T_ex_forecast_k,P_PV_forecast_k,Tk,SoCk,ck,UR_hvac_k,UR_abp_k,abp_varsk);
    
    
    % Compute battery limits
    temp_Pb_max = parms.Pnom_b;
    temp_Pb_min = -parms.Pnom_b;
    SoC_k1 = SoC(k) + parms.Dts/parms.Eb*(parms.eta_ch*hat_Pchk-1/parms.eta_dsc*hat_Pdsck);
    if SoC_k1 > 1 % charge command is wrong and it should be limited
        hat_Pchk = (1-SoC(k))*parms.Eb/(parms.Dts*parms.eta_ch);
        temp_Pb_max = (1-SoC(k))*parms.Eb/(parms.Dts*parms.eta_ch);
    elseif SoC_k1 < 0 % discharge command is wrong and it should be limited
        hat_Pdsck = SoC(k)*parms.Eb/(parms.Dts*parms.eta_dsc);
        temp_Pb_min = -SoC(k)*parms.Eb/(parms.Dts*parms.eta_dsc);
    else
        if SoC(k) + parms.Dts/parms.Eb*parms.eta_ch*parms.Pnom_b > 1 % charge cannot be nominal
            temp_Pb_max = (1-SoC(k))*parms.Eb/(parms.Dts*parms.eta_ch);
        end
        if SoC(k) - parms.Dts/parms.Eb*1/parms.eta_dsc*parms.Pnom_b < 0 % discharge cannot be nominal
            temp_Pb_min = -SoC(k)*parms.Eb/(parms.Dts*parms.eta_dsc);
        end
    end
    
    
    % simulate real system
    Ppv(k) = parms.Pnom_PV*Ir(k,2);
    Pi(k) = Pc(k) + Ph(k) + hat_Pchk - hat_Pdsck - Ppv(k) + Pabp(k);
    
    if battery_compensation
        if Pi(k) ~= hat_Pik % forecast error must be compensated
            Pbk = hat_Pik -Pc(k) - Ph(k) + Ppv(k) - Pabp(k);% Try to compensate with battery
            Pbk = max(Pbk,temp_Pb_min); % Saturate battery
            Pbk = min(Pbk,temp_Pb_max);
            Pch(k) = max(0,Pbk);
            Pdsc(k) = -min(0,Pbk);
            Pi(k) = Pc(k) + Ph(k) + Pch(k) - Pdsc(k) - Ppv(k) + Pabp(k);
        else
            Pch(k) = hat_Pchk;
            Pdsc(k) = hat_Pdsck;
        end
    else
        Pch(k) = hat_Pchk;
        Pdsc(k) = hat_Pdsck;
    end
    
    if Pi(k)<0
        Pcurt(k) = -Pi(k);
        Pscpv(k) = Ppv(k)-Pcurt(k);
        Pi(k) = 0;
    else
        Pscpv(k) = Ppv(k);
    end
    T(k+1) = parms.alpha*T(k)-parms.beta*parms.R*(parms.eta_c*Pc(k)-parms.eta_h*Ph(k))+parms.beta*T_ex(k,2);
    SoC(k+1) = SoC(k) + parms.Dts/parms.Eb*(parms.eta_ch*Pch(k)-1/parms.eta_dsc*Pdsc(k));
    
    % update ABP vars
    abp_varsk.Tabp_donek =  abp_varsk.Tabp_donek +  abp_varsk.d_abpk;
    abp_varsk.tabp_donek =  abp_varsk.tabp_donek +  abp_varsk.t_abpk;
    abp_varsk.Eabp_donek =  abp_varsk.Eabp_donek + parms.Dts*Pabp(k)*abp_varsk.d_abpk;
    
    % save current user requirements
    UR_hvac(k) = UR_hvac_k(1);
    UR_abp(k) = UR_abp_k(1);
    
end



%% Plot results
close all

% Temperatures
figure(1)
subplot(3,1,1)
plot(0:Tf,ones(Tf+1,1)*parms.Tsp,':','LineWidth',2)
hold on
plot(0:Tf,ones(Tf+1,1)*(parms.Tsp+parms.Delta),'LineWidth',1.5)
plot(0:Tf,ones(Tf+1,1)*(parms.Tsp-parms.Delta),'LineWidth',1.5)
plot(0:Tf,T,'k','LineWidth',1.5)
plot(0:Tf-1,T_ex(1:Tf,2),'--','LineWidth',1.5)
plot(0:Tf-1,T_ex(1:Tf,3),':','LineWidth',1.5)
xlabel('Time [h]')
ylabel('Temperature [°C]')
xlim([0 Tf])
grid on
legend('Set-point','Max Temp','Min Temp','Internal Temp','External Temp','Forecasted External Temp')
% Power HVAC
subplot(3,1,2)
stairs(0:Tf-1,Pc,'LineWidth',2)
hold on
stairs(0:Tf-1,Ph,':','LineWidth',1.5)
plot(0:Tf-1,parms.Pnom_hvac*ones(Tf,1),'LineWidth',1.5)
xlabel('Time [h]')
ylabel('Power [kW]')
xlim([0 Tf-1])
grid on
legend('Cooling Power','Heating Power','P^{nom}')
% UR HVAC
subplot(3,1,3)
stairs(0:Tf-1,UR_hvac,'LineWidth',1.5)
xlabel('Time [h]')
xlim([0 Tf-1])
grid on
legend('HVAC User Requirements')

% SOC
figure(2)
subplot(2,1,1)
plot(0:Tf,ones(Tf+1,1)*parms.SoCmax,'LineWidth',1.5)
hold on
plot(0:Tf,ones(Tf+1,1)*parms.SoCmin,'LineWidth',1.5)
plot(0:Tf,SoC,'k','LineWidth',1.5)
xlabel('Time [h]')
ylabel('SoC [pu]')
xlim([0 Tf])
ylim([-0.1 1.1])
grid on
legend('SoC max','SoC min','SoC')
% SOC Power
subplot(2,1,2)
stairs(0:Tf-1,Pch-Pdsc,'LineWidth',1.5)
hold on
plot(0:Tf-1,parms.Pnom_b*ones(Tf,1),'LineWidth',1.5)
plot(0:Tf-1,-parms.Pnom_b*ones(Tf,1),'LineWidth',1.5)
xlabel('Time [h]')
ylabel('Power [kW]')
xlim([0 Tf-1])
grid on
legend('Battery Power Exchange','Charging limit','Discharging limit')


% ABP Power
figure(5)
subplot(2,1,1)
stairs(0:Tf-1,Pabp,'LineWidth',1.5)
xlabel('Time [h]')
ylabel('Power [kW]')
xlim([0 Tf-1])
grid on
legend('ABP Power Consumption')
% UR ABP
subplot(2,1,2)
stairs(0:Tf-1,UR_abp,'LineWidth',1.5)
xlabel('Time [h]')
xlim([0 Tf-1])
grid on
legend('ABP User Requirements')

% Powers
figure(6)
hold on
plot(0:Tf-1,Pc+Ph,'LineWidth',1.5)
plot(0:Tf-1,Pch-Pdsc,'LineWidth',1.5)
plot(0:Tf-1,-Ppv,'LineWidth',2.5)
plot(0:Tf-1,-Pscpv,'LineWidth',1.5)
plot(0:Tf-1,Pabp,'LineWidth',1.5)
plot(0:Tf-1,Pi,'k','LineWidth',2)
grid on
box on
legend('HVAC Power Consumption','Battery Power Exchange','PV Generation','PV Self-Consumption','ABP Power Consumption','Imported power')
xlim([0 Tf-1])
xlabel('Time [h]')
ylabel('Power [kW]')

% Costs
figure(7)
subplot(3,1,1)
stairs(0:Tf-1,Pi.*c(1:Tf)*parms.Dts,'LineWidth',2.5)
%hold on
%stairs(0:Tf-1,-Pscpv*parms.tdSCPV*parms.Dts,'LineWidth',2.5)
%hold on
%stairs(0:Tf-1,Pi.*c(1:Tf)*parms.Dts-Pscpv*parms.tdSCPV*parms.Dts,'k','LineWidth',1)
xlim([0 Tf-1])
grid on
ylabel('Hourly Energy Cost [€]')
%legend('Imported energy cost','Self-consumed PV tarif discount','Total cost')
legend('Total cost')
subplot(3,1,2)
stairs(0:Tf-1,cumsum(Pi.*c(1:Tf)*parms.Dts-Pscpv*parms.tdSCPV*parms.Dts),'LineWidth',1.5)
xlim([0 Tf-1])
grid on
ylabel('Total Energy Cost [€]')
subplot(3,1,3)
stairs(0:Tf-1,c(1:Tf),'LineWidth',1.5)
hold on
stairs(0:Tf-1,ones(Tf,1)*parms.tdSCPV,'LineWidth',1.5)
xlim([0 Tf-1])
grid on
xlabel('Time [h]')
ylabel('[€/kWh]')
legend('Energy Price','Self-consumed PV tarif discount')



%% MPC step
function [Pck,Phk,Pchk,Pdsck,Pik,Pabpk,abp_varsk] = compute_control_step(parms,T_ex_forecast_k,P_PV_forecast_k,T_k,SoC_k,ck,UR_hvac_k,UR_abp_k,abp_varsk1)
%% AMPL SETUP
%setupAMPL %MATLAB path setup
ampl = AMPL('.\AMPL'); %open AMPL session
ampl.setOption('solver', 'cplex'); %set the solver
ampl.setOption('log_file', 'logfile.log');%make log file
ampl.read('progetto.mod') %read the ampl file (.mod)

%Communication of parameters to AMPL model

% time index
j = (1:parms.T)';
j = num2cell(j);
J = ampl.getSet('J');
J.setValues(j);

% abp phases index
i = (1:parms.n_abp_phases)';
i = num2cell(i);
I = ampl.getSet('I');
I.setValues(i);

% constant parameters
Pnom_i = ampl.getParameter('Pnom_i');
Pnom_i.setValues(parms.Pnom_i);

Pnom_hvac = ampl.getParameter('Pnom_hvac');
Pnom_hvac.setValues(parms.Pnom_hvac);

alpha = ampl.getParameter('alpha');
alpha.setValues(parms.alpha);

beta = ampl.getParameter('beta');
beta.setValues(parms.beta);

eta_c = ampl.getParameter('eta_c');
eta_c.setValues(parms.eta_c);

eta_h = ampl.getParameter('eta_h');
eta_h.setValues(parms.eta_h);

R = ampl.getParameter('R');
R.setValues(parms.R);

Pnom_b = ampl.getParameter('Pnom_b');
Pnom_b.setValues(parms.Pnom_b);

Eb = ampl.getParameter('Eb');
Eb.setValues(parms.Eb);

eta_ch = ampl.getParameter('eta_ch');
eta_ch.setValues(parms.eta_ch);

eta_dsc = ampl.getParameter('eta_dsc');
eta_dsc.setValues(parms.eta_dsc);

Tsp = ampl.getParameter('Tsp');
Tsp.setValues(parms.Tsp);

Delta = ampl.getParameter('Delta');
Delta.setValues(parms.Delta);

SoCmax = ampl.getParameter('SoCmax');
SoCmax.setValues(parms.SoCmax);

SoCmin = ampl.getParameter('SoCmin');
SoCmin.setValues(parms.SoCmin);

Eabp = ampl.getParameter('Eabp');
Eabp.setValues(parms.Eabp);

Tabp = ampl.getParameter('Tabp');
Tabp.setValues(parms.Tabp);

Dabp = ampl.getParameter('Dabp');
Dabp.setValues(parms.Dabp);

Pmax_abp = ampl.getParameter('Pmax_abp');
Pmax_abp.setValues(parms.Pmax_abp);

Pmin_abp = ampl.getParameter('Pmin_abp');
Pmin_abp.setValues(parms.Pmin_abp);

%tdSCPV = ampl.getParameter('tdSCPV');
%tdSCPV.setValues(parms.tdSCPV);

Dts = ampl.getParameter('Dts');
Dts.setValues(parms.Dts);

% time varying parameters
c = ampl.getParameter('c');
c.setValues(ck);

Tex_forecast  = ampl.getParameter('Tex_forecast');
Tex_forecast.setValues(T_ex_forecast_k);

P_PV_forecast  = ampl.getParameter('P_PV_forecast');
P_PV_forecast.setValues(P_PV_forecast_k);

Tk  = ampl.getParameter('Tk');
Tk.setValues(T_k);

SoCk  = ampl.getParameter('SoCk');
SoCk.setValues(SoC_k);
 
UR_hvac  = ampl.getParameter('UR_hvac');
UR_hvac.setValues(UR_hvac_k);

UR_abp  = ampl.getParameter('UR_abp');
UR_abp.setValues(UR_abp_k);

s_abp_k1  = ampl.getParameter('s_abp_k1');
s_abp_k1.setValues(abp_varsk1.s_abpk);

d_abp_k1  = ampl.getParameter('d_abp_k1');
d_abp_k1.setValues(abp_varsk1.d_abpk);

Tabp_done  = ampl.getParameter('Tabp_done');
Tabp_done.setValues(abp_varsk1.Tabp_donek);

t_abp_done  = ampl.getParameter('t_abp_done');
t_abp_done.setValues(abp_varsk1.tabp_donek);

Eabp_done  = ampl.getParameter('Eabp_done');
Eabp_done.setValues(abp_varsk1.Eabp_donek);

% AMPL SOLUTION 
ampl.solve;

% Get control trajectory
Pc=ampl.getVariable('Pc'); 
Pc=Pc.getValues;
Pc=Pc.getColumnAsDoubles('Pc.val');

Ph=ampl.getVariable('Ph'); 
Ph=Ph.getValues;
Ph=Ph.getColumnAsDoubles('Ph.val');

Pch=ampl.getVariable('Pch'); 
Pch=Pch.getValues;
Pch=Pch.getColumnAsDoubles('Pch.val');

Pdsc=ampl.getVariable('Pdsc'); 
Pdsc=Pdsc.getValues;
Pdsc=Pdsc.getColumnAsDoubles('Pdsc.val');

Pi=ampl.getVariable('Pi'); 
Pi=Pi.getValues;
Pi=Pi.getColumnAsDoubles('Pi.val');

Pabp=ampl.getVariable('Pabp_tot'); 
Pabp=Pabp.getValues;
Pabp=Pabp.getColumnAsDoubles('Pabp_tot.val');

s_abp=ampl.getVariable('s_abp'); 
s_abp=s_abp.getValues;
s_abp=s_abp.getColumnAsDoubles('s_abp.val');
s_abp=reshape(s_abp,parms.n_abp_phases,parms.T); % Consider that the order is inverted as the usual assumed for matrices

d_abp=ampl.getVariable('d_abp'); 
d_abp=d_abp.getValues;
d_abp=d_abp.getColumnAsDoubles('d_abp.val');
d_abp=reshape(d_abp,parms.n_abp_phases,parms.T); 

t_abp=ampl.getVariable('t_abp'); 
t_abp=t_abp.getValues;
t_abp=t_abp.getColumnAsDoubles('t_abp.val');
t_abp=reshape(t_abp,parms.n_abp_phases,parms.T); 

% Apply receding horizon principle
Pck = Pc(1); 
Phk = Ph(1); 
Pchk = Pch(1); 
Pdsck = Pdsc(1); 
Pik = Pi(1); 
Pabpk = Pabp(1);
abp_varsk = abp_varsk1;
abp_varsk.s_abpk = s_abp(:,1); 
abp_varsk.d_abpk = d_abp(:,1);  
abp_varsk.t_abpk = t_abp(:,1); 

ampl.close(); % close the AMPL Session

end

