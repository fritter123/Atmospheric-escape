clc;
clear;
close all;
%% Input parameters
M=59800; %E23 g
r=6351; %km
d=1; %orbital radius in AU
%Planet

P_H2O=130; %partial pressure of H2O in bar
P_CO2=12; %partial pressure of CO2 in bar
nebula_H = 0.00001; % (wt%) weight percent of nebula H in wt%(mass of H (captured from bebula at begining) in atmosphere / mass of Earth)
%Initial atmosphere

RT=2; %Type of rotator of the young Sun. 1,2 and 3 corresponds to slow, moderate and fast rotator
eta=0.15; %EUV heating efficiency
%EUV absorption

rmax=2.55*r*100000;
rmin=1.05*r*100000;
r0=0.7*rmax;
T=300; %atmospheric temperature
%Base of the hydrodynamic flow

t=50;%when escape starts
delta_time=0.01;%time step in Myr
%Starting time and time step

%% initialization
t_0=t;
Tbas=T;
beta=1;
gr=(6.67259*10^(-8)*M*10^23)/(r*100000)^2;

m_H=1.008*1.993/12*10^(-23);
m_O=16*1.993/12*10^(-23);

Amount_H_nebula=nebula_H/100*M*10^23/m_H;
H_1bar = (1 * 10 ^ 5 * 10 ^ (-4)) / (gr / 100) * 4 * pi * (r * 10 ^ 5) ^ 2 * 10 ^ 3 / m_H;

Mol_CO2=44.010;
m_CO2=Mol_CO2*1.993/12*10^(-23);
m_20Ne=19.99*1.993/12*10^(-23);
m_22Ne=21.99*1.993/12*10^(-23);
m_36Ar=35.97*1.993/12*10^(-23);
m_38Ar=37.96*1.993/12*10^(-23);

M_H2O=(P_H2O*10^5*10^(-4))/(gr/100)*4*pi*(r*10^5)^2*10^3;
M_CO2=(P_CO2*10^5*10^(-4))/(gr/100)*4*pi*(r*10^5)^2*10^3;
Amount_O=M_H2O/(2*m_H+m_O);

Amount_H=2*Amount_O+Amount_H_nebula;
Amount_CO2=M_CO2/m_CO2;
Amount_Ne_20=1.07e-4*Amount_H_nebula;
Amount_Ar_36=2.69e-6*Amount_H_nebula;
Amount_Ne_22=Amount_Ne_20/13.709;
Amount_Ar_38=Amount_Ar_36/5.804;

Amount_O_0=Amount_O;
Amount_H_0=Amount_H;
Amount_CO2_0=Amount_CO2;
Amount_Ar_36_0=Amount_Ar_36;
Amount_Ar_38_0=Amount_Ar_38;
Amount_Ne_20_0=Amount_Ne_20;
Amount_Ne_22_0=Amount_Ne_22;

f0=Amount_O/Amount_H;%note that 0 in 'f0' is written as number zero，represent for O (oxygen)
fco2=Amount_CO2/Amount_H;
fAr_36=Amount_Ar_36/Amount_H;
fAr_38=Amount_Ar_38/Amount_H;
fNe_20=Amount_Ne_20/Amount_H;
fNe_22=Amount_Ne_22/Amount_H;
time_gap=3.1536e+13*delta_time;   %  1Myr（in s）*delta_time
H_particle=[];
O_particle=[];
CO2_particle=[];
Ar_36_particle=[];
Ar_38_particle=[];
Ne_20_particle=[];
Ne_22_particle=[];

kB=1.380649e-23*1e7; %erg/K
Time=[];

mc_O=1;
mc_co2=1;
i=1;
t_total=t_0;
Amount_min=1;

%% calculation

while Amount_H > 1 * H_1bar && t_total < (t_0+2000) && Amount_min > 0
    
    if RT == 1
        if t>5.7
            FEUV=7.4*1e31/(4*pi*(1.496*1e13)^2)*(t^(-0.96)/d^2);  % erg.cm-2.s-1
        else
            FEUV=7.4*1e31/(4*pi*(1.496*1e13)^2)*(5.7^(-0.96)/d^2);
        end
        %slow rotator
    elseif RT == 2
        if t>23
            FEUV=4.8*1e32/(4*pi*(1.496*1e13)^2)*(t^(-1.22)/d^2);  % erg.cm-2.s-1
        else
            FEUV=4.8*1e32/(4*pi*(1.496*1e13)^2)*(23^(-1.22)/d^2);
        end
        %moderate rotator
    elseif RT ==3
        if t>226
            FEUV=1.2*1e36/(4*pi*(1.496*1e13)^2)*(t^(-2.15)/d^2);  % erg.cm-2.s-1
        else
            FEUV=1.2*1e36/(4*pi*(1.496*1e13)^2)*(226^(-2.15)/d^2);
        end
        %fast rotator
    end

    beta=1;
    delta_fei=(6.67259*10^(-8)*M*10^23)/r0;
    gr=(6.67259*10^(-8)*M*10^23)/(r*100000)^2;
    g0=gr*(((r*100000)/r0)^2);
    matm=1000000/gr;%（4+2：cm2→m2 + 980→9.8）
    H_1bar = (1 * 10 ^ 5 * 10 ^ (-4)) / (gr / 100) * 4 * pi * (r * 10 ^ 5) ^ 2 * 10 ^ 3 / m_H;
    
    if m_CO2 < mc_co2
        alpha_1 = (g0 * (m_O - m_H) * 4.8e17 * Tbas ^ 0.75 * 4 * delta_fei) / (kB * Tbas * (1 + f0) * beta ^ 2 * eta * FEUV);
        alpha_2 = (g0 * (m_CO2 - m_H) * 8.4e17 * Tbas ^ 0.6 * 4 * delta_fei) / (kB * Tbas * (1 + 8.4e17 / 7.86e16 * Tbas ^ (-0.176) * f0) * beta ^ 2 * eta * FEUV);
        beta_1 = (f0 * (8.4e17 / 7.86e16 * Tbas ^ (-0.176) - 8.4e17 / 4.8e17 * Tbas ^ (-0.15))) / (1 + 8.4e17 / 7.86e16 * Tbas ^ (-0.176) * f0);
        gamma = (1 + 8.4e17 / 4.8e17 * Tbas ^ (-0.15) * f0) / (1 + 8.4e17 / 7.86e16 * Tbas ^ (-0.176) * f0) - alpha_2 * m_H;
        
        x0 = (1 - alpha_1 * m_H - (alpha_1 * m_CO2 * fco2 * gamma) / (1 + alpha_2 * m_CO2 * fco2)) / (1 + alpha_1 * m_O * f0 + alpha_1 * m_CO2 * fco2 * (beta_1 - alpha_2 * m_O * f0) / (1 + alpha_2 * m_CO2 * fco2));
        xco2 = (gamma + (beta_1 - alpha_2 * m_O * f0) * (1 - -alpha_1 * m_H) / (1 + alpha_1 * m_O * f0)) / (1 + alpha_2 * m_CO2 * fco2 + (beta_1 - alpha_2 * m_O * f0) * (alpha_1 * m_CO2 * fco2) / (1 + alpha_1 * m_O * f0));
        FH = (beta ^ 2 * eta * FEUV) / (4 * delta_fei * (m_H + m_O * f0 * x0 + m_CO2 * fco2 * xco2));
    else
        alpha_1 = (g0 * (m_O - m_H) * 4.8e17 * Tbas ^ 0.75 * 4 * delta_fei) / (kB * Tbas * (1 + f0) * beta ^ 2 * eta * FEUV);
        x0 = (1 - alpha_1 * m_H) / (1 + alpha_1 * m_CO2 * fco2);
        xco2 = 0;
        FH = (beta ^ 2 * eta * FEUV) / (4 * delta_fei * (m_H + m_O * f0 * x0));
    end
    
    xAr_36=(1-g0*(m_36Ar-m_H)*1.06e18/1.38*1e16*T^(0.597)/(FH*T)+1.06e18/(4.8e17)*T^(-0.153)*f0*(1-x0)+1.06e18/(5.61e16)*T^(-0.244)*f0*x0)/(1+1.06e18/(5.61e16)*T^(-0.244)*f0);
    xAr_38=(1-g0*(m_38Ar-m_H)*1.06e18/1.38*1e16*T^(0.597)/(FH*T)+1.06e18/(4.8e17)*T^(-0.153)*f0*(1-x0)+1.06e18/(5.61e16)*T^(-0.244)*f0*x0)/(1+1.06e18/(5.61e16)*T^(-0.244)*f0);
    xNe_20=(1-g0*(m_20Ne-m_H)*7.9e17/1.38*1e16*T^(0.731)/(FH*T)+7.9e17/(4.8e17)*T^(-0.019)*f0*(1-x0)+7.9e17/(1.5e17)*f0*x0*T^(-0.019))/(1+7.9e17/(1.5e17)*f0*T^(-0.019));
    xNe_22=(1-g0*(m_22Ne-m_H)*7.9e17/1.38*1e16*T^(0.731)/(FH*T)+7.9e17/(4.8e17)*T^(-0.019)*f0*(1-x0)+7.9e17/(1.5e17)*f0*x0*T^(-0.019))/(1+7.9e17/(1.5e17)*f0*T^(-0.019));
    
    F0= FH*x0*f0;
    
    XH_O=Amount_H/(Amount_H+Amount_O);
    XH_co2=Amount_H/(Amount_H+Amount_CO2);
    XH_Ar_36=Amount_H/(Amount_H+Amount_Ar_36);
    XH_Ar_38=Amount_H/(Amount_H+Amount_Ar_38);
    XH_Ne_20=Amount_H/(Amount_H+Amount_Ne_20);
    XH_Ne_22=Amount_H/(Amount_H+Amount_Ne_22);
    
    mc_O=m_H+kB*T*FH/(4.8e17*T^(0.75)*g0*XH_O);
    mc_co2=m_H+kB*T*FH/(8.4e17*T^(0.6)*g0*XH_co2);
    mc_Ar_36=m_H+kB*T*FH/(1.06e18*T^(0.597)*g0*XH_Ar_36);
    mc_Ar_38=m_H+kB*T*FH/(1.06e18*T^(0.597)*g0*XH_Ar_38);
    mc_Ne_20=m_H+kB*T*FH/(7.9e17*T^(0.731)*g0*XH_Ne_20);
    mc_Ne_22=m_H+kB*T*FH/(7.9e17*T^(0.731)*g0*XH_Ne_22);
    
    if m_36Ar>mc_Ar_36 || xAr_36<0
        xAr_36=0;
    end
    
    if m_38Ar>mc_Ar_38 || xAr_38<0
        xAr_38=0;
    end
    
    if m_20Ne>mc_Ne_20 || xNe_20<0
        xNe_20=0;
    end
    
    if m_22Ne>mc_Ne_22 || xNe_22<0
        xNe_22=0;
    end
    
    FCO2= FH* xco2*fco2;
    FAr_36= FH* xAr_36*fAr_36;
    FAr_38= FH* xAr_38*fAr_38;
    FNe_20= FH* xNe_20*fNe_20;
    FNe_22= FH* xNe_22*fNe_22;
    
    H_particle(i)= Amount_H;
    escape_H=time_gap*4*pi*(r0)^2*FH;
    Amount_H= Amount_H- escape_H;
    Escape_H(i)=escape_H;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    O_particle(i)= Amount_O;
    escape_O=time_gap*4*pi*(r0)^2*F0;
    Amount_O= Amount_O- escape_O;
    Escape_O(i)=escape_O;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CO2_particle(i)= Amount_CO2;
    escape_CO2=time_gap*4*pi*(r0)^2*FCO2;
    Amount_CO2= Amount_CO2- escape_CO2;
    Escape_CO2(i)=escape_CO2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ar_36_particle(i)= Amount_Ar_36;
    escape_Ar_36=time_gap*4*pi*(r0)^2*FAr_36;
    Amount_Ar_36= Amount_Ar_36- escape_Ar_36;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ar_38_particle(i)= Amount_Ar_38;
    escape_Ar_38=time_gap*4*pi*(r0)^2*FAr_38;
    Amount_Ar_38= Amount_Ar_38- escape_Ar_38;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ne_20_particle(i)= Amount_Ne_20;
    escape_Ne_20=time_gap*4*pi*(r0)^2*FNe_20;
    Amount_Ne_20= Amount_Ne_20- escape_Ne_20;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ne_22_particle(i)= Amount_Ne_22;
    escape_Ne_22=time_gap*4*pi*(r0)^2*FNe_22;
    Amount_Ne_22= Amount_Ne_22- escape_Ne_22;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Amount_min_1 = min(Amount_O, Amount_CO2);
    Amount_min_2 = min(Amount_Ar_36,Amount_Ar_38);
    Amount_min_3 = min(Amount_Ne_20,Amount_Ne_22);
    Amount_min_4 = min(Amount_min_1,Amount_min_2);
    Amount_min = min(Amount_min_4,Amount_min_3);
    
    f0=Amount_O/Amount_H;
    
    fco2=Amount_CO2/Amount_H;
    fAr_36=Amount_Ar_36/Amount_H;
    fAr_38=Amount_Ar_38/Amount_H;
    fNe_20=Amount_Ne_20/Amount_H;
    fNe_22=Amount_Ne_22/Amount_H;
    
    H_inpressure=H_particle.*1.008./(matm*6.022e23*4*pi*(r*100000)^2);
    O_inpressure=O_particle.*15.999./(matm*6.022e23*4*pi*(r*100000)^2);
    CO2_inpressure=CO2_particle.*Mol_CO2./(matm*6.022e23*4*pi*(r*100000)^2);
    Time(i)=t;
    t=t+delta_time;
    
    i=i+1;
    t_total=t_total+delta_time;
    
end

ratio_Ne=Ne_20_particle./Ne_22_particle;
ratio_Ar=Ar_36_particle./Ar_38_particle;

%% Figure
LW = 1.5;
LW2 = 1.2;
plot(Time,O_inpressure,'LineWidth',LW,'Color',[0.00,0.78,0.55]);
hold on;
plot(Time,H_inpressure,'LineWidth',LW,'Color',[0.42,0.35,0.80]);
hold on;
plot(Time,CO2_inpressure,'LineWidth',LW,'Color',[0.25,0.88,0.82]);
axis([t_0,t_total,1,10000])
legend('O','H','CO2');
set(gca,'YScale','log')
xlabel('Time (Myr)');
ylabel('P(bar)');
title('evolution of H,O,CO2');
