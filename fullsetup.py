# Calculating Mass Flow Rate 
import math
import csv

#from decimal import *

#clear all; close all;

######################################Containment constants####################

#Givens
#Units are given in square brackets.
R_He = 2076.9; # Helium gas constant [J/Kg/K]
R_Air = 287.05; # Air gas constant [J/Kg/K]
c_pHe = 5190.1; # Specific Heat capacity under constant P for He [J/Kg/K]
c_vHe = 3122.0; # Specific Heat capacity under constant V for He[J/Kg/K]
c_pAir = 1.01e3; # Specific Heat capacity under constant P for He [J/Kg/K]
c_vAir = 718; # Specific Heat capacity under constant V for He[J/Kg/K]

c_wPV = 500; # Specific Heat capacity for PVessel and internals [J/Kg/K]
V_wPV = 0.018; # Volume of Pvessel structure & internals [m^3] 
rho_wPV = 8000;# density of Pvessel structure and internals [kg/m^3]

As_PV= 1.413;# PVTotal surface area [m^2]
#V_PV=265.9;
V_PV = 0.120;# Vessel free volume [m^3]
V_PVC = 0.232;# Vessel free volume [m^3]

Po = 2.5e8;# steady stay reactor power [W]

gamma_PV = 5.0/3; # ratio of specific heats for He [*unitless*]
gamma_PVC = c_pAir/c_vAir; # ratio of specific heats for Air [*unitless*]
T_amb = 300.0;# Ambiant or confinment Temp [K]
rho_amb = 0.13359;# ambiant He density[kg/m^3]
mu_amb = 1.96e-5;# ambiant He viscosity [kg/m/s or Pa*s]
a_PV = math.sqrt(gamma_PV*R_He*T_amb); #speed of sound at exit of pipe [m/s]
a_PVC = math.sqrt(gamma_PVC*R_Air*T_amb); #speed of sound at exit of pipe [m/s]
p_amb = 101325.0; #ambiant pressure [Pa]

mu_i_PV = 3.84e-5;#initial He viscosity [Pa*s]
m_i_PV = 0.173;# intial He Mass in vessel [kg]
mu_i_PVC = 35.47e-6;#initial Air viscosity [Pa*s]
m_i_PVC = 0.1059;# intial Air Mass in vessel [kg]
rho_i_PV = m_i_PV/V_PV;# intial helium density [kg/m^3]
rho_i_PVC = m_i_PVC/V_PVC;# intial Air density [kg/m^3]
p_i_PV = 1.724e6;# initial vessel pressure [Pa]
p_i_PVC = p_amb;# initial PVC pressure [Pa]

Ti_PV = 773.0; # initial helium temp [K]
To_PV = 273.0; # absolute temp [K]
Twi_PV = 773.0;# initial wall temp [K]
Ti_PVC = 500.0; # initial helium temp [K]
To_PVC = 273.0; # absolute temp [K]
Twi_PVC = 500.0;# initial wall temp [K]

eps = 4.5e-5;# pipe equivalant roughness for a new steel pipe [m]
L_PV = 0.0127; # length of pipe [m]
L_PVC = 0.0; # length of pipe [m]
d_PV = 0.0127;# diameter of hot duct [m]
d_PVC = 0.0929;# diameter of hot duct [m]
# d2 = 1.5;
# d3 = 2.29;

Ac_PV = math.pi/4*(d_PV**2);# flow area of duct [m^2]
Pw_PV =math.pi*(d_PV);# Wetted perimeter of duct [m]
D_PV = 4*Ac_PV/Pw_PV;#hydrolic diameter [m_PV]
Ac_PVC = math.pi/4*(d_PVC**2);# flow area of duct [m^2]
Pw_PVC =math.pi*(d_PVC);# Wetted perimeter of duct [m]
D_PVC = (24.0/13.0)/39.37;#hydrolic diameter [m]

K_c = 0.5;# Loss coefficent for sudden contraction [unitless]
f_g = 0.002;#friction factor guess []
alpha = 1;# kentic energy coefficents []

tol=1e-6;#tolerance
t_max = 10.0;# duration of blowdown transient [s]
delta_t =0.001;#time step [s]
n = int(math.floor(t_max/delta_t));# number of time steps
err_PV = 1;
err1_PV = 1;



h = 2;# heat transfer coefficent [W/m^2/K]
#h = 10;
#h=40;

t=0 ;
#PV initial conditions
T_PV = Ti_PV;
Tw_PV = Twi_PV;
p_PV = p_i_PV;
m_PV = m_i_PV;
rho_PV = rho_i_PV;
mu_PV = mu_i_PV;
rho_avg_PV=0;
mu_avg_PV = 0;
v_g_PV = 0 ;
Re_g_PV = 0 ;
f1_PV = 0 ;
f_PV = 0 ;
v_PV= 0;
m_dot_PV=0 ;
Re_PV=0 ;
m_dot_array_PV= [0];
T_array_PV = [Ti_PV];
Tw_array_PV = [Twi_PV];
p_array_PV = [p_i_PV];
v_g_array_PV = [0];
v_array_PV = [0];


#PVC initial conditions
T_PVC = Ti_PVC;
Tw_PVC = Twi_PVC;
p_PVC = p_i_PVC;
m_PVC = m_i_PVC;
rho_PVC = rho_i_PVC;
mu_PVC = mu_i_PVC;
rho_avg_PVC=0;
mu_avg_PVC = 0;
v_g_PVC = 0 ;
Re_g_PVC = 0 ;
f1_PVC = 0 ;
f_PVC = 0 ;
v_PVC= 0;
m_dot_PVC=0 ;
Re_PVC=0 ;
m_dot_array_PVC= [0];
T_array_PVC = [Ti_PVC];
Tw_array_PVC = [Twi_PVC];
p_array_PVC = [p_i_PVC];
v_g_array_PVC = [0];
v_array_PVC = [0];


#User defined
h = 2;

err_PV = 1;
err1_PV = 1;
t_array = [0];

for i in range(1,n):
# T_PV = [];
# Tw_PV = [];
# p_PV = [];
# m_PV = [];
# rho_PV = [];
# mu_PV = [];
# t= ;
# rho_avg_PV=0;
# mu_avg_PV = 0;
# v_g_PV =  ;
# Re_g_PV =  ;
# f1_PV =  ;
# f_PV =  ;
# v_PV= ;

# Re_PV= ;


#User defined

    
#    print "t=  ", t 
    t = ((i-1)*delta_t);
    rho_avg_PV = (rho_PV + rho_amb)/2;
    mu_avg_PV = (mu_PV  + mu_amb)/2;
    rho_avg_PVC = (rho_PVC + rho_amb)/2;
    mu_avg_PVC = (mu_PVC  + mu_amb)/2;

  
    while abs(err_PV) >= tol:
        v_g_PV =  math.sqrt((p_PV - p_amb)/(rho_avg_PV*(alpha / 2.0  + 0.5*(K_c + f_g*L_PV/D_PV))));
	v_g_PVC =  math.sqrt((p_PVC - p_amb)/(rho_avg_PVC*(alpha / 2.0  + 0.5*(K_c + f_g*L_PVC/D_PVC))));
        if v_g_PV  >= a_PV: 
            Re_g_PV = (rho_avg_PV * a_PV * D_PV/mu_avg_PV);
        else :
            Re_g_PV = (rho_avg_PV * v_g_PV  * D_PV/mu_avg_PV);
        #end
	if v_g_PVC  >= a_PVC: 
            Re_g_PVC = (rho_avg_PVC * a_PVC * D_PVC/mu_avg_PVC);
        else :
            Re_g_PVC = (rho_avg_PVC * v_g_PVC  * D_PVC/mu_avg_PVC);
        #end


        if Re_g_PV  >= 4000:
            f1_PV = ( 0.184*Re_g_PV  **(-0.2));
            while abs(err1_PV) >=tol:
                f_PV = ( (1/(-2*math.log(eps/D_PV)/3.7 + 2.51/(Re_g_PV *math.sqrt(f1_PV)))) **2);
                err1_PV = f1_PV  - f_PV ;
                f1_PV  = f_PV ;

            #end
        else:
               f_PV  = 64.0/Re_g_PV ;
        #end
                
        #print "t=  ", t 
        v_PV = (math.sqrt((p_PV - p_amb)/(rho_avg_PV *(alpha/2.0 + 0.5*(K_c + f_PV *L_PV/D_PV)))));
        err_PV = abs(v_g_PV  - v_PV) ;
        f_g_PV = f_PV ;

        
    #end

 	if Re_g_PVC  >= 4000:
            f1_PVC = ( 0.184*Re_g_PV  **(-0.2));
            while abs(err1_PVC) >=tol:
                f_PVC = ( (1/(-2*math.log(eps/D_PVC)/3.7 + 2.51/(Re_g_PVC *math.sqrt(f1_PVC)))) **2);
                err1_PVC = f1_PVC  - f_PVC ;
                f1_PVC  = f_PVC ;

            #end
        else:
               f_PVC  = 64.0/Re_g_PVC ;
        #end
                
        #print "t=  ", t 
        v_PVC = (math.sqrt((p_PVC - p_amb)/(rho_avg_PVC *(alpha/2.0 + 0.5*(K_c + f_PV *L_PVC/D_PVC)))));
        err_PVC = abs(v_g_PVC  - v_PVC) ;
        f_g_PVC = f_PVC ;

        
    #end


    if v_PV  >= a_PV:
        m_dot_PV = ( rho_avg_PV * Ac_PV * a_PV);
        Re_PV = ( rho_avg_PV * a_PV * D_PV / mu_avg_PV);
        v_PV = a_PV;
    else :
        m_dot_PV = ( rho_avg_PV * Ac_PV * v_PV);
        Re_PV = ( rho_avg_PV * v_PV * D_PV / mu_avg_PV);
    #end

    if v_PVC  >= a_PVC:
        m_dot_PVC = ( rho_avg_PVC * Ac_PVC * a_PVC);
        Re_PVC = ( rho_avg_PVC * a_PVC * D_PVC / mu_avg_PVC);
        v_PVC = a_PVC;
    else :
        m_dot_PVC = ( rho_avg_PVC * Ac_PVC * v_PVC);
        Re_PVC = ( rho_avg_PVC * v_PVC * D_PVC / mu_avg_PVC);
    #end


    T_PV = ( delta_t * h * As_PV / (c_vHe * m_PV )*Tw_PV  + (1 - delta_t*(h * As_PV + m_dot_PV *R_He)/(c_vHe * m_PV ))*T_PV );
    Tw_PV = (  delta_t * 0.066 * Po *((i-1)*delta_t + (delta_t **2) **(-0.2))/(rho_wPV * V_wPV *c_wPV) + delta_t * h * As_PV/(rho_wPV * V_wPV *c_wPV)*T_PV  + ( 1 - delta_t * h * As_PV / (rho_wPV * V_wPV *c_wPV))*Tw_PV);
    t = ( t  + delta_t);
    m_PV = ( m_PV  - m_dot_PV  * delta_t);
    mu_PV = ( 1.865e-5 * (T_PV/To_PV) **(0.7));
    rho_PV = ( rho_PV  - m_dot_PV  * delta_t/V_PV);
    p_PV = (rho_PV  * R_He * T_PV);
#### Stuff needs to change as He flows into the PVC.
    T_PVC = ( delta_t * h * As_PVC / (c_vAir * m_PVC  )*Tw_PVC  + (1 - delta_t*(h * As_PVC + m_dot_PVC *R_Air)/(c_vAir * m_PVC ))*T_PVC );

    Tw_PVC = ( delta_t * 0.066 * Po *((i-1)*delta_t + (delta_t **2) **(-0.2))/(rho_wPVC * V_wPVC *c_wPVC) + delta_t * h * As_PVC/(rho_wPVC * V_wPVC *c_wPVC)*T_PVC  + ( 1 - delta_t * h * As_PVC / (rho_wPVC * V_wPVC *c_wPVC))*Tw_PVC);
   
    m_PVC = ( m_PVC  - m_dot_PVC  * delta_t);
    mu_PVC = ( 1.865e-5 * (T_PVC/To_PVC) **(0.7));
    rho_PVC = ( rho_PVC  - m_dot_PVC  * delta_t/V_PVC);
    p_PVC = (rho_PVC  * R_Air * T_PVC);
    
    #store wanted values
    t_array.append(t);
    m_dot_array_PV.append(m_dot_PV);
    T_array_PV.append(T_PV);
    Tw_array_PV.append(Tw_PV);
    p_array_PV.append(p_PV);
    v_g_array_PV.append(v_g_PV);
    v_array_PV.append(v_PV);

    m_dot_array_PVC.append(m_dot_PVC);
    T_array_PVC.append(T_PVC);
    Tw_array_PVC.append(Tw_PVC);
    p_array_PVC.append(p_PVC);
    v_g_array_PVC.append(v_g_PVC);
    v_array_PVC.append(v_PVC);
    
    if p_PV  <= p_amb:
        break
    #end
    err_PV = err_PV + 1;
    err1_PV += 1;
    err_PVC += 1;
    err1_PVC += 1;



#end		
			
			

#h2 = 10;
#h3 = 40;

#(t,T_PV,Tw_PV,m_dot1) = mass_transfer(h)
#(T2,Tw2) = mass_transfer(h2)
#(T3,Tw3) = mass_transfer(h3)

# print "T_PV =  ", T_PV 
   
#print "Re_g_PV =  ", Re_g_PV 
#print "t=  ", t
# print "Tw_PV =  ", Tw_PV
# for i in m_dot_array_PV:
#     print "m_dot_PV =  ", i
# for i in t_array:
#     print "t =  ", i
# for i in T_array_PV:
#     print "T_PV =  ", i
# for i in Tw_array_PV:
#     print "tw =  ", i
# print "Length of m_dot_PV array = ", len(m_dot_array_PV)
# print "Length of t array = ", len(t_array)
# print "Length of T_PV array = ", len(T_array_PV)
# print "Length of Tw_PV array = ", len(Tw_array_PV)
# print "Length of p_PV array = ", len(p_array_PV)
# print "Length of v_g_PV array = ", len(v_g_array_PV)
#print "t[333] = ", t_array[333]


#write data to csv file
PVmdotdata = [t_array,m_dot_array_PV];
PVvgdata = [t_array,v_g_array_PV];
PVvdata = [t_array,v_array_PV];
PVpdata = [t_array,p_array_PV];

PVCmdotdata = [t_array,m_dot_array_PVC];
PVCvgdata = [t_array,v_g_array_PVC];
PVCvdata = [t_array,v_array_PVC];
PVCpdata = [t_array,p_array_PVC];

with open('mdot_PV.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVmdotdata)

csvFile.close()

with open('vg_PV.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVvgdata)
    
csvFile.close()

with open('p_PV.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVpdata)
    
csvFile.close()

with open('v_PV.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVvdata)
    
csvFile.close()



with open('mdot_PVC.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVCmdotdata)

csvFile.close()

with open('vg_PVC.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVCvgdata)
    
csvFile.close()

with open('p_PVC.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVCpdata)
    
csvFile.close()

with open('v_PVC.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(PVCvdata)
    
csvFile.close()




# print "rho_avg_PV =  ", rho_avg_PV

print "max of v_g_array_PV =",max(v_g_array_PV)
print "max of v_array_PV =",max(v_array_PV)
# print "max of t_array =",max(t_array)
print "a_PV = ", a_PV
			
			
			
