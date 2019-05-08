# Calculating Mass Flow Rate 
import math
import csv

#from decimal import *

#clear all; close all;

#############################Containment constants#############################



#Givens
#Units are given in square brackets.



############################# Needed Constants ###############################

c_w = 500; # Specific Heat capacity for Vessel and internals [J/Kg/K]
rho_w = 8000;# density of vessel structure and internals [kg/m_v1^3]
Po = 0.0 #2.5e8;# steady stay reactor power [W]
tol=1e-6;#tolerance
t_max = 10.0;# duration of blowdown transient [s]
delta_t =0.001;#time step [s]
n = int(math.floor(t_max/delta_t));# number of time steps
alpha = 1;# kentic energy coefficents []
K_c = 0.5;# Loss coefficent for sudden contraction [unitless]
eps = 4.5e-5;# pipe equivalant roughness for a_v1 new steel pipe [m]
To = 273.0; # absolute temp [K]
p_amb = 101325.0; #ambiant pressure [Pa]
T_amb = 300.0;# Ambiant or confinment Temp [K]
h = 2;# heat transfer coefficent [W/m_v1^2/K]
t = 0;
t_array = [0];
#h = 10;
#h=40;




########################## For Vessel 1: The PV ########################################

# Gas properties
R_v1 =287.05 # 2076.9; # Helium gas constant [J/Kg/K]
c_p_v1 = 1.01e3 #5190.1; # Specific Heat capacity under constant P for He [J/Kg/K]
c_v_v1 = 718.0 #3122.0; # Specific Heat capacity under constant V_v1 for He[J/Kg/K]

gamma_v1 = c_p_v1/c_v_v1; # ratio of specific heats for He [*unitless*]
rho_amb_v1 = 0.13359;# ambiant He density[kg/m_v1^3]
mu_amb_v1 = 1.96e-5;# ambiant He viscosity [kg/m_v1/s or Pa*s]


# physical properties of vessel
V_w_v1 = 0.02; # Volume of vessel structure & internals [m^3] 
As_v1= 1.413;# Total surface area [m^2]
V_v1 = 0.120;# Vessel free volume [m^3]



# Mass, Mu, Pressure, and Rho
mu_i_v1 = 3.84e-5;#initial He viscosity [Pa*s]
m_i_v1 = 0.173;# intial He Mass in vessel [kg]
rho_i_v1 = m_i_v1/V_v1;# intial helium density [kg/m_v1^3]
p_i_v1 = 1.724e6;# initial vessel pressure [Pa]


# Tempurature
Ti_v1 = 773.15; # initial helium temp [K]
Twi_v1 = 773.15;# initial wall temp [K]


# Duct properties
L_v1 = 0.0127; # length of pipe [m]
d_v1 = 0.0127;# diameter of hot duct [m]
Ac_v1 = math.pi/4*(d_v1**2);# flow area of duct [m^2]
Pw_v1 =math.pi*(d_v1);# Wetted perimeter of duct [m]
D_v1 = 4*Ac_v1/Pw_v1;#hydrolic diameter [m]
f_g_v1 = 0.002;#friction factor guess []
a_v1 = math.sqrt(gamma_v1*R_v1*T_amb); #speed of sound at exit of pipe [m/s]



# Error 
err_v1 = 1;
err1_v1 = 1;


# Initial Conditions
T_v1 = Ti_v1;
T_v1temp = Ti_v1;
Tw_v1 = Twi_v1;
Tw_v1temp = Tw_v1
p_v1 = p_i_v1;
m_v1 = m_i_v1;
m_v1temp = m_v1
rho_v1 = rho_i_v1;
mu_v1 = mu_i_v1;
rho_avg_v1= 0;
mu_avg_v1 = 0;
v_g_v1 = 0 ;
Re_g_v1 = 0 ;
f1_v1 = 0 ;
f_v1 = 0 ;
v_v1= 0;
m_dot_v1=0 ;
Re_v1=0 ;
delta_u_v1 =0;


# Initialize arrays for storage
mdot_array_v1= [0];
T_array_v1 = [Ti_v1];
Tw_array_v1 = [Twi_v1];
p_array_v1 = [p_i_v1];
vg_array_v1 = [0];
v_array_v1 = [0];
du_array_v1 = [delta_u_v1];
m_array_v1 = [m_i_v1]




########################## Vessel 2: The PVC ########################################


# Gas properties
R_v2 = 287.05; # Helium gas constant [J/Kg/K]
c_p_v2 =  1.01e3; # Specific Heat capacity under constant P for He [J/Kg/K]
c_v_v2 = 718.0; # Specific Heat capacity under constant V_v1 for He[J/Kg/K]
gamma_v2 = c_p_v2/c_v_v2; # ratio of specific heats for He [*unitless*]
rho_amb_v2 = 1.225;# ambiant He density[kg/m_v1^3]
mu_amb_v2 = 1.813e-5;# ambiant He viscosity [kg/m_v1/s or Pa*s]
a_v2 = math.sqrt(gamma_v2*R_v2*T_amb); #speed of sound at exit of pipe [m/s]



# Physical properties of vessel
V_w_v2 = 0.02; # Volume of vessel structure & internals [m^3] 
As_v2= 24971.0/1550.02;# Total surface area [m^2]
V_v2 = 0.232;# Vessel free volume [m^3]



# Mass, Mu, Pressure, and Rho
mu_i_v2 = 3.64e-5;#initial He viscosity [Pa*s]
m_i_v2 = 0.173;# intial He Mass in vessel [kg]
rho_i_v2 = 0.4567;# intial helium density [kg/m_v1^3]
p_i_v2 = 101325;# initial vessel pressure [Pa]



# Tempurature
Ti_v2 = 773.15; # initial helium temp [K]
Twi_v2 = 773.15;# initial wall temp [K]


# Duct properties
L_v2 = 0.0; # length of pipe [m]
d_v2 = 0.0929;# diameter of hot duct [m]
Ac_v2 = math.pi/4*(d_v2**2);# flow area of duct [m^2]
Pw_v2 =math.pi*(d_v2);# Wetted perimeter of duct [m]
D_v2 = (24.0/13.0)/39.37;#hydrolic diameter [m]
f_g_v2 = 0.002;#friction factor guess []


# Initial Conditions
T_v2 = Ti_v2;
T_v2temp = Ti_v2;
Tw_v2 = Twi_v2;
Tw_v2temp = Tw_v2;
p_v2 = p_i_v2;
m_v2 = m_i_v2;
m_v2temp = m_v2;
rho_v2 = rho_i_v2;
mu_v2 = mu_i_v2;
rho_avg_v2= 0;
mu_avg_v2 = 0;
v_g_v2 = 0 ;
Re_g_v2 = 0 ;
f1_v2 = 0 ;
f_v2 = 0 ;
v_v2= 0;
m_dot_v2=0 ;
Re_v2=0 ;
delta_u_v2 = 0; 

rho_amb_v1 = rho_i_v2;# ambiant He density[kg/m_v1^3]
mu_amb_v1 = mu_i_v2;# ambiant He viscosity [kg/m_v1/s or Pa*s]


# Error 
err_v2 = 1;
err1_v2 = 1;


# Initialize arrays for storage
mdot_array_v2= [0];
T_array_v2 = [Ti_v2];
Tw_array_v2 = [Twi_v2];
p_array_v2 = [p_i_v2];
vg_array_v2 = [0];
v_array_v2 = [0];
du_array_v2 = [delta_u_v2];
m_array_v2 = [m_i_v2]


################################## Run Blowdown  ############################


for i in range(1,n+1):

    t += (delta_t);# update total time
    
####################### Vessel 1 calculations#########################
    
    a_v1 = math.sqrt(gamma_v1*R_v1*T_v2); # Recalulate new sound speed.
    rho_avg_v1 = (rho_v1 + rho_amb_v1)/2; # Recalculate new rho avg
    mu_avg_v1 = (mu_v1  + mu_amb_v1)/2;
  
    while abs(err_v1) >= tol:
        v_g_v1 =  math.sqrt((p_v1  - p_amb)/(rho_avg_v1*(alpha / 2.0  + 0.5*(K_c + f_g_v1*L_v1/D_v1))));
        if v_g_v1  >= a_v1: 
            Re_g_v1 = (rho_avg_v1 * a_v1 * D_v1/mu_avg_v1);
        else :
            Re_g_v1 = (rho_avg_v1 * v_g_v1  * D_v1/mu_avg_v1);
        #end
        if Re_g_v1  >= 4000:
            f1_v1 = ( 0.184*Re_g_v1  **(-0.2));
            while abs(err1_v1) >=tol:
                f_v1 = ( (1/(-2*math.log(eps/D_v1)/3.7 + 2.51/(Re_g_v1 *math.sqrt(f1_v1)))) **2);
                err1_v1 = f1_v1  - f_v1 ;
                f1_v1  = f_v1 ;

            #end
        else:
               f_v1  = 64.0/Re_g_v1 ;
        #end
                
        #print "t=  ", t 
        v_v1 = (math.sqrt((p_v1 - p_amb)/(rho_avg_v1 *(alpha/2.0 + 0.5*(K_c + f_v1 *L_v1/D_v1)))));
        err_v1 = abs(v_g_v1  - v_v1) ;
        f_g_v1 = f_v1 ;

        
    #end
    if v_v1  >= a_v1:
        m_dot_v1 = ( rho_avg_v1 * Ac_v1 * a_v1);
        Re_v1 = ( rho_avg_v1 * a_v1 * D_v1 / mu_avg_v1);
        v_v1 = a_v1;
    else :
        m_dot_v1 = ( rho_avg_v1 * Ac_v1 * v_v1);
        Re_v1 = ( rho_avg_v1 * v_v1 * D_v1 / mu_avg_v1);
    #end
    
    m_v1 -=  m_dot_v1  * delta_t;
    # T =  delta_t * h * As / (c_v * m_temp )*Twtemp  + (1 - delta_t*(h * As + m_dot *R)/(c_v * m_temp ))*Ttemp ;
    T_v1 =  delta_t * h * As_v1 / (c_v_v1 * m_v1temp )*Tw_v1temp  + (1 - delta_t*(h * As_v1 + m_dot_v1 *R_v1)/(c_v_v1 * m_v1temp ))*T_v1temp ;
    delta_u_v1 +=  m_v1  * c_p_v1 * (T_v1 - T_v1temp)
    Tw_v1 =   (delta_t * 0.066 * Po *( t - delta_t + delta_t ** 2) ** (-0.2))/(rho_w * V_w_v1 *c_w) + delta_t * h * As_v1/(rho_w * V_w_v1 *c_w)*T_v1temp  + ( 1 - delta_t * h * As_v1 / (rho_w * V_w_v1 *c_w))*Tw_v1temp;

    mu_v1 = ( mu_amb_v1 * (T_v1/To) **(0.7));
    rho_v1 = m_v1/V_v1;
    p_v1 = (rho_v1  * R_v1 * T_v1);
    
    T_v1temp = T_v1;
    Tw_v1temp = Tw_v1;
    m_v1temp = m_v1;
    
    ####################### Vessel 2 calculations#########################
    
    
    rho_avg_v2 = (rho_v2 + rho_amb_v2)/2;
    mu_avg_v2 = (mu_v2  + mu_amb_v2)/2;
  
    # if p_v2 < p_amb:
    #     p_v2 = p_amb
        
    v_g_v2 =  math.sqrt((p_v2  - p_amb)/(rho_avg_v2*(alpha + K_c)/2.0 ));

    
    if v_g_v2  > a_v2: 
        v_g_v2 = a_v2;
        
    else :
        pass
    #end
    v_v2 = v_g_v2;
    #end
    v_v2 = 0
    m_dot_v2 = 0#( rho_avg_v2 * Ac_v2 * v_v2);
    #end
    
    m_v2 +=  m_dot_v1 * delta_t ; #- m_dot_v2 * delta_t 
    T_v2 =  delta_t * h * As_v2 / (c_v_v2 * m_v2temp )*Tw_v2temp  + T_v2temp #(1 - delta_t*(h * As_v2 + m_dot_v2 *R_v2)/(c_v_v2 * m_v2temp ))*T_v2temp ;
    delta_u_v2 +=  m_v2  * c_p_v2 * (T_v2 - T_v2temp);
    Tw_v2 = (  delta_t * 0.066 * Po *((i-1)*delta_t + (delta_t **2) **(-0.2))/(rho_w * V_w_v2 *c_w) + delta_t * h * As_v2/(rho_w * V_w_v2 *c_w)*T_v2temp  + ( 1 - delta_t * h * As_v2 / (rho_w * V_w_v2 *c_w))*Tw_v2temp);

    mu_v2 = ( mu_amb_v2 * (T_v2/To) **(0.7));
    rho_v2 = m_v2/V_v2;
    p_v2 = (rho_v2  * R_v2 * T_v2);
    
    T_v2temp = T_v2;
    Tw_v2temp = Tw_v2;
    m_v2temp = m_v2;
    rho_amb_v1 = rho_v2
    mu_amb_v1 = mu_v2
    
    
    
    #store wanted values
    t_array.append(t);
    
    
    mdot_array_v1.append(m_dot_v1);
    m_array_v1.append(m_v1);
    T_array_v1.append(T_v1);
    Tw_array_v1.append(Tw_v1);
    p_array_v1.append(p_v1);
    vg_array_v1.append(v_g_v1);
    v_array_v1.append(v_v1);
    du_array_v1.append(delta_u_v1);
    
    
    mdot_array_v2.append(m_dot_v2);
    m_array_v2.append(m_v2);
    T_array_v2.append(T_v2);
    Tw_array_v2.append(Tw_v2);
    p_array_v2.append(p_v2);
    vg_array_v2.append(v_g_v2);
    v_array_v2.append(v_v2);
    du_array_v2.append(delta_u_v2);
    
    if p_v1  <= p_amb:
        print "p_v1 <= p_amb"
        print "t = ",t
        break
    #end
    if p_v1  <= p_v2:
        print "p_v1 <= p_v2"
        print "t = ",t
        break
    #end
    
    # if p_v2  <= p_amb:
    #     print "p_v2 <= p_amb"
    #     print "t = ",t
    #     break
    #end
    
    if T_v1 <0 or T_v2<0:
        print "error. Temp below 0 K"
        break
    
    if t> t_max:
        print "t > tmax"
        break
    
    # if p_v2  <= p_amb:
    #     print "p_v2 <= p_amb"
    #     print "t = ",t
    #     break
    # # end
    
    # if p_v1  <= p_v2:
    #     print "p_v1 <= p_v2"
    #     print "t = ",t
    #     break
    
    err_v1 = err_v1 + 1;
    err1_v1 += 1;



#end		

# print "T_v1 =  ", T_v1 
   
#print "Re_g_v1 =  ", Re_g_v1 
#print "t=  ", t
# print "Tw_v1 =  ", Tw_v1
# for i in mdot_array_v1:
#     print "m_dot_v1 =  ", i
# for i in t_array_v1:
#     print "t =  ", i
# for i in T_array_v1:
#     print "T_v1 =  ", i
# for i in Tw_array_v1:
#     print "tw =  ", i
# print "Length of m_dot_v1 array = ", len(mdot_array_v1)
# print "Length of t array = ", len(t_array_v1)
# print "Length of T_v1 array = ", len(T_array_v1)
# print "Length of Tw_v1 array = ", len(Tw_array_v1)
# print "Length of p_v1 array = ", len(p_array_v1)
# print "Length of v_g_v1 array = ", len(vg_array_v1)
#print "t[333] = ", t_array_v1[333]


###################### write data to csv file ########################

#                 Vessel 1

mdot_Data_v1 = [t_array,mdot_array_v1];
m_Data_v1 = [t_array,m_array_v1];
du_Data_v1 = [t_array,du_array_v1];
vData_v1 = [t_array,v_array_v1];
pData_v1 = [t_array,p_array_v1];
T_Data_v1 = [t_array,T_array_v1];


with open('mdotV1.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(mdot_Data_v1)

csvFile.close()

with open('duV1.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(du_Data_v1)
    
csvFile.close()

with open('pV1.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(pData_v1)
    
csvFile.close()

with open('vV1.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(vData_v1)
    
csvFile.close()


with open('TV1.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(T_Data_v1)
    
csvFile.close()

with open('mV1.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(m_Data_v1)
    
csvFile.close()


#             Vessel 2

mdot_Data_v2 = [t_array,mdot_array_v2];
du_Data_v2 = [t_array,du_array_v2];
vData_v2 = [t_array,v_array_v2];
pData_v2 = [t_array,p_array_v2];
m_Data_v2 = [t_array,m_array_v2];
T_Data_v2 = [t_array,T_array_v2];

with open('mdotV2.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(mdot_Data_v2)

csvFile.close()

with open('duV2.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(du_Data_v2)
    
csvFile.close()

with open('pV2.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(pData_v2)
    
csvFile.close()

with open('vV2.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(vData_v2)
    
csvFile.close()

with open('TV2.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(T_Data_v2)
    
csvFile.close()

with open('mV2.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(m_Data_v2)
    
csvFile.close()

# print "rho_avg_v1 =  ", rho_avg_v1

#print "max of vg_array_v1 =",max(vg_array_v1)
#print "max of v_array_v1 =",max(v_array_v1)
# print "max of t_array_v1 =",max(t_array_v1)
#print "a_v1 = ", a_v1
			
			
			