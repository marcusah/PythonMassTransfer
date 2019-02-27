# Calculating Mass Flow Rate 
import math
import csv

#from decimal import *

#clear all; close all;

######################################Containment constants####################

#Givens
#Units are given in square brackets.
R = 2076.9; # Helium gas constant [J/Kg/K]
c_p = 5190.1; # Specific Heat capacity under constant P for He [J/Kg/K]
c_v = 3122.0; # Specific Heat capacity under constant V for He[J/Kg/K]

c_w = 500; # Specific Heat capacity for Vessel and internals [J/Kg/K]
V_w = 0.018; # Volume of vessel structure & internals [m^3] 
rho_w = 8000;# density of vessel structure and internals [kg/m^3]

As= 1.413;# Total surface area [m^2]
#V=265.9;
V = 0.120;# Vessel free volume [m^3]

Po = 2.5e8;# steady stay reactor power [W]

gamma = 5.0/3; # ratio of specific heats for He [*unitless*]
T_amb = 300.0;# Ambiant or confinment Temp [K]
rho_amb = 0.13359;# ambiant He density[kg/m^3]
mu_amb = 1.96e-5;# ambiant He viscosity [kg/m/s or Pa*s]
a = math.sqrt(gamma*R*T_amb); #speed of sound at exit of pipe [m/s]
p_amb = 101325.0; #ambiant pressure [Pa]

mu_i = 3.84e-5;#initial He viscosity [Pa*s]
m_i = 0.173;# intial He Mass in vessel [kg]
rho_i = m_i/V;# intial helium density [kg/m^3]
p_i = 1.724e6;1# initial vessel pressure [Pa]

Ti = 773.0; # initial helium temp [K]
To = 273.0; # absolute temp [K]
Twi = 773.0;# initial wall temp [K]

eps = 4.5e-5;# pipe equivalant roughness for a new steel pipe [m]
L = 0.0127; # length of pipe [m]
d = 0.0127;# diameter of hot duct [m]
# d2 = 1.5;
# d3 = 2.29;
Ac = math.pi/4*(d**2);# flow area of duct [m^2]
Pw =math.pi*(d);# Wetted perimeter of duct [m]
D = 4*Ac/Pw;#hydrolic diameter [m]
K_c = 0.5;# Loss coefficent for suden contraction [unitless]
f_g = 0.002;#friction factor guess []
alpha = 1;# kentic energy coefficents []

tol=1e-6;#tolerance
t_max = 10.0;# duration of blowdown transient [s]
delta_t =0.001;#time step [s]
n = int(math.floor(t_max/delta_t));# number of time steps
err = 1;
err1 = 1;



h = 2;# heat transfer coefficent [W/m^2/K]
#h = 10;
#h=40;



T = Ti;
Tw = Twi;
p = p_i;
m = m_i;
rho = rho_i;
mu = mu_i;
t=0 ;
rho_avg=0;
mu_avg = 0;
v_g = 0 ;
Re_g = 0 ;
f1 = 0 ;
f = 0 ;
v= 0;
m_dot=0 ;
Re=0 ;

#User defined
h = 2;

err = 1;
err1 = 1;

m_dot_array= [0];
t_array = [0];
T_array = [Ti];
Tw_array = [Twi];
p_array = [p_i];
v_g_array = [0];
v_array = [0];

for i in range(1,n):
# T = [];
# Tw = [];
# p = [];
# m = [];
# rho = [];
# mu = [];
# t= ;
# rho_avg=0;
# mu_avg = 0;
# v_g =  ;
# Re_g =  ;
# f1 =  ;
# f =  ;
# v= ;

# Re= ;


#User defined

    
#    print "t=  ", t 
    t = ((i-1)*delta_t);
    rho_avg = (rho + rho_amb)/2;
    mu_avg = (mu  + mu_amb)/2;
  
    while abs(err) >= tol:
        v_g =  math.sqrt((p  - p_amb)/(rho_avg*(alpha / 2.0  + 0.5*(K_c + f_g*L/D))));
        if v_g  >= a: 
            Re_g = (rho_avg * a * D/mu_avg);
        else :
            Re_g = (rho_avg * v_g  * D/mu_avg);
        #end
        if Re_g  >= 4000:
            f1 = ( 0.184*Re_g  **(-0.2));
            while abs(err1) >=tol:
                f = ( (1/(-2*math.log(eps/D)/3.7 + 2.51/(Re_g *math.sqrt(f1)))) **2);
                err1 = f1  - f ;
                f1  = f ;

            #end
        else:
               f  = 64.0/Re_g ;
        #end
                
        #print "t=  ", t 
        v = (math.sqrt((p - p_amb)/(rho_avg *(alpha/2.0 + 0.5*(K_c + f *L/D)))));
        err = abs(v_g  - v) ;
        f_g = f ;

        
    #end
    if v  >= a:
        m_dot = ( rho_avg * Ac * a);
        Re = ( rho_avg * a * D / mu_avg);
        v = a;
    else :
        m_dot = ( rho_avg * Ac * v);
        Re = ( rho_avg * v * D / mu_avg);
    #end
    T = ( delta_t * h * As / (c_v * m )*Tw  + (1 - delta_t*(h * As + m_dot *R)/(c_v * m ))*T );
    Tw = (  delta_t * 0.066 * Po *((0)*delta_t + (delta_t **2) **(-0.2))/(rho_w * V_w *c_w) + delta_t * h * As/(rho_w * V_w *c_w)*T  + ( 1 - delta_t * h * As / (rho_w * V_w *c_w))*Tw);
    t = ( t  + delta_t);
    m = ( m  - m_dot  * delta_t);
    mu = ( 1.865e-5 * (T/To) **(0.7));
    rho = ( rho  - m_dot  * delta_t/V);
    p = (rho  * R * T);
    
    #store wanted values
    t_array.append(t);
    m_dot_array.append(m_dot);
    T_array.append(T);
    Tw_array.append(Tw);
    p_array.append(p);
    v_g_array.append(v_g);
    v_array.append(v);
    
    if p  <= p_amb:
        break
    #end
    err = err + 1;
    err1 += 1;



#end		
			
			

#h2 = 10;
#h3 = 40;

#(t,T,Tw,m_dot1) = mass_transfer(h)
#(T2,Tw2) = mass_transfer(h2)
#(T3,Tw3) = mass_transfer(h3)

# print "T =  ", T 
   
#print "Re_g =  ", Re_g 
#print "t=  ", t
# print "Tw =  ", Tw
# for i in m_dot_array:
#     print "m_dot =  ", i
# for i in t_array:
#     print "t =  ", i
# for i in T_array:
#     print "T =  ", i
# for i in Tw_array:
#     print "tw =  ", i
# print "Length of m_dot array = ", len(m_dot_array)
# print "Length of t array = ", len(t_array)
# print "Length of T array = ", len(T_array)
# print "Length of Tw array = ", len(Tw_array)
# print "Length of p array = ", len(p_array)
# print "Length of v_g array = ", len(v_g_array)
#print "t[333] = ", t_array[333]


#write data to csv file
mdotData = [t_array,m_dot_array];
vgData = [t_array,v_g_array];
vData = [t_array,v_array];
pData = [t_array,p_array];

with open('mdot.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(mdotData)

csvFile.close()

with open('vg.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(vgData)
    
csvFile.close()

with open('p.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(pData)
    
csvFile.close()

with open('v.csv', 'w') as csvFile:
    writer = csv.writer(csvFile)
    writer.writerows(vData)
    
csvFile.close()

# print "rho_avg =  ", rho_avg

print "max of v_g_array =",max(v_g_array)
print "max of v_array =",max(v_array)
# print "max of t_array =",max(t_array)
print "a = ", a
			
			
			