import numpy as np
from math import *
import matplotlib.pyplot as plt
k = 1.4 ; a = k+1 ; b = k-1
M_infi = 3; p_infi = 1 ; n = 100001 ; 
comp_angle = float(input("Enter a comp angle: "))
def down_mach(M):
    return sqrt((2 + b*(M**2)) / (2*k*(M**2) - b))
def pressure_ratio(M):
    return 1 + ((2*k / (a)) *((M**2) - 1))
def mu(M):
    x = 1/M
    y = np.arcsin(x)
    return y
def beta(M):
    return np.linspace(mu(M),np.pi/2,n )
def isen_stag_PR(M):
    return ((1 + 0.5*b*(M**2))**(k/b))
def prandtl_meyer(M):
    return (180/np.pi)* (sqrt(a/b) *np.arctan(sqrt((b/a)*((M**2)-1)))  - np.arctan(sqrt((M**2)-1)))
M = M_infi
beta_array = beta(M)    #np.empty(1001,float) # this is better
theta_array = np.empty(n,float)
for i in range(0,n):
    theta_array[i] = np.arctan((2/tan(beta_array[i]))*((M*sin(beta_array[i]))**2 - 1) / ( (M**2)*(k+ cos(2*beta_array[i])) + 2 ))   # M is M_infi
for i in range(0,n):
    beta_array[i] = beta_array[i] * (180/np.pi)
    theta_array[i] = theta_array[i] * (180/np.pi)
for i in range(0,n):
    if abs(theta_array[i] - comp_angle) <= 0.001:
        print(i, theta_array[i], beta_array[i])

Mn1_lower = M_infi*sin((np.pi/180)*beta_array[5193])
Mn1_upper = 1
Mn2_lower = down_mach(Mn1_lower)
M2_lower = Mn2_lower / (sin((np.pi/180)*(beta_array[5193] - theta_array[5193])))
M2_upper = 3
M3_upper = 3.27341
M3_lower = 2.99627
# print(M2_lower)
p2_upper = p_infi*pressure_ratio(Mn1_upper)
p2_lower = p_infi*pressure_ratio(Mn1_lower)
p3_upper = p2_upper*((isen_stag_PR(M2_upper))/isen_stag_PR(M3_upper))
p3_lower = p2_lower*((isen_stag_PR(M2_lower))/isen_stag_PR(M3_lower))


p_upper = [p2_upper, p2_upper,p3_upper,p3_upper]
cp_upper = np.empty(4,float)
x = [0,1/2,1/2, 1]
p_lower = [p2_lower, p2_lower,p3_lower,p3_lower]
cp_lower = np.empty(4,float)
for i in range(0,4):
    cp_upper[i] = (2/(9*1.4*287))*(p_upper[i] - 1)*1000
for i in range(0,4):
    cp_lower[i] = (2/(9*1.4*287))*(p_lower[i] - 1)*1000
A, = plt.plot(x,cp_upper, color = 'blue', label = 'Upper')
B, = plt.plot(x,cp_lower, color = 'red', label = 'Lower')
plt.xlabel("Chord Length ")
plt.ylabel("Cp *(1000)")
plt.legend(['Upper', "Lower"])
plt.title("Cp variation wrt x for alpha = 2.5 deg")
plt.grid()
plt.show()