import numpy as np
from math import *
from scipy.optimize import newton
k = 1.4 ; a = k+1 ; b = k-1 ; n = 100001
# M3_upper = float(input("Enter M3_upper : "))
# M3_lower = float(input("Enter M3_lower : "))
M3_upper = 3.27341
M3_lower = 2.99627
phi = 0.01

comp_angle = 5 + phi
expan_angle = 0 + phi
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

M = M3_upper
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
d = int(input("Enter a beta array index of your choice: "))

Mn3_upper = M3_upper*sin((np.pi/180)*beta_array[d])
p3_upper = 0.6673
p3_lower  = 1.00026
# Mn4_upper = down_mach(Mn3_upper)
nu_M4_lower = phi + prandtl_meyer(M3_lower)
f = lambda M: (180/np.pi)* (sqrt(a/b) *np.arctan(sqrt((b/a)*((M**2)-1)))  - np.arctan(sqrt((M**2)-1))) - nu_M4_lower
M4_lower = newton(f,3)
p4_upper = p3_upper*pressure_ratio(Mn3_upper)
p4_lower = p3_lower*((isen_stag_PR(M3_lower))/isen_stag_PR(M4_lower))



print( "Upper pressure : " + str(p4_upper))
print("Lower Pressure : " + str(p4_lower))
# print("value  of nu_M4_lower is : " + str(nu_M4_lower))


