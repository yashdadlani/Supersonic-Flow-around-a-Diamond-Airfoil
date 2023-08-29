import numpy as np
from math import *
from scipy.optimize import newton
import matplotlib.pyplot as plt
k = 1.4 ; a = k+1 ; b = k-1
M_infi = 3; p_infi = 1 ; n = 100001 ; R = 287
comp_angle = float(input("Enter a comp angle: "))
expan_angle = float(input("Enter a expan angle: "))
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
nu_M2_upper = prandtl_meyer(M_infi) + expan_angle

f = lambda M: (180/np.pi)* (sqrt(a/b) *np.arctan(sqrt((b/a)*((M**2)-1)))  - np.arctan(sqrt((M**2)-1))) - nu_M2_upper
M2_upper = newton(f,3)
print("M2_upper is : " + str(M2_upper))
p2_upper = p_infi*((isen_stag_PR(M_infi))/isen_stag_PR(M2_upper))
print("p2_upper is :" + str(p2_upper))
nu_M3_upper = nu_M2_upper + 5
f = lambda M: (180/np.pi)* (sqrt(a/b) *np.arctan(sqrt((b/a)*((M**2)-1)))  - np.arctan(sqrt((M**2)-1))) - nu_M3_upper
M3_upper = newton(f,3.2)
print("M3_upper is : " + str(M3_upper))
p3_upper = p2_upper*((isen_stag_PR(M2_upper))/isen_stag_PR(M3_upper))
print("p3_upper is :" + str(p3_upper))
M = M_infi
beta_array = beta(M)                        #np.empty(1001,float) # this is better
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
Mn1_lower = M_infi*sin((np.pi/180)*beta_array[d])
Mn2_lower = down_mach(Mn1_lower)
M2_lower = Mn2_lower / (sin((np.pi/180)*(beta_array[d] - theta_array[d])))
print("M2_lower is : " + str(M2_lower))
p2_lower = p_infi*pressure_ratio(Mn1_lower)
print("p2_lower is :" + str(p2_lower))
nu_M3_lower = prandtl_meyer(M2_lower) + 5
f = lambda M: (180/np.pi)* (sqrt(a/b) *np.arctan(sqrt((b/a)*((M**2)-1)))  - np.arctan(sqrt((M**2)-1))) - nu_M3_lower
M3_lower = newton(f,1.5)
print("M3_lower is : " + str(M3_lower))
p3_lower = p2_lower*((isen_stag_PR(M2_lower))/isen_stag_PR(M3_lower))
print("p3_lower is :" + str(p3_lower))

phi = float(input("Assume flow angle: "))
comp_angle_trailing_edge = phi + comp_angle
expan_angle_trailing_edge = expan_angle + phi
M = M3_upper
beta_array1 = beta(M)                        #np.empty(1001,float) # this is better
theta_array1 = np.empty(n,float)
for i in range(0,n):
    theta_array1[i] = np.arctan((2/tan(beta_array1[i]))*((M*sin(beta_array1[i]))**2 - 1) / ( (M**2)*(k+ cos(2*beta_array1[i])) + 2 ))   # M is M3_upper
for i in range(0,n):
    beta_array1[i] = beta_array1[i] * (180/np.pi)
    theta_array1[i] = theta_array1[i] * (180/np.pi)
for i in range(0,n):
    if abs(theta_array1[i] - comp_angle_trailing_edge) <= 0.001:
        print(i, theta_array1[i], beta_array1[i])
d1 = int(input("Enter a beta array index of your choice: "))
Mn3_upper = M3_upper*sin((np.pi/180)*beta_array[d1])
p4_upper = p3_upper*pressure_ratio(Mn3_upper)
print("p4_upper is :" + str(p4_upper))
nu_M4_lower = prandtl_meyer(M3_lower) + expan_angle_trailing_edge
f = lambda M: (180/np.pi)* (sqrt(a/b) *np.arctan(sqrt((b/a)*((M**2)-1)))  - np.arctan(sqrt((M**2)-1))) - nu_M4_lower
M4_lower = newton(f,1.2)
print("M4_lower is : " + str(M4_lower))
p4_lower = p3_lower * ((isen_stag_PR(M3_lower))/isen_stag_PR(M4_lower))
print("p4_lower is : " + str(p4_lower))

p_upper = [p2_upper, p2_upper,p3_upper,p3_upper]
cp_upper = np.empty(4,float)
x = [0,1/2,1/2, 1]
p_lower = [p2_lower, p2_lower,p3_lower,p3_lower]
cp_lower = np.empty(4,float)
for i in range(0,4):
    cp_upper[i] = (2/(9*1.4*287))*(p_upper[i] - 1)*1000
for i in range(0,4):
    cp_lower[i] = (2/(9*1.4*287))*(p_lower[i] - 1)*1000
print(cp_upper)
print(cp_lower)
A, = plt.plot(x,cp_upper, color = 'blue', label = 'Upper')
B, = plt.plot(x,cp_lower, color = 'red', label = 'Lower')
plt.xlabel("Chord Length ")
plt.ylabel("Cp *(1000)")
plt.legend(['Upper', "Lower"])
plt.title("Cp variation wrt x for alpha = 5 deg")
plt.grid()
plt.show()


