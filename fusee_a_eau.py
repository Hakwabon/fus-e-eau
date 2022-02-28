import matplotlib.pyplot as plt
from numpy import *
import math
from fractions import Fraction
import numpy as np
pi = math.pi

V_eau = 1/2000
P_int = 10*101325
t = 0

global R,rho_eau,T_air,V_fusee,dt,S_ejection,P_ext,m_vide,rho_air_0,M_molaire,coef_adi,g,coef_f,n,rho_air_i

#constante :
R = 8.314
T_air = 300
P_ext = 101325
S_ejection = (0.0125**2)*pi
m_vide = 0.1
rho_eau = 1000
V_fusee = 1/1000 #int(input("volume de la fusee en L : "))
dt = 0.00001
rho_air_0 = 1.2
M_molaire = 28.965/1000
coef_adi = 1.4
g = 9.806
coef_f = 1/2 * 1.292 * 0.41 * pi*(0.07*0.06+ 0.06*2*0.3)		# surface total exposé
n = 1/2000*5*101325/(R*T_air)
rho_air_i = P_int*M_molaire/(R*T_air)

def vitesse_eau(V_eau,P_int):
	
	ve = sqrt(2*(P_int-P_ext)/rho_eau)
	
	q_eau = ve*S_ejection
	
	V_eau = V_eau - q_eau*dt
	
	V_air = V_fusee - V_eau
	
	P_int = n*R*T_air/V_air
	
	return ve, V_eau, P_int

def vitesse_air(P_int,t): 
	
	ve = np.sqrt(T_air*R/M_molaire*2*coef_adi/(coef_adi-1)*(1-(P_ext/P_int)**((coef_adi-1)/coef_adi)))
	
	q_air = ve*S_ejection
	
	rho_air = V_fusee*rho_air_i/(V_fusee+q_air*t)
	
	P_int = rho_air*R*T_air/M_molaire
	
	return ve,rho_air,P_int


def acc_1(V_eau,P_int,v):
	
	liste = vitesse_eau(V_eau,P_int)
	ve = liste[0]
	V = liste[1]
	P_int = liste[2]
	
	a = - g + ve**2*rho_eau*S_ejection/(m_vide+V*rho_eau) - coef_f*v**2/(m_vide+V*rho_eau) 
	
	return a, V, P_int

def acc_2(P_int,v,t):
	
	liste = vitesse_air(P_int,t)
	ve = liste[0]
	rho_air = liste[1]
	P_int = liste[2]
	
	a = -g - coef_f*abs(v)*v/(m_vide+V_fusee*rho_air) + ve**2*rho_air*S_ejection/(m_vide+V_fusee*rho_air)
	
	return a, rho_air, P_int

def acc_3(v):
	
	a = -g - coef_f*abs(v)*v
	
	return a



def PFD_1_euler(V_eau,P_int):
	
	vitesse = []
	hauteur = []
	temps = []
	
	t = 0
	h= 0
	v = 0
	
	while V_eau > 0:
		
		li = acc_1(V_eau,P_int,v)
		
		a = li[0]
		V_eau = li[1]
		
		v += a*dt 
		h += v*dt 
		t += dt
		
		vitesse.append(v)
		hauteur.append(h)
		temps.append(t)
	
	P_int = acc_1(V_eau,P_int,v)[2]
	
	return temps,vitesse,hauteur,P_int

def PFD_2_euler(V_eau,P_int):
	
	li = PFD_1_euler(V_eau,P_int)
	
	vitesse = li[1]
	hauteur = li[2]
	temps = li[0]
	
	v = vitesse[len(vitesse)-1]
	h = hauteur[len(vitesse)-1]
	t = temps[len(vitesse)-1]
	
	while P_int > P_ext*1.4 :
		
		li = acc_2(P_int,v,t)
		
		a = li[0]
		P_int = li[2]
		
		v += a*dt 
		h += v*dt 
		t += dt
		
		vitesse.append(v)
		hauteur.append(h)
		temps.append(t)
	
	return temps, vitesse, hauteur

def PFD_3_euler(V_eau,P_int):
	
	li = PFD_2_euler(V_eau,P_int)
	
	vitesse = li[1]
	hauteur = li[2]
	temps = li[0]
	
	v = vitesse[len(vitesse)-1]
	h = hauteur[len(vitesse)-1]
	t = temps[len(vitesse)-1]
	
	while h > 0 :
		
		a = acc_3(v)
		v += a*dt
		h += v*dt
		t += dt
		
		vitesse.append(v)
		hauteur.append(h)
		temps.append(t)
		
	return temps, vitesse, hauteur



def PFD_1_RK4(V_eau,P_int):	
	"oui"

"""
def test(P_int):
	t = 0
	V = []
	P = []
	v = []
	T = []
	
	vitesse = 0
	
	while P_int > P_ext*2:
		
		
		P.append(P_int)
		liste = acc_2(P_int,vitesse,t)
		rho_air = liste[1]
		P_int = liste[2]
		print(P_int)
		vitesse += liste[0]*dt
		v.append(vitesse)
		V.append(rho_air)
		t+= dt
		
		T.append(t)
	
	li = [T,v,V,P]
	return li


def test_air(P_int):
	
	t = 0
	V = []
	P = []
	v = []
	T = []
	
	while P_int > P_ext*1.4:
		
		P.append(P_int)
		liste = vitesse_air(P_int,t)
		ve = liste[0]
		rho_air = liste[1]
		P_int = liste[2]
		
		V.append(rho_air)
		v.append(ve)
		
		t += dt
		
		T.append(t)
		
	return T,v,P,V
"""	


V = Fraction(input("volume d'eau en litre : "))/1000
P = float(input("Pression en bar : "))*P_ext
li = PFD_3_euler(V,P)

v_max = round(max(li[1]),2)
h_max = round(max(li[2]),1)
v_final = round(li[1][len(li[2])-1],2)
Ec_f = round(1/2*m_vide*v_final**2,2)


print("la vitesse maximale est ",v_max,"m/s")
print("la hauteur maximale est ",h_max,"m")
print("la vitesse finale est",v_final,"pour une energie cinetique de",Ec_f,"J")

plt.subplot(2,1,1)
plt.title('vitesse')
plt.plot(li[0],li[1])
plt.xlabel('temps')
plt.ylabel('vitesse')
plt.grid()

plt.subplot(2,1,2)
plt.title("altitude" )
plt.plot(li[0],li[2])
plt.xlabel('temps')
plt.ylabel('altitude')
plt.grid()

plt.show()
plt.close()
