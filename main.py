# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 16:36:30 2020

@author: Widmark and Luana
"""

import numpy as np
import matplotlib.pylab as plt
import EDOs as ed

#%% Implementando função e setando as condições iniciais

w0 = 50 #frequência angular inicial
theta0 = 0 #crank angle inicial
t0 = 0 #tempo inicial
h = [0.1, 0.05, 0.01, 0.001] #passos tempos utilizados
 
             


def f(t,r):
    w, theta = r #criando os vetores de ângulo e velocidade angular
    
    R = 0.305 #Raio do flywheel em metros
    x0 = 0.0254 #Espaçamento mínimo quando o pistão está em compressão
    I = 3.171 #momento de inércia do flywheel em kgm^2
    C = 0.0113 #momento de inércia em kgm^2
    theta1 = np.pi #crank angle inicial para o processo de expansão em rad
    theta3 = 0 #crank angle inicial para o processo de compressão em rad
    n = 1.3 #Identidficador do processo politrópico
    A = 0.00188 #área de seção transversal do pistão em m^2
    P1 = 0.1E6 #pressão inicial de compressão em Pa
    P3 = 10.3E6 #pressão inicial de expansão em Pa
    
    T1 = (A*R/I)*np.sin(theta)
    T2 = (C/I)*(w**2)
    
    
    while(theta > 2*np.pi):
        theta = theta - 2*np.pi
            
    if (theta >= 0 and theta <= np.pi):
        T3 = P3*((R - R*np.cos(theta3) + x0)/(R - R*np.cos(theta) + x0))**n
       
        
    elif(theta > np.pi and theta <= 2*np.pi):
        T3 = P1*((R - R*np.cos(theta1) + x0)/(R - R*np.cos(theta) + x0))**n
    
    
    return np.array([T3*T1 - T2, w] )



#%%                             Método de Euler



# =============================================================================
# #Comparação dos passos de tempo até 12 segundos
# =============================================================================

for i in range (len(h)):
    t, r = ed.odeEulerSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(13/h[i]), h = h[i])
    plt.title("Método de Euler")
    plt.xlabel("tempo (s)")
    plt.ylabel("Velovidade angular (rad/s)")
    plt.plot(t,r[:,0], label = "h =" + str(h[i]))
    plt.legend()

plt.grid()
plt.show()

# =============================================================================
# # Comparação dos passos de tempo até 0.5 segundos
# =============================================================================
    
for i in range (len(h)):
    t, r = ed.odeEulerSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(0.5/h[i]), h = h[i])
    plt.title("Método de Euler")
    plt.xlabel("tempo (s)")
    plt.ylabel("Velovidade angular (rad/s)")
    plt.plot(t,r[:,0], label = "h =" + str(h[i]))
    plt.legend()

plt.grid()   
plt.show()

# =============================================================================
# Gráficos da velocidade angular com o passo de 0.01
# =============================================================================

t, r = ed.odeEulerSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(13/h[2]), h = h[2])
plt.title("Método de Euler")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0], "r--" ,label = "h =" + str(h[2]))
plt.legend()

plt.grid()
plt.show()


t, r = ed.odeEulerSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(0.5/h[2]), h = h[2])
plt.title("Método de Euler")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0], "r--" ,label = "h =" + str(h[2]))
plt.legend()

plt.grid()
plt.show()


#%%                            Método de Heun

# =============================================================================
# Gráficos da velocidade angular com o passo de 0.01
# =============================================================================

t, r = ed.odeHeunSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(13/h[2]), h = h[2])
plt.title("Metodo de Heun ")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0],"b--", label = "h =" + str(h[2]))
plt.legend()

plt.grid()
plt.show()

t, r = ed.odeHeunSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(0.5/h[2]), h = h[2])
plt.title("Método de Heun")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0],"b--", label = "h =" + str(h[2]))
plt.legend()

plt.grid()   
plt.show()    

    
#%%                        Método de Runge Kutta de 4º ordem

# =============================================================================
# Gráficos da velocidade angular com o passo de 0.01
# =============================================================================

t, r = ed.odeRunge_KuttaSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(12/h[2]), h = h[2])
plt.title("Método de Runge-Kutta de 4º ordem ")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0],"g--", label = "h =" + str(h[2]))
plt.legend()
plt.grid()
plt.show()
    
t, r = ed.odeRunge_KuttaSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(0.5/h[2]), h = h[2])
plt.title("Método de Runge-Kutta de 4º ordem ")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0],"g--", label = "h =" + str(h[2]))
plt.legend()

plt.grid()   
plt.show() 
    
#%%                        Método de Runge-Kutta-Fehlberg

# =============================================================================
# Gráficos da velocidade angular com o passo inicial de 0.1
# =============================================================================


t, r, H = ed.odeRunge_Kutta_FehlbergSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(38/h[0]), h = h[0])
plt.title("Método de Runge-Kutta-Fehlberg ")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0],"y--", label = "h0 =" + str(h[0]))
plt.legend()
plt.grid()
plt.show() 

t, r, H = ed.odeRunge_Kutta_FehlbergSys(f, (w0, theta0), t0, NUMBER_OF_STEPS = int(1.3/h[0]), h = h[0])
plt.title("Método de Runge-Kutta-Fehlberg ")
plt.xlabel("tempo (s)")
plt.ylabel("Velovidade angular (rad/s)")
plt.plot(t,r[:,0],"y--", label = "h0 =" + str(h[0]))
plt.legend()
plt.grid()
plt.show() 
    



#%%

