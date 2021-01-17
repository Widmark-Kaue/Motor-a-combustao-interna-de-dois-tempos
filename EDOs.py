#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 16:25:45 2020

@author: luana and Widmark
"""
import numpy as np
#%%
############################Métodos de Resolução de EDOs####################################

"""
Implementação do método de Euler    
"""
def odeEuler(f, y0, t0, NUMBER_OF_STEPS = 100, h=0.01):
    """
    Resolve EDOs de primeira ordem pelo método de Euler.
    
    Parameters
    ----------
    f : Function
        EDO na forma padrão "y' = f(t, y)".
    y0 : Float
        Valor da função para o P.V.I.
    t0 : Float
        Ponto inicial em que a função é conhecida.
    NUMBER_OF_STEPS : Int, optional
        Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo em que o vetor t é atualizado. The default is 0.01.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    y : Array
        Vetor com as imagens correspondentes de t.

    """
    y = np.zeros(NUMBER_OF_STEPS, dtype = np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32)
    
    #Passandos as condições de contorno
    y[0] = y0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS-1):
        k1 = f(t[n], y[n])
        y[n+1] = y[n] + k1*h
        t[n+1] = t[n]+h
        
    return t, y

def odeEulerSys(f, r0, t0, NUMBER_OF_STEPS = 100, h=0.01):
    """
    Resolve sistema de EDOs de primeira ordem pelo método de Euler.

    Parameters
    ----------
    f :  Function
        EDOs escritas na forma padrão "r' =  f(t,r)", em que r é um vetor com "n" EDOs do sistema.
    r0 : Array
        Condições de contorno para as "n" equações.
    t0 : Float
        Ponto inicial em que as funções são conhecidas.
    NUMBER_OF_STEPS : Int, optional
         Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo em que o vetor t é atualizado. The default is 0.01.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    r : Array
        Vetor com as imagens correspondentes de t para cada função.

    """
    NUMBER_OF_EQUATIONS =  len(r0) #verificando o número de equações diferenciais
    
   #matriz com as funções x, y, z ...
    r = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_EQUATIONS], dtype = np.float32) 
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #vetor com valores de tempo
    
    r[0] = r0 #passando as condições de contorno para a linha zero da matriz de funções
    t[0] = t0 #passando a condição de contorno para o vetor de tempo

    for n in range(0, NUMBER_OF_STEPS-1):
       t[n+1] = t[n]+h #Atualiza o vetor de tempo de acordo com o passo
       K1 = f(t[n], r[n]) #retorna k1x, k1y, k1z,... dentro de um vetor K
       r[n+1] = r[n] + K1*h #armazena os novos valores das funções na próxima linha da matriz
       
    return t, r

"""
Implementação do método de Euler melhorado, ou método de Heun    
"""

def odeHeun(f, y0, t0, NUMBER_OF_STEPS = 100, h=0.01):
    """
    Resolve EDOs de primeira ordem pelo método de Euler melhorado ou método de Heun.
    
    Parameters
    ----------
    f : Function
        EDO na forma padrão "y' = f(t, y)".
    y0 : Float
        Valor da função para o P.V.I.
    t0 : Float
        Ponto inicial em que a função é conhecida.
    NUMBER_OF_STEPS : Int, optional
        Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo em que o vetor t é atualizado. The default is 0.01.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    y : Array
        Vetor com as imagens correspondentes de t.

    """
    y = np.zeros(NUMBER_OF_STEPS, dtype = np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32)
    
    #Passandos as condições de contorno
    y[0] = y0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS-1):
        t[n+1] = t[n]+h
        k1 = f(t[n], y[n])
        k2 = f(t[n+1], y[n] + k1*h)
        y[n+1] = y[n] + 0.5*(k1+k2)*h
        
    return t, y


def odeHeunSys(f, r0, t0, NUMBER_OF_STEPS = 100, h=0.01):
    """
    Resolve sistema de EDOs de primeira ordem pelo método de Euler melhorado ou método de Heun.

    Parameters
    ----------
    f :  Function
        EDOs escritas na forma padrão "r' =  f(t,r)", em que r é um vetor com "n" EDOs do sistema.
    r0 : Array
        Condições de contorno para as "n" equações.
    t0 : Float
        Ponto inicial em que as funções são conhecidas.
    NUMBER_OF_STEPS : Int, optional
         Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo em que o vetor t é atualizado. The default is 0.01.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    r : Array
        Vetor com as imagens correspondentes de t para cada função.

    """
    NUMBER_OF_EQUATIONS =  len(r0) #verificando o número de equações diferenciais
    
    #matriz com as funções x, y, z ...
    r = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_EQUATIONS], dtype = np.float32) 
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #vetor com valores de tempo
    

    
    r[0] = r0 #passando as condições de contorno para a linha zero da matriz de funções
    t[0] = t0 #passando a condição de contorno para o vetor de tempo

    for n in range(0, NUMBER_OF_STEPS-1):
        t[n+1] = t[n]+h #Atualiza o vetor de tempo de acordo com o passo
        K1 = f(t[n], r[n]) #retorna k1x, k1y, k1z,... dentro de um vetor K1
        K2 = f(t[n+1], r[n] + K1*h) #retorna k2x, k2y, k2z, ... dentro de vetor K2
        r[n+1] = r[n] +0.5*(K1 + K2)*h #armazena os novos valores das funções na próxima linha da matriz
        
    return t, r

"""
Implementação do método de Runge-Kutta    
"""
def odeRunge_Kutta(f, y0, t0, NUMBER_OF_STEPS = 100, h=0.01):
    """
    Resolve EDOs de primeira ordem pelo método de Runge-Kutta de quarta ordem ou método de clássico de Runge-Kutta.
    
    Parameters
    ----------
    f : Function
        EDO na forma padrão "y' = f(t, y)".
    y0 : Float
        Valor da função para o P.V.I.
    t0 : Float
        Ponto inicial em que a função é conhecida.
    NUMBER_OF_STEPS : Int, optional
        Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo em que o vetor t é atualizado. The default is 0.01.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    y : Array
        Vetor com as imagens correspondentes de t.

    """
    y = np.zeros(NUMBER_OF_STEPS, dtype = np.float32)
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32)
    
    #Passandos as condições de contorno
    y[0] = y0
    t[0] = t0

    for n in range(0, NUMBER_OF_STEPS-1):
        t[n+1] = t[n] + h
        k1 = f(t[n], y[n])
        k2 = f(t[n] + 0.5*h , y[n] + 0.5*k1*h)
        k3 = f(t[n] + 0.5*h , y[n] + 0.5*k2*h)
        k4 = f(t[n+1], y[n] + k3*h)
        y[n+1] = y[n] + (k1 + 2*k2 + 2*k3 + k4)*(h/6)
        
    return t, y

def odeRunge_KuttaSys(f, r0, t0, NUMBER_OF_STEPS = 100, h=0.01):
    """
    Resolve sistema de EDOs de primeira ordem pelo método de Runge-Kutta de quarta ordem.
    Parameters
    ----------
    f :  Function
        EDOs escritas na forma padrão "r' =  f(t,r)", em que r é um vetor com "n" EDOs do sistema.
    r0 : Array
        Condições de contorno para as "n" equações.
    t0 : Float
        Ponto inicial em que as funções são conhecidas.
    NUMBER_OF_STEPS : Int, optional
         Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo em que o vetor t é atualizado. The default is 0.01.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    r : Array
        Vetor com as imagens correspondentes de t para cada função.

    """
    NUMBER_OF_EQUATIONS =  len(r0) #verificando o número de equações diferenciais
    
    #matriz com as funções x, y, z ...
    r = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_EQUATIONS], dtype = np.float32) 
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #vetor com valores de tempo
    

    
    r[0] = r0 #passando as condições de contorno para a linha zero da matriz de funções
    t[0] = t0 #passando a condição de contorno para o vetor de tempo

    for n in range(0, NUMBER_OF_STEPS-1):
        t[n+1] = t[n]+h #Atualiza o vetor de tempo de acordo com o passo
        K1 = f(t[n], r[n]) #retorna k1x, k1y, k1z,... dentro de um vetor K1
        K2 = f(t[n] + 0.5*h, r[n] + 0.5*K1*h) #retorna k2x, k2y, k2z, ... dentro de vetor K2
        K3 = f(t[n] + 0.5*h, r[n] + 0.5*K2*h) #retorna k3x, k3y, k3z,... dentro de um vetor K3
        K4 = f(t[n+1], r[n] + K3*h) #retorna k4x, k4y, k4z, ... dentro de vetor K4
        r[n+1] = r[n] + (K1 + 2*K2 + 2*K3 + K4)*(h/6) #armazena os novos valores das funções na próxima linha da matriz
        
    return t, r

"""
Implementação do método de Runge-Kutta-Fehlberg    
"""
def odeRunge_Kutta_Fehlberg(f, y0, t0, NUMBER_OF_STEPS = 100, h=0.01, alpha = 0.9, e = 1E-3):
    """
     Resolve EDOs de primeira ordem pelo método de Runge-Kutta-Fehlberg.
     Método adaptativo que consiste em máximizar h de acordo com o comportamento da solução
     mantendo o erro dentro de uma tolerância.
    
    Parameters
    ----------
    f : Function
        EDO na forma padrão "y' = f(t, y)".
    y0 : Float
        Valor da função para o P.V.I.
    t0 : Float
        Ponto inicial em que a função é conhecida.
    NUMBER_OF_STEPS : Int, optional
        Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo inicial em que o vetor t é atualizado. The default is 0.01.
    alpha: Float optional
        Magnitude na qual q deve ser menor que uma dada expressão para que h seja adpatado.
    e: Float, optional
        Tolerância máxima de erro admissível para um determinado passo.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    y : Array
        Vetor com as imagens correspondentes de t.

    """
    y = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #vetor para a função de 4º ordem
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #passo de tempo
    H = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #variável auxiliar para monitorar variável h
    
    #Passandos as condições de contorno
    y[0] = y0
    t[0] = t0


    for n in range(0, NUMBER_OF_STEPS - 1):
        H[n] = h
        
        while (True):
            #Cálculos dos parâmetros Ks
            k1 = f(t[n], y[n])
            k2 = f(t[n] + h/4, y[n] + (k1/4)*h)
            k3 = f(t[n] + (3/8)*h , y[n] + h*((3*k1 + 9*k2)/32))
            k4 = f(t[n] + (12/13)*h, y[n] + h*((1932*k1 - 7200*k2 + 7296*k3)/2197))
            k5 = f(t[n] + h, y[n] + h*((439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4))
            k6 = f(t[n] + 0.5*h, y[n] + h*(-(8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5))
            
            #Cálculando o valor das funções no tempo n+1
            y5 = y[n] + ((16/135)*k1 + (6656/12825)*k3 + (28561/56430)*k4 - (9/50)*k5 + (2/55)*k6)*h #Runge-Kutta de 5º ordem
        
            y4 = y[n] + ((25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - (1/5)*k5)*h #Runge-Kutta de 4º ordem
            
              
            #Verificando se o erro está dentro da tolerância
            Num = e*h
               
            Den = abs(y4 - y5)
            
            if (Den != 0):
                q = alpha*((Num/Den)**(1/4)) #Problema: O valor não está sendo atualizado corretamente
                
                if (q < 1):
                    h = q*h  #diminui o passo h para uma fração de h que esteja dentro da tolerância.

                else:
                    break
            else:
             break
        
        y[n+1] = y5 #Atualiza a função.
        
        t[n+1] = t[n] + h #Atualização do vetor de tempo.
        
        if (Den !=0): #verifica se o q foi calculado
            h = q*h #aumenta o passso h para o próximo passo
        
        H[n] = h
        
    return t, y, H

def odeRunge_Kutta_FehlbergSys(f, r0, t0, NUMBER_OF_STEPS = 100, h=0.01, alpha = 0.9, e = 1E-3):
    """
     Resolve Sistema de EDOs de primeira ordem pelo método de Runge-Kutta-Fehlberg.
     Método adaptativo que consiste em máximizar h de acordo com o comportamento da solução
     mantendo o erro dentro de uma tolerância.
    
    Parameters
    ----------
    f : Function
        EDO na forma padrão "y' = f(t, y)".
    r0 : Array, lista ou tupla
        Condições de contorno para as "n" equações.
    t0 : Float
        Ponto inicial em que a função é conhecida.
    NUMBER_OF_STEPS : Int, optional
        Número de passos executados pelo método. The default is 100.
    h : Float, optional
        Passo inicial em que o vetor t é atualizado. The default is 0.01.
    alpha: Float optional
        Magnitude na qual q deve ser menor que uma dada expressão para que h seja adpatado.
    e: Float, optional
        Tolerância máxima de erro admissível para um determinado passo.

    Returns
    -------
    t : Array
        Vetor com valores da abcissa.
    y : Array
        Vetor com as imagens correspondentes de t.

    """
    NUMBER_OF_EQUATIONS =  len(r0) #verificando o número de equações diferenciais
    
    
    #matriz com as funções x, y, z ...
    r = np.zeros([NUMBER_OF_STEPS, NUMBER_OF_EQUATIONS], dtype = np.float32) 
    t = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #vetor com valores de tempo
    H = np.zeros(NUMBER_OF_STEPS, dtype = np.float32) #variável auxiliar para monitorar variável h
    
    #Passandos as condições de contorno
    r[0] = r0
    t[0] = t0


    for n in range(0, NUMBER_OF_STEPS - 1):
     
        H[n] = h
        
        
        while (True):
            
            #Cálculos dos parâmetros Ks
            K1 = f(t[n], r[n])
            K2 = f(t[n] + h/4, r[n] + (K1/4)*h)
            K3 = f(t[n] + (3/8)*h , r[n] + h*((3*K1 + 9*K2)/32))
            K4 = f(t[n] + (12/13)*h, r[n] + h*((1932*K1 - 7200*K2 + 7296*K3)/2197))
            K5 = f(t[n] + h, r[n] + h*((439/216)*K1 - 8*K2 + (3680/513)*K3 - (845/4104)*K4))
            K6 = f(t[n] + 0.5*h, r[n] + h*(-(8/27)*K1 + 2*K2 - (3544/2565)*K3 + (1859/4104)*K4 - (11/40)*K5))
            
            #Cálculando o valor das funções no tempo n+1
            r5 = r[n] + ((16/135)*K1 + (6656/12825)*K3 + (28561/56430)*K4 - (9/50)*K5 + (2/55)*K6)*h #Runge-Kutta de 5º ordem
        
            r4 = r[n] + ((25/216)*K1 + (1408/2565)*K3 + (2197/4104)*K4 - (1/5)*K5)*h #Runge-Kutta de 4º ordem
            

              
            #Verificando se o erro está dentro da tolerância
            Num = e*h
               
            Den = abs(r4 - r5)
            
# =============================================================================
#  Para o caso em que algum dos valores de Den são igual a zero, esse valor é substituido
#  por um valor específico para que quando q seja calculado o valor dê igual 459.
# =============================================================================
            
            for i in range (len (Den)): #verifica se existe algum valor nulo no vetor Den.
               if (Den[i] == 0):
                   Den[i] = (e*Num)/((459/alpha)**(4)) #Possibilita o cálculo do q.
                
            q = alpha*((Num/Den)**(1/4)) #calcula o q
                
            if any(i < 1 for i in q): #verifica se algum valor do vetor q é menor que 1.
                
                h = min(q)*h  #utiliza o menor valor de q para diminuir o passo h para uma fração de h que esteja dentro da tolerância.
                
            else:
                break
            
        
        r[n+1] = r5 #Atualiza a função.
        
        t[n+1] = t[n] + h #Atualização do vetor de tempo.
        
        if not all (i == 459 for i in q): #verifica se o passo h pode ser aumentado
            h = min(q)*h #aumenta o passso h para o próximo passo
        
    return t, r, H

#%%
################################################################################################




