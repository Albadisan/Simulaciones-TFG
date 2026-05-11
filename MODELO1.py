"""
Modelo 1 del biochip alvéolo--capilar: advección-reacción puro

Este script resuelve, mediante diferencias finitas la ecuación:

    W_x = Da_C (C - W),

con condiciones de Dirichlet en la entrada. 

Salidas gráficas:
    1) Comparación de solución analítica con dos numéricas para N=100 y N=1000.
    2) Gráfica de convergencia a escala log-log

"""
import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# DEFINICIÓN DE PARÁMETROS ADIMENSIONALES
# ==========================================
DaC = 3.0      
CA0 = 1.0       # Concentración alveolar fija
Win = 0.15      # Condición inicial

# ==========================================
# FUNCIÓN SOLUCIÓN ANALÍTICA
# ==========================================
def solucion_analitica(x):
    return CA0 + (Win - CA0) * np.exp(-DaC * x)

# ==========================================
# FUNCIÓN SOLUCIÓN NUMÉRICA
# ==========================================

'''
   Se reescribe la ecuación usando esquema upwind:
    
       (W[i] - W[i-1])/dx = Da_C (CA0 - W[i])

   Se despeja W[i]:
       
       W[i] = (W[i-1] + DaC * dx * CA0) / (1 + DaC * dx)

    Se hace un bucle iterativo sobre i almacenando las soluciones en el vector solución W.
'''

def solucion_numerica(N):
    dx = 1.0 / N  
    
    W = np.zeros(N + 1)
    W[0] = Win    # Condición inicial 
    
    x = np.linspace(0, 1.0, N + 1)
    
    # Bucle simple iterativo para i = 1, ..., N
    for i in range(1, N + 1):
        W[i] = (W[i-1] + DaC * dx * CA0) / (1 + DaC * dx)
        
    return x, W

# ==========================================
# Comparación para mallas N=100 y N=1000
# ==========================================
'''
   Se grafica la solución analítica para la misma malla que la numérica y se saca por pantalla una gráfica con las tres funciones.
'''

x_100, W_100 = solucion_numerica(100)
x_1000, W_1000 = solucion_numerica(1000)

# Malla fina para graficar la solución analítica suave
x_exact = np.linspace(0, 1.0, 1000)
W_exact = solucion_analitica(x_exact)

plt.figure(figsize=(7, 5))

# Gráfica solución exacta
plt.plot(
    x_exact, 
    W_exact, 
    'k-', 
    label='Solución analítica de Ec. 39', 
    linewidth=2)

# Gráfica solución para N=100
plt.plot(
    x_100, 
    W_100, 
    'r--', 
    label='Solución numérica N=100')

# Gráfica solución para N=100
plt.plot(
    x_1000, 
    W_1000, 
    'b:', 
    label='Solución numérica N=1000')

plt.title('Comparación de Soluciones')
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{W}(\hat{x})$')
plt.legend()
plt.grid(True)

# ==========================================
# ESTUDIO DE CONVERGENCIA
# ==========================================

'''
   Para estudiar la convergencia se hace una gráfica log-log de E(N) vs N:

      Se asume que el error es proporcional a 1/N por usar upwind

         E(N) = C/N -> log_10 (E) = log_10 (C) - log_10 (N)  (y = b - x    función del error a escala logarítmica con pendiente esperada -1)
'''

N_vals = np.array([10, 20, 40, 80, 160, 320, 640, 1280, 2560])
errores = []

# Generamos vector de errores
for N in N_vals:
    x_num, W_num = solucion_numerica(N)
    # Evaluamos la analítica en los mismos puntos de la malla numérica
    W_exact = solucion_analitica(x_num)
    
    # Error máximo (Norma L-infinito)
    error = np.max(np.abs(W_num - W_exact))
    errores.append(error) #meto en un vector los errores que van saliendo conforme aumenta N

errores = np.array(errores)

# Ajuste lineal en escala log-log para calcular la pendiente
log_N = np.log10(N_vals)
log_E = np.log10(errores)

# Queremos la función log(E) = m*log(N) + log(C)
pendiente, _ = np.polyfit(log_N, log_E, 1) #esta función nos da la pendiente y la b, solo queremos la pendiente

# GRÁFICA DE CONVERGENCIA

plt.figure(figsize=(7,5))

#dibujamos log-log y comentamos la pendiente obtenida
plt.loglog(
    N_vals, 
    errores, 
    'o-', 
    label=f'Error numérico\n(Pendiente obtenida: {pendiente:.3f})') 

# Línea de referencia teórica con pendiente -1
ref_line = errores[0] * (N_vals[0] / N_vals)

#Añado la funcion de pendiente -1
plt.loglog(
    N_vals, 
    ref_line, 
    'k--', 
    label='Pendiente esperada = -1')

plt.title('Estudio de Convergencia')
plt.xlabel('Número de Nodos ($N$)')
plt.ylabel('Error Máximo $E(N)$')
plt.legend()
plt.grid(True, which="both", ls="--")

plt.tight_layout()
plt.show()
