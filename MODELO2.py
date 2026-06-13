"""
Modelo 2 del biochip alvéolo--capilar: advección-reacción-difusión

Este script resuelve, mediante diferencias finitas la ecuación:

    W_x = (1/Pe_C)W_xx + Da_C (C - W),

con condiciones de Dirichlet en la entrada y Neumann en la salida.

Salidas gráficas:
    1) Comparación de solución analítica con la numérica.
    2) Gráfica de convergencia a escala log-log

"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded
from pathlib import Path

DIRECTORIO_SCRIPT = Path(__file__).resolve().parent
DIRECTORIO_SALIDA_MODELO_2 = DIRECTORIO_SCRIPT / "figuras_modelo2"

# GUARDAR FIGURAS

salida = DIRECTORIO_SALIDA_MODELO_2

def guardar_figura(fig: plt.Figure, ruta_base: Path) -> None:
    """Guarda una figura en PNG y PDF con resolución alta."""

    ruta_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(ruta_base.with_suffix(".pdf"), bbox_inches="tight")


plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 16,      # Etiquetas de ejes (ej: "$\hat{x}$")
    'axes.titlesize': 18,      # Títulos
    'legend.fontsize': 13,     # Leyendas
    'xtick.labelsize': 14,     # ← NÚMEROS en eje X
    'ytick.labelsize': 14,     # ← NÚMEROS en eje Y
})

# ============================================================
# PARÁMETROS ADIMENSIONALES
# ============================================================

PeC = 20.0          
DaC = 3.0          
CA0 = 1.0           # Concentración alveolar fija
Win = 0.15          # Condición de entrada

# ============================================================
# FUNCIÓN SOLUCIÓN ANALÍTICA
# ============================================================

'''
   La ecuación se reescribe como:
    
     (1/PeC)W_xx - W_x - DaC W = -DaC CA0 
     
     W(0) = Win  Dirichlet en la entrada
     W_x(1) = 0  Neumann en la salida

   EDO lineal de segundo orden con coeficientes constantes y término inhomogéneo. 

    - Solución particular: W_p = CA0
    - Solución parte homogénea: las raíces de la ecuación característica de la parte homogénea son 
                                λ_mas = (PeC + sqrt(PeC**2 + 4 * PeC * DaC)) / 2
                                λ_menos = (PeC - sqrt(PeC**2 + 4 * PeC * DaC)) / 2

   La solución imponiendo condiciones de contorno queda:

      W(x) = CA0 + (Win - CA0)(λ_mas e^(λ_mas)e^(λ_menos x) - λ_menos e^(λ_menos)e^(λ_mas x))/ (λ_mas e^(λ_mas) - λ_menos e^(λ_menos))
     
'''

def solucion_analitica(x, PeC, DaC, CA0, Win):

    # Calculamos las raíces del polinomio característico
    sqrt_term = np.sqrt(PeC**2 + 4 * PeC * DaC)
    lambda_mas = (PeC + sqrt_term) / 2
    lambda_menos = (PeC - sqrt_term) / 2

    # Escribimos numerador de la solución
    denom = (
        lambda_mas * np.exp(lambda_mas)
        - lambda_menos * np.exp(lambda_menos)
    )

    # Escribimos denominador de la solución
    term = (
        lambda_mas * np.exp(lambda_mas)
        * np.exp(lambda_menos * x)
        - lambda_menos * np.exp(lambda_menos)
        * np.exp(lambda_mas * x)
    )

    # Escribimos solución
    W = (
        CA0 + (Win - CA0) * term / denom
    )

    return W

# ============================================================
# FUNCIÓN SOLUCIÓN NUMÉRICA
# ============================================================

'''
   Se reescribe la ecuación usando diferencias finitas: 
  
      (1/PeC) (W[i+1] -2W[i] + W[i-1])/dx^2 - (W[i] - W[i-1])/dx - DaC W = -DaC CA0 

   Se reescribe la ecuación como a_menos W[i-1] + a_0 W[i] + a_mas W[i+1] = b donde:
      
       a_menos = 1/PeC + dx/2
       a_0 = -2/PeC - DaC dx^2
       a_mas = 1/PeC - dx/2
       b = - DaC dx^2 CA0
      
   Se quiere resolver el sistema AW = b donde 

       W = (W_1, ... , W_N)^T 
       A es la matriz de diagonal principal a_0, diagonal superior a_mas y diagonal inferior a_menos.
       b es el vector de términos independientes constante para todos los nodos salvo para el primero que añade la condición de contorno en la entrada

'''

def solucion_numerica(N, PeC, DaC, CA0, Win):

    dx = 1.0 / N

    x = np.linspace(0, 1, N + 1)

    # --------------------------------------------------------
    # Coeficientes del esquema
    # --------------------------------------------------------

    a_menos = 1/PeC + dx/2
    a_mas   = 1/PeC - dx/2
    a_0     = -2/PeC - DaC * dx**2

    b = -DaC * dx**2 * (CA0)

    # --------------------------------------------------------
    # Sistema tridiagonal
    # --------------------------------------------------------

    # Número de incógnitas:
    # W_1, W_2, ..., W_N
    M = N

    diag_infer = np.zeros(M)
    diag  = np.zeros(M)
    diag_super = np.zeros(M)

    rhs = np.zeros(M)

    # ---------------------------------------------------------------------------------------------------------------------------------
    # Nodo i = 1 => a_menos W_0 + a_0 W_1 + a_mas W_2 = -DaC CA0 dx^2 => (usando W_0=Win) a_0 W_1 + a_mas W_2 = -DaC CA0 dx^2 - a_menos Win
    # ---------------------------------------------------------------------------------------------------------------------------------

    diag[0] = a_0
    diag_super[0] = a_mas

    rhs[0] = b - a_menos * Win # añadimos condición de contorno al vector de términos independientes al primer nodo

    # --------------------------------------------------------
    # Nodos interiores
    # --------------------------------------------------------

    for i in range(1, M-1):

        diag_infer[i] = a_menos
        diag[i]  = a_0
        diag_super[i] = a_mas

        rhs[i] = b

    # ---------------------------------------------------------------------------------------------------------------------------------------------------
    # Nodo i = N => a_menos W_{N-1} + a_0 W_N + a_mas W_{N+1} = -DaC CA0 dx^2 => (usando W_{N-1}=W_{N+1}) (a_menos + a_mas)W_{N-1} + a_0 W_N = -DaC CA0 dx^2
    # ---------------------------------------------------------------------------------------------------------------------------------------------------

    diag_infer[M-1] = a_menos + a_mas
    diag[M-1]  = a_0

    rhs[M-1] = b

    # ---------------------------------------------------------------------------------------------------------------------------
    # Formato banded para scipy: guarda las diagonales en una matriz 3xN y opera sabiendo que forman las diagonales de una matriz
    # ---------------------------------------------------------------------------------------------------------------------------

    ab = np.zeros((3, M))

    ab[0, 1:] = diag_super[:-1]
    ab[1, :]  = diag
    ab[2, :-1] = diag_infer[1:]

    # Resolver sistema
    W_internal = solve_banded((1,1), ab, rhs)

    # Añadir condición de Dirichlet
    W= np.zeros(N + 1)

    W[0] = Win
    W[1:] = W_internal

    return x, W

# ============================================================
# COMPARACIÓN ANALÍTICA VS NUMÉRICA
# ============================================================

N = 100

x, W_num = solucion_numerica(
    N, PeC, DaC, CA0, Win)

W_exact = solucion_analitica(
    x, PeC, DaC, CA0,Win)


# ============================================================
# GRÁFICA DE PERFILES
# ============================================================

plt.figure(figsize=(7,5))

plt.plot(
    x,
    W_exact,
    label='Solución analítica',
    linewidth=2)

plt.plot(
    x,
    W_num,
    '--',
    label='Solución numérica')

plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{W}(\hat{x})$')

plt.title('Modelo 2: advección–difusión–reacción')

plt.legend()
plt.grid(True)

guardar_figura(plt, salida/"01_Solucion" )

# ============================================================
# ESTUDIO DE CONVERGENCIA
# ============================================================

'''
   Para estudiar la convergencia se hace una gráfica log-log de E(N) vs N:

      Se asume que el error es proporcional a 1/N por usar upwind

         E(N) = C/N^2 -> log_10 (E) = log_10 (C) - 2log_10 (N)  (y = b - x    función del error a escala logarítmica con pendiente esperada -2)
'''

N_vals = np.array([10, 20, 40, 80, 160, 320, 640, 1280, 2560])

errores = []

for N in N_vals:

    x_num, W_num = solucion_numerica(
        N,
        PeC,
        DaC,
        CA0,
        Win
    )
    # Evaluamos solucion analitica en los mismos puntos que la numerica
    W_exact = solucion_analitica(
        x_num,
        PeC,
        DaC,
        CA0,
        Win
    )

    error = np.max(np.abs(W_num - W_exact))
    errores.append(error) #voy metiendo en un vector los errores conforme avanza N

errores = np.array(errores)

# Ajuste lineal en escala log-log para calcular la pendiente
log_N = np.log10(N_vals)
log_E = np.log10(errores)
pendiente, _ = np.polyfit(log_N, log_E, 1) #esta función nos da la pendiente y la b, solo queremos la pendiente


# GRÁFICA DE CONVERGENCIA

plt.figure(figsize=(7,5))

plt.loglog(
    N_vals,
    errores,
    'o-',
    label=f'Error numérico\n(Pendiente log-log ≈ {pendiente:.2f})')

# Línea de referencia teórica con pendiente -2
ref_line = errores[0] * (np.array(N_vals)/N_vals[0])**(-2)

#Añado la función de pendiente -2
plt.loglog(
    N_vals,
    ref_line,
    '--',
    label='Pendiente esperada = -2')

plt.title('Estudio de Convergencia')
plt.xlabel('Número de Nodos ($N$)')
plt.ylabel('Error Máximo $E(N)$')
plt.legend()
plt.grid(True, which="both", ls="--")

plt.tight_layout()
guardar_figura(plt, salida/"02_Convergencia" )



