"""
Modelo 3 del biochip alvéolo--capilar: sistema estacionario lineal acoplado

Este script resuelve el sistema:

    rv C_x = (1/Pe_A)C_xx + Da_A (C - W),
    W_x = (1/Pe_C)W_xx + Da_C (C - W),

con condiciones de Dirichlet en la entrada y Neumann en la salida.

Salidas gráficas:
    1) Comparación de solución analítica con la numérica.

"""

from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csc_matrix, eye, lil_matrix
from typing import Dict, Iterable, List, Tuple
from scipy.sparse.linalg import splu
from pathlib import Path


DIRECTORIO_SCRIPT = Path(__file__).resolve().parent
DIRECTORIO_SALIDA_MODELO_3 = DIRECTORIO_SCRIPT / "figuras_modelo3"

# GUARDAR FIGURAS

salida = DIRECTORIO_SALIDA_MODELO_3

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
# Definición de parámetros y estructuras de datos
# ============================================================


@dataclass(frozen=True)
class CasoModelo:
    """Parámetros adimensionales de un caso de simulación."""

    nombre: str
    rv: float = 1.0
    PeA: float = 20.0
    PeC: float = 20.0
    DaA: float = 3.0
    DaC: float = 3.0
    C_in: float = 1.0
    W_in: float = 0.15


# ============================================================
# Ensamblaje del operador espacial común a los Modelos 3 y 4
# ============================================================

def ensamblar_operador_espacial(N: int, dx: float, caso: CasoModelo) -> csc_matrix:
    """
     Ensambla el operador espacial lineal L para las incógnitas internas

        U = (C_1, ..., C_N, W_1, ..., W_N)^T.

     El Modelo 3 estacionario satisface

        L U_est = b.

        donde 
                   ( L_C  -DaA*I )
          L_st =   (             )    donde L_C es la matriz asociada a la ecuación de coeficientes de C y L_W a la de W.
                   (-DaA*I   L_W )
     Se usa:
        - upwind de primer orden para los términos advectivos;
        - diferencias centradas de segundo orden para la difusión;
        - condición de Neumann homogénea en x = 1 mediante punto fantasma.
    """

    if N < 3:
        raise ValueError("N debe ser al menos 3 para construir una malla útil.")

    DA = 1.0 / caso.PeA
    DC = 1.0 / caso.PeC
    matriz = lil_matrix((2 * N, 2 * N), dtype=float)

    for j in range(N):
        # j = 0,...,N-1 representa el nodo físico i = 1,...,N.
        i = j + 1

        # ----------------------------------------------------
        # Ecuación para C:
        # r_v C_x - DA C_xx + DaA(C - W) = término de contorno.
        # ----------------------------------------------------
        fila_C = j

        if i < N:
            matriz[fila_C, j] += caso.rv / dx + 2.0 * DA / dx**2 + caso.DaA
            matriz[fila_C, j + 1] += -DA / dx**2
            if j - 1 >= 0:
                matriz[fila_C, j - 1] += -caso.rv / dx - DA / dx**2
        else:
            # En x = 1: C_{N+1} = C_{N-1}, por la condición C_x(1)=0.
            matriz[fila_C, j] += caso.rv / dx + 2.0 * DA / dx**2 + caso.DaA
            matriz[fila_C, j - 1] += -caso.rv / dx - 2.0 * DA / dx**2

        # Acoplamiento con W_i.
        matriz[fila_C, N + j] += -caso.DaA

        # ----------------------------------------------------
        # Ecuación para W:
        # W_x - DC W_xx - DaC C + DaC W = término de contorno.
        # ----------------------------------------------------
        fila_W = N + j

        # Acoplamiento con C_i.
        matriz[fila_W, j] += -caso.DaC

        if i < N:
            matriz[fila_W, N + j] += 1.0 / dx + 2.0 * DC / dx**2 + caso.DaC
            matriz[fila_W, N + j + 1] += -DC / dx**2
            if j - 1 >= 0:
                matriz[fila_W, N + j - 1] += -1.0 / dx - DC / dx**2
        else:
            # En x = 1: W_{N+1} = W_{N-1}, por la condición W_x(1)=0.
            matriz[fila_W, N + j] += 1.0 / dx + 2.0 * DC / dx**2 + caso.DaC
            matriz[fila_W, N + j - 1] += -1.0 / dx - 2.0 * DC / dx**2

    return matriz.tocsc()


def vector_contorno(N: int, dx: float, caso: CasoModelo, C_entrada: float, W_entrada: float) -> np.ndarray:
    """
    Construye el vector de contribuciones de contorno para la entrada x=0.

    Las condiciones de Dirichlet son:
        C(0,t) = C_entrada,
        W(0,t) = W_entrada.
    """

    DA = 1.0 / caso.PeA
    DC = 1.0 / caso.PeC
    b = np.zeros(2 * N, dtype=float)

    # Contribuciones del nodo C_0 en la ecuación de C_1.
    b[0] += (caso.rv / dx + DA / dx**2) * C_entrada # representa el nodo físico 1

    # Contribuciones del nodo W_0 en la ecuación de W_1.
    b[N] += (1.0 / dx + DC / dx**2) * W_entrada # representa el nodo físico N+1

    return b


def reconstruir_solucion(U: np.ndarray, N: int, C_entrada: float, W_entrada: float) -> Tuple[np.ndarray, np.ndarray]:
    """Reconstruye los vectores completos C y W, incluyendo el nodo de entrada."""

    C = np.empty(N + 1, dtype=float)
    W = np.empty(N + 1, dtype=float)
    C[0] = C_entrada
    W[0] = W_entrada
    C[1:] = U[:N]
    W[1:] = U[N:]
    return C, W


# ============================================================
#  Modelo 3: problema estacionario lineal acoplado
# ============================================================


def resolver_modelo3_estacionario(N: int, dx: float, caso: CasoModelo) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Resuelve el Modelo 3 estacionario:

        L U_est = b.

    Devuelve C_est, W_est y U_est.
    """

    L = ensamblar_operador_espacial(N, dx, caso)
    b = vector_contorno(N, dx, caso, caso.C_in, caso.W_in)
    U_est = splu(L).solve(b)
    C_est, W_est = reconstruir_solucion(U_est, N, caso.C_in, caso.W_in)
    return C_est, W_est, U_est

# ============================================================
#  GRÁFICA
# ============================================================

N = 100
dx = 1.0 / N


C_est, W_est, _ = resolver_modelo3_estacionario(N, dx, CasoModelo)
x = np.linspace(0, 1, N + 1)

plt.figure(figsize=(8, 5))
plt.plot(x, C_est, 'b-', linewidth=2, label=r'Fase Alveolar ($\hat{C}_A$)')
plt.plot(x, W_est, 'r--', linewidth=2, label=r'Fase Capilar ($\hat{W}$)')

plt.title('Modelo 3: Sistema Estacionario Lineal Acoplado', fontsize=13)
plt.xlabel(r'$\hat{x}$', fontsize=12)
plt.ylabel('Concentración', fontsize=12)
plt.grid(True, alpha=0.3)
plt.legend()
plt.ylim(0, 1.1)

guardar_figura(plt, salida/"01_Solucion" )
