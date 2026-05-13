"""
Modelo 4 del biochip alvéolo--capilar: sistema transitorio lineal acoplado

Este script resuelve el sistema:

   C_t + rv C_x = (1/Pe_A)C_xx + Da_A (C - W),
   W_t +    W_x = (1/Pe_C)W_xx + Da_C (C - W),

con condiciones de Dirichlet en la entrada y condiciones de Neumann homogéneas
en la salida. Además, calcula el Modelo 3 estacionario correspondiente y verifica
la convergencia del Modelo 4 hacia dicho estado estacionario cuando t es grande.

Salidas gráficas principales:
    1) Mapas espacio-temporales C(t,x) en matriz 3x1, con barra de color común.
    2) Mapas espacio-temporales W(t,x) en matriz 3x1, con barra de color común.
    3) Perfiles espaciales C(t_k,x) para tiempos fijos.
    4) Perfiles espaciales W(t_k,x) para tiempos fijos.
    5) Series temporales C(t,x_j) para posiciones fijas.
    6) Series temporales W(t,x_j) para posiciones fijas.
    7) Error relativo de convergencia Modelo 4 -> Modelo 3.
    8) Error espacial final respecto al estado estacionario.

"""

import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass
from scipy.sparse import csc_matrix, eye, lil_matrix
from scipy.sparse.linalg import splu
from typing import List, Tuple, Iterable
import argparse

# Constantes
TIEMPOS_PERFILES_RELATIVOS = (0.0, 0.10, 0.30, 0.60, 1.0)
POSICIONES_SERIES = (0.25, 0.50, 0.75, 1.00)

# ============================================================
# Definición de parámetros y estructuras de datos
# ============================================================

@dataclass(frozen=True)
class CasoModelo:
    nombre: str
    rv: float = 1.0
    PeA: float = 20.0
    PeC: float = 20.0
    DaA: float = 3.0
    DaC: float = 3.0
    C_in: float = 1.0
    W_in: float = 0.15

@dataclass
class ResultadoModelo:
    caso: CasoModelo
    x: np.ndarray
    t: np.ndarray
    C: np.ndarray
    W: np.ndarray
    C_est: np.ndarray
    W_est: np.ndarray
    error_relativo: np.ndarray

# ============================================================
# Utilidades generales
# ============================================================

def indices_mas_cercanos(malla: np.ndarray, valores: Iterable[float]) -> List[int]:
    return [int(np.argmin(np.abs(malla - valor))) for valor in valores]

def datos_variable(res: ResultadoModelo, variable: str) -> Tuple[np.ndarray, np.ndarray]:
    if variable == "C":
        return res.C, res.C_est
    return res.W, res.W_est

# ============================================================
# Ensamblaje del operador espacial común a los Modelos 3 y 4
# ============================================================

def ensamblar_operador_espacial(N: int, dx: float, caso: CasoModelo) -> csc_matrix:

    """
    Ensambla el operador espacial lineal L para las incógnitas internas

        U = (C_1, ..., C_N, W_1, ..., W_N)^T.

    La ecuación transitoria se escribe como

        dU/dt + L U = b(t),

    mientras que el Modelo 3 estacionario satisface

        L U_est = b.

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
    b[0] += (caso.rv / dx + DA / dx**2) * C_entrada

    # Contribuciones del nodo W_0 en la ecuación de W_1.
    b[N] += (1.0 / dx + DC / dx**2) * W_entrada

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
# Modelo 3: problema estacionario lineal acoplado
# ============================================================


def resolver_modelo3_estacionario(N, dx, caso):
    L = ensamblar_operador_espacial(N, dx, caso)
    b = vector_contorno(N, dx, caso, caso.C_in, caso.W_in)
    U_est = splu(L).solve(b)
    C_est, W_est = reconstruir_solucion(U_est, N, caso.C_in, caso.W_in)
    return C_est, W_est, U_est

# ============================================================
# Modelo 4: problema transitorio lineal acoplado
# ============================================================

def condicion_inicial(x: np.ndarray, caso: CasoModelo) -> Tuple[np.ndarray, np.ndarray]:
    """
    Define una condición inicial suave y no estacionaria.

    La elección no pretende representar un caso clínico concreto; su objetivo es
    generar una relajación transitoria clara hacia el Modelo 3 estacionario.
    """

    C0 = caso.C_in * (0.20 + 0.80 * np.exp(-4.0 * x))
    W0 = caso.W_in * np.ones_like(x)

    # Se impone compatibilidad exacta con las condiciones de entrada.
    C0[0] = caso.C_in
    W0[0] = caso.W_in
    return C0, W0


def simular_modelo4_transitorio(N: int, T: float, dt: float, caso: CasoModelo) -> ResultadoModelo:
    """
    Resuelve el Modelo 4 mediante Euler implícito en tiempo.

    En cada paso temporal se resuelve

        (I/dt + L) U^{n+1} = U^n/dt + b.

    Como los coeficientes son constantes, la matriz se factoriza una sola vez.
    """

    if T <= 0:
        raise ValueError("El tiempo final T debe ser positivo.")
    if dt <= 0:
        raise ValueError("El paso temporal dt debe ser positivo.")
    
    x = np.linspace(0.0, 1.0, N + 1)
    dx = 1.0 / N
    nt = int(np.round(T / dt))
    t = np.linspace(0.0, nt * dt, nt + 1)

    L = ensamblar_operador_espacial(N, dx, caso)
    A = eye(2 * N, format="csc") / dt + L
    solver = splu(A)

    C_est, W_est, U_est = resolver_modelo3_estacionario(N, dx, caso)

    C0, W0 = condicion_inicial(x, caso)
    U = np.concatenate([C0[1:], W0[1:]])

    C_hist = np.empty((nt + 1, N + 1), dtype=float)
    W_hist = np.empty_like(C_hist)
    error_relativo = np.empty(nt + 1, dtype=float)

    C_hist[0, :] = C0
    W_hist[0, :] = W0
    denominador = max(1.0, np.linalg.norm(U_est, ord=np.inf))
    error_relativo[0] = np.linalg.norm(U - U_est, ord=np.inf) / denominador

    b = vector_contorno(N, dx, caso, caso.C_in, caso.W_in)

    for n in range(nt):
        rhs = U / dt + b
        U = solver.solve(rhs)

        C_hist[n+1], W_hist[n+1] = reconstruir_solucion(U, N, caso.C_in, caso.W_in)

        error_relativo[n + 1] = np.linalg.norm(U - U_est, ord=np.inf) / denominador

    return ResultadoModelo(
        caso=caso,
        x=x,
        t=t,
        C=C_hist,
        W=W_hist,
        C_est=C_est,
        W_est=W_est,
        error_relativo=error_relativo,
    )

###############################
# --- FUNCIONES DE GRÁFICOS ---
###############################

def plot_mapa_3x1(resultados: List[ResultadoModelo], variable: str):
    """
    Genera una figura 3x1 con mapas espacio-temporales y una barra de color común.

    variable debe ser "C" o "W".
    """

    if variable not in {"C", "W"}:
        raise ValueError("La variable debe ser 'C' o 'W'.")
    if len(resultados) != 3:
        raise ValueError("Esta figura está diseñada para exactamente tres casos.")
    
    datos = [getattr(res, variable) for res in resultados]
    vmin = min(np.min(z) for z in datos)
    vmax = max(np.max(z) for z in datos)

    titulo_variable = r"$\hat C_A$" if variable == "C" else r"$\hat W$"
    etiqueta_barra = r"$\hat C_A(\hat t, \hat x)$" if variable == "C" else r"$\hat W (\hat t, \hat x)$"

    fig, ejes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    for ax, res, z in zip(ejes, resultados, datos):
        m = ax.pcolormesh(res.t, res.x, z.T, shading="auto", vmin=vmin, vmax=vmax)
        ax.set_ylabel(r"$\hat x$")
        ax.set_title( f"{titulo_variable}: {res.caso.nombre} "
                      f"(PeA={res.caso.PeA:g}, PeC={res.caso.PeC:g}, "
                      f"DaA={res.caso.DaA:g}, DaC={res.caso.DaC:g})"
                     )
    plt.colorbar(m, ax=ejes, label=etiqueta_barra)

def plot_perfiles_espaciales(res: ResultadoModelo, variable: str):
    """Dibuja perfiles espaciales para varios tiempos fijos."""

    if variable not in {"C", "W"}:
        raise ValueError("La variable debe ser 'C' o 'W'.")
    
    datos, est = datos_variable(res, variable)

    etiqueta = r"$\hat C_A(\hat t_k,\hat x)$" if variable == "C" else r"$\hat W(\hat t_k,\hat x)$"
    etiqueta_est = r"Modelo 3: $\hat C_{A,\mathrm{st}}$" if variable == "C" else r"Modelo 3: $\hat W_{\mathrm{st}}$"

    plt.figure()

    etiqueta_est = r"Modelo 3: $\hat C_{A,\mathrm{st}}$" if variable == "C" else r"Modelo 3: $\hat W_{\mathrm{st}}$"

    tiempos_objetivo = [factor * res.t[-1] for factor in TIEMPOS_PERFILES_RELATIVOS]
    indices_t = indices_mas_cercanos(res.t, tiempos_objetivo)

    for idx in indices_t:
        plt.plot(res.x, datos[idx, :], label=f"t={res.t[idx]:.2f}")

    plt.plot(res.x, est , linestyle="--", linewidth=2.0, label=etiqueta_est)
    plt.xlabel(r"$\hat x$")
    plt.ylabel(etiqueta)
    plt.title(f"Perfiles espaciales de {variable} - Caso: {res.caso.nombre}")
    plt.legend()


def plot_series_temporales(res: ResultadoModelo, variable: str):
    """Dibuja la evolución temporal para varias posiciones fijas."""

    if variable not in {"C", "W"}:
        raise ValueError("La variable debe ser 'C' o 'W'.")
    
    datos, est = datos_variable(res, variable)
    etiqueta = r"$\hat C_A(\hat t,\hat x_j)$" if variable == "C" else r"$\hat W(\hat t,\hat x_j)$"

    plt.figure()

    indices = indices_mas_cercanos(res.x, POSICIONES_SERIES)

    for idx in indices:
        plt.plot(res.t, datos[:, idx], label=f"x={res.x[idx]:.2f}")

    plt.xlabel(r"$\hat t$")
    plt.ylabel(etiqueta)
    plt.title(f"Evolución temporal de {variable} - Caso: {res.caso.nombre}")
    plt.legend()


def plot_convergencia(resultados: List[ResultadoModelo]):
    """Dibuja el error relativo de convergencia."""
    plt.figure(figsize=(7.2, 4.8))

    for res in resultados:
        plt.semilogy(res.t, res.error_relativo, label=res.caso.nombre)

    plt.xlabel(r"$\hat t$")
    plt.ylabel(r"$E_{\mathrm{rel}}(\hat t)$")
    plt.title("Convergencia del Modelo 4 hacia el Modelo 3")
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()


def plot_error_final(res: ResultadoModelo):
    """Dibuja el error espacial final respecto al estado estacionario."""

    error_C = np.abs(res.C[-1, :] - res.C_est)
    error_W = np.abs(res.W[-1, :] - res.W_est)

    fig, ejes = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)

    ejes[0].semilogy(res.x, error_C, color='blue')
    ejes[0].semilogy(res.x, error_C)
    ejes[0].set_ylabel(r"$|\hat C_A(T,\hat x)-\hat C_{A,\mathrm{st}}(\hat x)|$")
    ejes[0].set_title(f"Error espacial final: {res.caso.nombre}")
    ejes[0].grid(True, which="both", alpha=0.3)

    ejes[1].semilogy(res.x, error_W)
    ejes[1].set_xlabel(r"$\hat x$")
    ejes[1].set_ylabel(r"$|\hat W(T,\hat x)-\hat W_{\mathrm{st}}(\hat x)|$")
    ejes[1].grid(True, which="both", alpha=0.3)

    for ax in ejes: ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()

# ============================================================
# Ejecución de casos y guardado de datos
# ============================================================


def definir_casos():
    """Define tres casos de comparación para las figuras 3x1."""
    
    return [
        CasoModelo("Base", PeA=20, PeC=20, DaA=3, DaC=3),
        CasoModelo("Acoplamiento Alto", PeA=20, PeC=20, DaA=10, DaC=10),
        CasoModelo("Difusión Mayor", PeA=5, PeC=5, DaA=3, DaC=3)
    ]

def ejecutar(N: int, T: float, dt: float) -> None:
    casos = definir_casos()
    resultados = [simular_modelo4_transitorio(N, T, dt, c) for c in casos]
    
    # Mapas de calor (3 casos juntos)
    plot_mapa_3x1(resultados, "C")
    plot_mapa_3x1(resultados, "W")
    
    # Convergencia (Todos los casos en una gráfica)
    plot_convergencia(resultados)
    
    # Gráficos específicos del primer caso (Base)
    res0 = resultados[0]
    plot_perfiles_espaciales(res0, "C")
    plot_perfiles_espaciales(res0, "W")
    plot_series_temporales(res0, "C")
    plot_series_temporales(res0, "W")
    plot_error_final(res0)
    
    print("\nSimulación terminada. Abriendo ventanas...")
    plt.show() 

def main():
    args = parsear_argumentos()
    ejecutar(N=args.N, T=args.T, dt=args.dt)

def parsear_argumentos():
    parser = argparse.ArgumentParser(
        description="Resuelve el Modelo 4 transitorio lineal acoplado del biochip por diferencias finitas."
    )
    parser.add_argument("--N", type=int, default=160, help="Número de subintervalos espaciales. Valor por defecto: 160.")
    parser.add_argument("--T", type=float, default=4.0, help="Tiempo final adimensional. Valor por defecto: 4.0.")
    parser.add_argument("--dt", type=float, default=0.0025, help="Paso temporal adimensional. Valor por defecto: 0.0025.")
    return parser.parse_args()

if __name__ == "__main__":
    main()
