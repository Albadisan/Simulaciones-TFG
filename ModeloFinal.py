"""
Modelo Final del biochip alvéolo--capilar

Este script resuelve el sistema:

            C_t + rv C_x = (1/Pe_A)C_xx + Da_A (C - W),
   alfa(W) (W_t +    W_x) = (1/Pe_C)(D(W)W_x)_x + Da_C (C - W),

con condiciones de Dirichlet en la entrada y condiciones de Neumann homogéneas
en la salida. 

Salidas gráficas principales:
    1) Mapas espacio-temporales C(t,x) en matriz 3x1, con barra de color común.
    2) Mapas espacio-temporales W(t,x) en matriz 3x1, con barra de color común.
    3) Perfiles espaciales C(t_k,x) para tiempos fijos.
    4) Perfiles espaciales W(t_k,x) para tiempos fijos.
    5) Series temporales C(t,x_j) para posiciones fijas.
    6) Series temporales W(t,x_j) para posiciones fijas.
    

"""
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass
from scipy.sparse import csc_matrix, eye, lil_matrix, diags
from scipy.sparse.linalg import splu
from typing import List, Tuple, Iterable
import argparse

# Constantes
TIEMPOS_PERFILES_RELATIVOS = (0.0, 0.10, 0.30, 0.60, 1.0)
POSICIONES_SERIES = (0.25, 0.50, 0.75, 1.00)

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
    nHill: float = 2.7
    kHill: float = 4.2e-5
    DY: float = 1.4e-11
    DW: float = 2.4e-9
    delta: float = DY/DW
    Z0: float = 2.0

@dataclass
class ResultadoModelo:
    caso: CasoModelo
    x: np.ndarray
    t: np.ndarray
    C: np.ndarray
    W: np.ndarray

# --- FUNCIONES AUXILIARES ---
def indices_mas_cercanos(malla: np.ndarray, valores: Iterable[float]) -> List[int]:
    return [int(np.argmin(np.abs(malla - valor))) for valor in valores]

def datos_variable(res: ResultadoModelo, variable: str) -> Tuple[np.ndarray, np.ndarray]:
    if variable == "C":
        return res.C
    return res.W

###############################
# --- Funciones No Lineales ---
###############################

def derivada_funcion_Hill(W: np.ndarray, caso: CasoModelo) -> np.ndarray:
    '''
       Derivada de la curva de Hill: 
        
         d/dW [ W^n / (W^n + k^n) ] 
    
    '''
    numerador = caso.nHill * (caso.kHill**caso.nHill) * (W**(caso.nHill - 1))
    denominador = (W**caso.nHill + caso.kHill**caso.nHill)**2
    derivada= numerador / denominador
    return derivada

def calcular_coeficientes_no_lineales(W_vector, caso):
    """Calcula alpha_i y D_i para cada nodo."""

    # alpha_i = 1 + 4 * Z0 * f'(W_i)
    alphas = 1.0 + 4.0 * caso.Z0 * derivada_funcion_Hill(W_vector,caso)
    
    # D_i = 1 + 4 * Z0 * delta * f'(W_i)
    Ds = 1.0 + 4.0 * caso.Z0 * caso.delta * derivada_funcion_Hill(W_vector,caso)
    return alphas, Ds

# ============================================================
# Ensamblaje del operador espacial 
# ============================================================


def ensamblar_sistema_picard(N, dx, caso, W_eval):
    
    """
     Ensambla el operador espacial lineal L para las incógnitas 

        U = (C_0,C_1, ..., C_N, W_1, ..., W_N, W_N+1)^T.

    La ecuación se escribe como

         alfa_i dU/dt + L(D) U = b(t).

    mientras que el Modelo 4 estacionario satisface (alfa_i=1, D=1)

        dU/dt + L U = b(t).

    Se usa:
        - upwind de primer orden para los términos advectivos;
        - diferencias centradas de segundo orden para la difusión;
        - condición de Neumann homogénea en x = 1 mediante punto fantasma.
    """
    
    if N < 3:
        raise ValueError("N debe ser al menos 3 para construir una malla útil.")


    DA = 1.0 / caso.PeA
    DC = 1.0 / caso.PeC
    alphas, Ds = calcular_coeficientes_no_lineales(W_eval, caso)

    matriz = lil_matrix((2 * N, 2 * N), dtype=float)


    for j in range(N):
        i = j + 1  # Índice físico (1 a N)
        fila_W = N + j
        
        # ----------------------------------------------------
        # Ecuación para C:
        # r_v C_x - DA C_xx + DaA(C - W) = término de contorno.
        # ----------------------------------------------------   
        fila_C = j
 
        # Generamos la matriz 
      
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

        # -------------------------------------------------------------------------------
        # Ecuación para W:
        # alfa_i W_x - DC( D(W)W_x)_x - DaC C + DaC W = término de contorno.
        # -------------------------------------------------------------------------------

        fila_W = N + j

        # Acoplamiento con C_i.
        matriz[fila_W, j] += -caso.DaC

        alpha_i = alphas[i]
        D_mas = (Ds[i] + Ds[i+1])/2.0 if i < N else Ds[i]
        D_menos = (Ds[i-1] + Ds[i])/2.0

        if i < N:
            matriz[fila_W, N + j] += alpha_i/dx + (DC/dx**2)*(D_mas + D_menos) + caso.DaC
            matriz[fila_W, N + j + 1] += -(DC / dx**2)*D_mas
            if j - 1 >= 0:
                matriz[fila_W, N + j - 1] += -alpha_i/dx - (DC/dx**2)*D_menos
        else:
            # En x = 1: W_{N+1} = W_{N-1}, por la condición W_x(1)=0.
            matriz[fila_W, N + j] += alpha_i/dx + (DC/dx**2)*(D_mas + D_menos) + caso.DaC
            matriz[fila_W, N + j - 1] += -alpha_i/ dx - (DC/dx**2)*(D_mas + D_menos)


    return matriz.tocsc(), alphas[1:] # Retornamos alphas de nodos internos

def vector_contorno(N, dx, caso, C_in, W_in, W_eval):
    b = np.zeros(2 * N)

    alphas, Ds = calcular_coeficientes_no_lineales(W_eval, caso)

    b[0] = (caso.rv / dx + (1.0/caso.PeA) / dx**2) * C_in
    b[N] = (alphas[0]/ dx + ((1.0/caso.PeC) / dx**2)*((Ds[1] + Ds[0])/2.0)) * W_in
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
# Modelo Final
# ============================================================


def condicion_inicial(x: np.ndarray, caso: CasoModelo) -> Tuple[np.ndarray, np.ndarray]:
    """Define una condición inicial suave y no estacionaria."""
    C0 = caso.C_in * (0.20 + 0.80 * np.exp(-4.0 * x))
    W0 = caso.W_in * np.ones_like(x)
    C0[0] = caso.C_in
    W0[0] = caso.W_in
    return C0, W0


def simular_modelo_final(N: int, T: float, dt: float, caso: CasoModelo,tol_picard=1e-6) -> ResultadoModelo:

    """
    Resuelve el Modelo final

    En cada paso temporal se resuelve
                                                           ( I      0  )
        (A/dt + L(D)) U^{n+1} = A/dt*U^n + b.    Donde A = (           )
                                                           ( 0   alfa*I)
    """

    if T <= 0:
        raise ValueError("El tiempo final T debe ser positivo.")
    if dt <= 0:
        raise ValueError("El paso temporal dt debe ser positivo.")

    x = np.linspace(0.0, 1.0, N + 1)
    dx = 1.0 / N
    nt = int(np.round(T / dt))
    t = np.linspace(0.0, nt * dt, nt + 1)

    C0, W0 = condicion_inicial(x, caso)
    U = np.concatenate([C0[1:], W0[1:]])

    C_hist = np.empty((nt + 1, N + 1), dtype=float)
    W_hist = np.empty_like(C_hist)
    

    C_hist[0, :] = C0
    W_hist[0, :] = W0

    
    for n in range(nt):
        U_old = U.copy()
        U_iter = U.copy()
        
        # Picard para acoplamiento no lineal
        for k in range(15):
            _, W_eval = reconstruir_solucion(U_iter, N, caso.C_in, caso.W_in)
            L, alphas_nodos = ensamblar_sistema_picard(N, dx, caso, W_eval)
            b = vector_contorno(N, dx, caso,caso.C_in, caso.W_in, W_eval)
            
            M = diags(np.concatenate([np.ones(N), alphas_nodos]), format="csc")
            A_mat = M / dt + L
            rhs = (M / dt) @ U_old + b
            
            U_new = splu(A_mat).solve(rhs)
            
            if np.linalg.norm(U_new - U_iter, ord=np.inf) < tol_picard:
                U_iter = U_new
                break
            U_iter = U_new
        
        U = U_iter
        C_hist[n+1], W_hist[n+1] = reconstruir_solucion(U, N, caso.C_in, caso.W_in)


    return ResultadoModelo(
        caso=caso,
        x=x,
        t=t,
        C=C_hist,
        W=W_hist,
    )

# ============================================================
# Figuras profesionales
# ============================================================

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
    vmin = min(float(np.min(z)) for z in datos)
    vmax = max(float(np.max(z)) for z in datos)

    etiqueta = r"$\hat C_A(\hat t,\hat x)$" if variable == "C" else r"$\hat W(\hat t,\hat x)$"
    titulo_variable = r"$\hat C_A$" if variable == "C" else r"$\hat W$"

    fig, ejes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

    for ax, res, z in zip(ejes, resultados, datos):
        m = ax.pcolormesh(res.t, res.x, z.T, shading="auto", vmin=vmin, vmax=vmax)
        ax.set_ylabel(r"$\hat x$")
        ax.set_title(
            f"{titulo_variable}: {res.caso.nombre} "
            f"(PeA={res.caso.PeA:g}, PeC={res.caso.PeC:g}, "
            f"DaA={res.caso.DaA:g}, DaC={res.caso.DaC:g})"
        )
    
    plt.colorbar(m, ax=ejes, label=etiqueta)

def plot_perfiles_espaciales(res: ResultadoModelo, variable: str):
    """Dibuja perfiles espaciales para varios tiempos fijos."""

    if variable not in {"C", "W"}:
        raise ValueError("La variable debe ser 'C' o 'W'.")
    
    datos = datos_variable(res, variable)
    etiqueta = r"$\hat C_A(\hat t_k,\hat x)$" if variable == "C" else r"$\hat W(\hat t_k,\hat x)$"

    tiempos_objetivo = [factor * res.t[-1] for factor in TIEMPOS_PERFILES_RELATIVOS]
    indices_t = indices_mas_cercanos(res.t, tiempos_objetivo)

    plt.figure()

    for idx in indices_t:
        plt.plot(res.x, datos[idx, :], label=f"$\hat t$={res.t[idx]:.2f}")
    plt.xlabel(r"$\hat x$")
    plt.ylabel(etiqueta)
    plt.title(f"Perfiles espaciales de {variable}: {res.caso.nombre}")
    plt.legend()

def plot_series_temporales(res: ResultadoModelo, variable: str):
    """Dibuja perfiles espaciales para varios tiempos fijos."""

    if variable not in {"C", "W"}:
        raise ValueError("La variable debe ser 'C' o 'W'.")
    
    datos = datos_variable(res, variable)
    etiqueta = r"$\hat C_A(\hat t,\hat x_j)$" if variable == "C" else r"$\hat W(\hat t,\hat x_j)$"

    plt.figure()

    indices = indices_mas_cercanos(res.x, POSICIONES_SERIES)

    for idx in indices:
        plt.plot(res.t, datos[:, idx], label=f"$\hat x$ ={res.x[idx]:.2f}")

    plt.xlabel(r"$\hat t$")
    plt.ylabel(etiqueta)
    plt.title(f"Evolución temporal de {variable} - Caso: {res.caso.nombre}")
    plt.legend()


# ============================================================
# Ejecución de casos y guardado de datos
# ============================================================

def definir_casos():
    return [
        CasoModelo("Base", PeA=20, PeC=20, DaA=3, DaC=3),
        CasoModelo("Acoplamiento Alto", PeA=20, PeC=20, DaA=10, DaC=10),
        CasoModelo("Difusión Mayor", PeA=5, PeC=5, DaA=3, DaC=3)
    ]

def ejecutar(N: int, T: float, dt: float):
    casos = definir_casos()
    resultados = [simular_modelo_final(N, T, dt, c) for c in casos]
    
    # Mapas de calor (3 casos juntos)
    plot_mapa_3x1(resultados, "C")
    plot_mapa_3x1(resultados, "W")
    
   
    # Gráficos específicos del primer caso (Base)
    res0 = resultados[0]
    plot_perfiles_espaciales(res0, "C")
    plot_perfiles_espaciales(res0, "W")
    plot_series_temporales(res0, "C")
    plot_series_temporales(res0, "W")
   
    
    print("\nSimulación terminada. Abriendo ventanas...")
    plt.show()

def main():
    args = parsear_argumentos()
    ejecutar(N=args.N, T=args.T, dt=args.dt)

def parsear_argumentos():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=100)
    parser.add_argument("--T", type=float, default=4.0)
    parser.add_argument("--dt", type=float, default=0.01)
    return parser.parse_args()

if __name__ == "__main__":
    main()
