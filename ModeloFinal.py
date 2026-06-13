"""
Modelo Final del biochip alvéolo-capilar

Este script resuelve el sistema:

             C_t + rv C_x = (1/Pe_A)C_xx + Da_A (C - W),
   alfa(W)  (W_t +    W_x) = (1/Pe_C)(D(W)W_x)_x + Da_C (C - W),

con condiciones de Dirichlet en la entrada y condiciones de Neumann homogéneas
en la salida. 

Salidas gráficas principales:
    1) Mapas espacio-temporales C(t,x) en matriz 3x1, con barra de color común.
    2) Mapas espacio-temporales W(t,x) en matriz 3x1, con barra de color común.
    3) Perfiles espaciales C(t_k,x) para tiempos fijos.
    4) Perfiles espaciales W(t_k,x) para tiempos fijos.
    5) Series temporales C(t,x_j) para posiciones fijas.
    6) Series temporales W(t,x_j) para posiciones fijas.
    7) Barrido de la hemoglobina Z_0 en W(4,x_j) en tiempo estacionario.
    8) Validación con modelo 4.

Requisitos:
    numpy, scipy, matplotlib
    

"""
import matplotlib.pyplot as plt
import numpy as np
from dataclasses import dataclass
from scipy.sparse import csc_matrix, eye, lil_matrix, diags
from scipy.sparse.linalg import splu
from typing import List, Tuple, Iterable
import argparse
from pathlib import Path
import re
from dataclasses import replace

plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 14,      # Etiquetas de ejes (ej: "$\hat{x}$")
    'axes.titlesize': 14,      # Títulos
    'legend.fontsize': 13,     # Leyendas
    'xtick.labelsize': 14,     # ← NÚMEROS en eje X
    'ytick.labelsize': 14,     # ← NÚMEROS en eje Y
})

# Constantes
DIRECTORIO_SCRIPT = Path(__file__).resolve().parent
DIRECTORIO_SALIDA_MODELO_FINAL = DIRECTORIO_SCRIPT / "figuras_modelo_final"
TIEMPOS_PERFILES_RELATIVOS = (0.0, 0.10, 0.30, 0.60, 1.0) # Adimensionalizado
POSICIONES_SERIES = (0.25, 0.50, 0.75, 1.00) # Adimensionalizado

# Constantes físicas
DY = 1.4e-11             # m^2/s
DW = 2.4e-9              # m^2/s
DELTA = DY / DW          # Adimensional 
KHILL_MOL = 4.2e-5       # mol/L  
Z0_MOL = 2.33e-3         # mol/L  (hemoglobina total típica)

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
    nHill: float = 2.7           # Adimensional
    kHill: float = 1.0           # Adimensional
    delta: float = DELTA         # Adimensional
    Z0: float = Z0_MOL/KHILL_MOL # Adimensional

@dataclass
class ResultadoModelo:
    caso: CasoModelo
    x: np.ndarray
    t: np.ndarray
    C: np.ndarray
    W: np.ndarray


# ============================================================
# Utilidades generales
# ============================================================

def normalizar_nombre(texto: str) -> str:
    """Convierte un nombre de caso en una cadena segura para nombres de archivo."""

    texto = texto.lower().strip()
    texto = re.sub(r"[^a-z0-9áéíóúñü]+", "_", texto)
    texto = texto.strip("_")
    return texto or "caso"


def guardar_figura(fig: plt.Figure, ruta_base: Path) -> None:
    """Guarda una figura en PNG y PDF con resolución alta."""

    ruta_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(ruta_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(ruta_base.with_suffix(".pdf"), bbox_inches="tight")


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
    '''
       Calcula alpha_i y D_i para cada nodo.
    '''

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
     Ensambla el operador espacial lineal A para las incógnitas 

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
        # j = 0,...,N-1 representa el nodo físico i = 1,...,N.
        i = j + 1  
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


    return matriz.tocsc(), alphas[1:] # Devolvemos alphas de nodos internos

def vector_contorno(N: int, dx: float, caso: CasoModelo, C_entrada: float, W_entrada: float, W_eval)-> np.ndarray:
    """
     Construye el vector de contribuciones de contorno para la entrada x=0.

     Las condiciones de Dirichlet son:
        C(0,t) = C_entrada,
        W(0,t) = W_entrada.
    """

    DA = 1.0 / caso.PeA
    DC = 1.0 / caso.PeC
    b = np.zeros(2 * N)

    # Contribuciones del nodo C_0 en la ecuación de C_1.
    b[0] = (caso.rv / dx + DA / dx**2) * C_entrada

    # Contribuciones del nodo W_0 en la ecuación de W_1 (necesitamos alpha(0) y Ds(0)).
    alphas, Ds = calcular_coeficientes_no_lineales(W_eval, caso)
    b[N] = (alphas[0]/ dx + (DC / dx**2)*((Ds[1] + Ds[0])/2.0)) * W_entrada
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
            b = vector_contorno(N, dx, caso,caso.C_in, caso.W_in, W_eval) # se deja dentro del bucle de Picard por si las condiciones de entrada se cambian a dependientes de t.
            
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
# VALIDACIÓN CON MODELO 4
# ============================================================

# Calculamos primero modelo 4

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

def vector_contorno_lin(N: int, dx: float, caso: CasoModelo, C_entrada: float, W_entrada: float) -> np.ndarray:
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


    C0, W0 = condicion_inicial(x, caso)
    U = np.concatenate([C0[1:], W0[1:]])

    C_hist = np.empty((nt + 1, N + 1), dtype=float)
    W_hist = np.empty_like(C_hist)

    C_hist[0, :] = C0
    W_hist[0, :] = W0

    b = vector_contorno_lin(N, dx, caso, caso.C_in, caso.W_in)

    for n in range(nt):
        rhs = U / dt + b
        U = solver.solve(rhs)

        C, W = reconstruir_solucion(U, N, caso.C_in, caso.W_in)
        C_hist[n + 1, :] = C
        W_hist[n + 1, :] = W


    return ResultadoModelo(
        caso=caso,
        x=x,
        t=t,
        C=C_hist,
        W=W_hist,
    )

# Gráfica de errores

def validacion_limite_Z0(N: int, T: float, dt: float, caso_base: CasoModelo, salida: Path):
    """
    Validación del modelo final en el caso límite Z0=0.

    Cuando Z0=0 la hemoglobina desaparece y el modelo final debe recuperar
    exactamente el Modelo 4. Se superponen ambas soluciones y se calcula
    el error máximo en todo el dominio espacio-temporal.
    """

    # Simular modelo final con Z_0=0. 
    caso_Z0 = replace(caso_base, nombre="Modelo Final Z0=0", Z0=0.0)
    res_final = simular_modelo_final(N, T, dt, caso_Z0)

    # Simular modelo 4.
    res_M4 = simular_modelo4_transitorio(N, T, dt, caso_base)

    # Errores máximos de C y W.
    error_C = float(np.max(np.abs(res_final.C - res_M4.C)))
    error_W = float(np.max(np.abs(res_final.W - res_M4.W)))

    print(f"\n[Validación Z0=0]")
    print(f"  Error máximo en Ĉ_A : {error_C:.2e}")
    print(f"  Error máximo en Ŵ   : {error_W:.2e}")

    # Figura: perfiles estacionarios superpuestos 
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    x = res_final.x

    # Canal alveolar
    axes[0].plot(x, res_final.C[-1, :], color="#2166AC", lw=2,
                 label=r"Modelo final, $\hat{Z}_0=0$")
    axes[0].plot(x, res_M4.C[-1, :], color="#D73027", lw=1.5,
                 linestyle="--", label="Modelo 4")
    axes[0].set_xlabel(r"$\hat{x}$")
    axes[0].set_ylabel(r"$\hat{C}_A$")
    axes[0].set_title(r"Canal alveolar: perfil estacionario $(\hat{t}=4)$")
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    axes[0].text(0.05, 0.10,
                 rf"$\varepsilon_{{max}}={error_C:.1e}$",
                 transform=axes[0].transAxes, fontsize=9,
                 bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    # Canal capilar
    axes[1].plot(x, res_final.W[-1, :], color="#2166AC", lw=2,
                 label=r"Modelo final, $\hat{Z}_0=0$")
    axes[1].plot(x, res_M4.W[-1, :], color="#D73027", lw=1.5,
                 linestyle="--", label="Modelo 4")
    axes[1].set_xlabel(r"$\hat{x}$")
    axes[1].set_ylabel(r"$\hat{W}$")
    axes[1].set_title(r"Canal capilar: perfil estacionario $(\hat{t}=4)$")
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    axes[1].text(0.05, 0.10,
                 rf"$\varepsilon_{{max}}={error_W:.1e}$",
                 transform=axes[1].transAxes, fontsize=9,
                 bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    fig.suptitle(
        r"Validación: recuperación del Modelo 4 en el límite $\hat{Z}_0 \to 0$",
        fontsize=12
    )
    fig.tight_layout()
    guardar_figura(fig, salida / "05_validacion_Z0_cero")


# ============================================================
# Figuras profesionales
# ============================================================

def plot_mapa_3x1(resultados: List[ResultadoModelo], variable: str,  salida: Path):
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
            f"(PeA={res.caso.PeA:g}, PeC={res.caso.PeC:g},"
            f"DaA={res.caso.DaA:g}, DaC={res.caso.DaC:g})", fontsize=14
        )
    
    plt.colorbar(m, ax=ejes, label=etiqueta)
    guardar_figura(fig, salida / f"01_mapa_{variable}_3x1")


def plot_perfiles_espaciales(res: ResultadoModelo, variable: str, salida: Path):
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

    nombre = normalizar_nombre(res.caso.nombre)
    guardar_figura(plt, salida / f"02_perfiles_espaciales_{variable}_{nombre}")

def plot_series_temporales(res: ResultadoModelo, variable: str, salida: Path):
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

    nombre = normalizar_nombre(res.caso.nombre)
    guardar_figura(plt, salida / f"03_series_temporales_{variable}_{nombre}")


def barrido_hemoglobina(N: int, T: float, dt: float, caso_base: CasoModelo, salida: Path):
    """
    Barrido paramétrico en Z0 (concentración de hemoglobina).
    
    Simula desde anemia severa hasta policitemia severa y muestra:
      - Perfil estacionario de W(x) para cada Z0
    
    Los valores de Z0 son adimensionales (Z0_mol / Wref).
      - Z0 = 0     →  sin hemoglobina  (Modelo 4)
      - Z0 ~ 44    →  anemia severa    
      - Z0 ~ 50    →  anemia leve     
      - Z0 ~ 55    →  fisiológico     
      - Z0 ~ 60    →  policitemia leve   
      - Z0 ~ 66    →  policitemia severa    
    """
    from dataclasses import replace

    valores_Z0 = [0, 44, 50, 55, 60, 66]
    etiquetas  = [
        r"$\hat{Z}_0=0$ (sin Hb, Modelo 4)",
        r"$\hat{Z}_0=44$ (anemia severa)",
        r"$\hat{Z}_0=50$ (anemia leve)",
        r"$\hat{Z}_0=55$ (fisiológico)",
        r"$\hat{Z}_0=60$ (policitemia leve)",
        r"$\hat{Z}_0=66$ (policitemia severa)",
    ]
    colores = ["#888780", "#B5D4F4", "#378ADD", "#185FA5", "#0E4E87", "#0B1014"]

    resultados = []
    for z0_val in valores_Z0:
        caso = replace(caso_base, nombre=f"Z0={z0_val}", Z0=float(z0_val))
        res  = simular_modelo_final(N, T, dt, caso)
        resultados.append(res)

    plt.figure()

    # --- Perfil estacionario de W ---
    for res, etiq, col in zip(resultados, etiquetas, colores):
        ls = "--" if res.caso.Z0 == 0 else "-"
        plt.plot(res.x, res.W[-1, :], label=etiq, color=col, linestyle=ls)

    plt.xlabel(r"Posición axial $\hat{x}$")
    plt.ylabel(r"Oxígeno en sangre $\hat{W}$")
    plt.title("Efecto de la concentración de hemoglobina $\hat{Z}_0$ sobre $\hat{W}$")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.legend()

    guardar_figura(plt, salida / f"04_barrido_Hemoglobina")


# ============================================================
# Ejecución de casos 
# ============================================================

def definir_casos() -> List[CasoModelo]:
    """Define tres casos de comparación para las figuras 3x1."""

    return [
        CasoModelo(
            nombre="Caso base",
            rv=1.0,
            PeA=20.0,
            PeC=20.0,
            DaA=3.0,
            DaC=3.0,
            C_in=1.0,
            W_in=0.15,
        ),
        CasoModelo(
            nombre="Acoplamiento alto",
            rv=1.0,
            PeA=20.0,
            PeC=20.0,
            DaA=10.0,
            DaC=10.0,
            C_in=1.0,
            W_in=0.15,
        ),
        CasoModelo(
            nombre="Difusión axial mayor",
            rv=1.0,
            PeA=5.0,
            PeC=5.0,
            DaA=3.0,
            DaC=3.0,
            C_in=1.0,
            W_in=0.15,
        ),
    ]


def ejecutar(N: int, T: float, dt: float):
    """Ejecuta la simulación completa y genera todas las figuras."""

    salida = DIRECTORIO_SALIDA_MODELO_FINAL
    salida.mkdir(parents=True, exist_ok=True)
    casos = definir_casos()

    print("\n=== Modelo final: simulación transitoria no lineal acoplada ===")
    print(f"Malla espacial: N = {N}, dx = {1.0 / N:.5g}")
    print(f"Malla temporal: T = {T:.5g}, dt = {dt:.5g}, pasos = {int(round(T / dt))}")
    print(f"Carpeta de salida: {salida.resolve()}\n")

    resultados: List[ResultadoModelo] = []
    for caso in casos:
        print(f"Resolviendo: {caso.nombre}")
        res = simular_modelo_final(N=N, T=T, dt=dt, caso=caso)
        resultados.append(res)

    print("\nGenerando figuras...")
    
    # Mapas de calor (3 casos juntos)
    plot_mapa_3x1(resultados, "C", salida)
    plot_mapa_3x1(resultados, "W", salida)
    
   
    # Cortes 1D para cada caso.
    for res in resultados:
        plot_perfiles_espaciales(res, "C", salida)
        plot_perfiles_espaciales(res, "W", salida)
        plot_series_temporales(res, "C", salida)
        plot_series_temporales(res, "W", salida)

    # Gráfico del barrido hemoglobina
    barrido_hemoglobina(N, T, dt, casos[0],salida)

    # Gráfico de validación modelo final
    validacion_limite_Z0(N, T, dt, casos[0], salida)
    
    print("\nSimulación terminada. Abriendo ventanas...")
    plt.show()

# ============================================================
# Interfaz de línea de comandos
# ============================================================


def parsear_argumentos() -> argparse.Namespace:
    """Lee argumentos de línea de comandos."""

    parser = argparse.ArgumentParser(
        description="Resuelve el Modelo Finla transitorio no lineal acoplado del biochip por diferencias finitas."
    )
    parser.add_argument("--N", type=int, default=160, help="Número de subintervalos espaciales. Valor por defecto: 160.")
    parser.add_argument("--T", type=float, default=4.0, help="Tiempo final adimensional. Valor por defecto: 4.0.")
    parser.add_argument("--dt", type=float, default=0.0025, help="Paso temporal adimensional. Valor por defecto: 0.0025.")
    #parser.add_argument("--dt", type=float, default=0.1, help="Paso temporal adimensional. Valor por defecto: 0.1.")
    return parser.parse_args()


def main():
    """Punto de entrada principal."""

    args = parsear_argumentos()
    ejecutar(N=args.N, T=args.T, dt=args.dt)


if __name__ == "__main__":
    main()
