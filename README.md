
   ### Modelo 1D de Reacción-Advección-Difusión para Intercambio Gaseoso en Biochips Pulmonares: Análisis y Extensiones Fisiológicas

Repositorio de código científico que acompaña al Trabajo de Fin de Grado (TFG). En este proyecto se simula el sistema acoplado transitorio y no lineal para el transporte de gases en un biochip alvéolo-capilar, relajando hipótesis de forma progresiva en cuatro etapas para garantizar la validación del esquema numérico.

En el desarrollo del TFG se pretende simular en Python el sistema acoplado:

$$
\frac{\partial \hat{C}_{A}}{\partial \hat{t}} + r_{v}\ \frac{ \partial \hat{C}_A}{\partial \hat{x}}
= \frac{1}{\mathrm{Pe}_{A}}\ \frac{\partial^2 \hat{C}_{A}}{\partial \hat{x}^2} - \mathrm{Da}_A\ (\hat{C}_A - \hat{W})
$$

$$ 
 \hat{\alpha}(\hat{W})\left[\frac{\partial\hat{W}}{\partial {\hat{t}}} + \frac{\partial\hat{W}}{\partial \hat{x}} \right] = \frac{1}{\mathrm{Pe}_C}\frac{\partial}{\partial \hat{x}} \left(\hat{D}(\hat{W})\frac{\partial \hat{W}}{\partial \hat{x}} \right) + Da_C(\hat{C}_A -\hat{W})
$$

    Nota: Todas las variables son adimensionales; se omite el sombrero en lo sucesivo.

# Instrucciones y comandos de ejecución:
    Para preparar el entorno de Python e instalar todas las librerías científicas necesarias, ejecuta en tu terminal:
    ```bash
       pip install -r requirements.txt

    Para reproducir las simulaciones y regenerar todas las figuras que aparecen en la memoria del TFG, ejecuta los scripts en tu consola en     el siguiente orden secuencial:
    ```bash
      python MODELO1.py
      python MODELO2.py
      python MODELO3.py
      python MODELO4.py
      python ModeloFinal.py

# Estructura del repositorio:
     ├── figuras_modelo1              # Figuras modelo 1: advección–reacción puro.
     ├── figuras_modelo2              # Figuras modelo 2: advección–difusión–reacción.
     ├── figuras_modelo3              # Figuras modelo 3: sistema estacionario lineal acoplado.
     ├── figuras_modelo4              # Figuras modelo 4: sistema transitorio lineal acoplado.
     ├── figuras_modelo_final         # Figuras modelo final: sistema no lineal con hemoglobina.
     ├── MODELO1.py                   # Modelo 1: advección–reacción puro.
     ├── MODELO2.py                   # Modelo 2: advección–difusión–reacción.
     ├── MODELO3.py                   # Modelo 3: sistema estacionario lineal acoplado.
     ├── MODELO4.py                   # Modelo 4: sistema transitorio lineal acoplado.
     ├── ModeloFinal.py               # Modelo final: sistema no lineal con hemoglobina.
     ├── requirements.txt             # Dependencias Python
     └── README.md

# MODELOS

- MODELO 1 : Primera etapa. 
  
    Se imponen 4 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $\partial / \partial t = 0$
     3) Concentración alveolar fija: $C_A(x) = C_A^0$
     4) No hay difusión axial: $Pe_C \to \inf $

    Se resuelve por tanto la ecuación:

     $$\frac{d W}{dx} = \mathrm{Da}_C (C_A^0 − W)$$
  
     $$W(0) = W_{in}$$

- MODELO 2 : Segunda etapa. 
  
    Se imponen 4 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $\partial / \partial t = 0$
     3) Concentración alveolar fija: $C_A(x) = C_A^0$
     4) Hay difusión axial: $Pe_C$ finito

    Se resuelve por tanto la ecuación:

     $$\frac{d W}{dx} = \frac{1}{\mathrm{Pe}_{C}} \frac{d^2 W}{dx^2} + \mathrm{Da}_C (C_A^0 − W)$$
  
     En $x=0$: $\hspace{1cm} W(0) = W_{in} \hspace{1cm}$  (Dirichlet en la entrada)
  
     En $x=1$: $\hspace{1cm} \frac{d W}{dx}(1) = 0 \hspace{1cm}$ (Neumann en la salida)
 
- MODELO 3 : Tercera etapa. 
  
    Se imponen 2 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $\partial / \partial t = 0$

    Se resuelve por tanto el sistema:

     $$r_v \frac{d C}{d x}= \frac{1}{\mathrm{Pe}_{A}}\frac{d^2 C}{d x^2} - \mathrm{Da}_A(C_A − W)$$
     
     $$\frac{d W}{dx}= \frac{1}{\mathrm{Pe}_{C}} \frac{d^2 W}{dx^2} + \mathrm{Da}_C (C_A^0 − W)$$

     En $x=0$: $\hspace{1cm} C_A (0) = C_{A,in} \hspace{1cm} W(0) = W_{in} \hspace{1cm}$  (Dirichlet en la entrada)
  
     En $x=1$: $\hspace{1cm} \frac{d C_A}{dx}(1) = 0 \hspace{1cm} \frac{d W}{dx}(1) = 0 \hspace{1cm}$ (Neumann en la salida)

- MODELO 4 : Cuarta etapa. 
  
    Se impone 1 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$

    Se resuelve por tanto el sistema:

     $$\frac{\partial C}{\partial t} + r_v\frac{\partial C}{\partial x} = \frac{1}{\mathrm{Pe}_{A}}\frac{\partial^2 C}{\partial x^2} - \mathrm{Da}_A(C_A − W)$$

     $$\frac{\partial W}{\partial t} + \frac{\partial W}{\partial x} = \frac{1}{\mathrm{Pe}_{C}} \frac{\partial^2 W}{\partial x^2} + \mathrm{Da}_C (C_A^0 − W)$$

     Condición inicial : $C_A (x,0) = C_{A,0}(x)$
  
     En $x=0$: $\hspace{1cm} C_A (0,t) = C_{A,in} \hspace{1cm} W(0,t) = W_{in} \hspace{1cm}$ (Dirichlet en la entrada)
  
     En $x=1$: $\hspace{1cm} \frac{\partial C_A}{\partial x}(1,t) = 0 \hspace{1cm} \frac{\partial W}{\partial x}(1,t) = 0 \hspace{1cm}$ (Neumann en la salida)
  
 - MODELO FINAL : simulación final.

    Se añaden las funciones no lineales
   
    $\alpha(W) = 1 + 4 Z_0 f'(W)$
  
    $D(W) = 1 + 4 \delta Z_0 f'(W)$

    Donde $f'(W)$ representa la derivada de la función de Hill.
   
    Y se tiene el sistema original que pretendíamos simular.

# REQUERIMIENTOS

Python 3.14 o superior. Instalar dependencias con:

    pip install -r requirements.txt

# RESULTADOS GRÁFICOS

## MODELO 1

Genera en /figuras_modelo1:

| Archivo / Visualización | Contenido |
|---|---|
|![Solución modelo 1](figuras_modelo1/01_Solucion.png)|Comparación solución analítica vs numérica (N=100 y N=1000).|
|![Convergencia modelo 1](figuras_modelo1/02_Convergencia.png)|Estudio de convergencia log-log (pendiente obtenida ≈ −1).|


## MODELO 2

Genera en /figuras_modelo2:

| Archivo / Visualización | Contenido |
|---|---|
|![Solución modelo 2](figuras_modelo2/01_Solucion.png)|Comparación solución analítica vs numérica.|
|![Convergencia modelo 2](figuras_modelo2/02_Convergencia.png)|Estudio de convergencia log-log (pendiente obtenida ≈ −2).|


## MODELO 3

Genera en /figuras_modelo3:

| Archivo / Visualización | Contenido |
|---|---|
|![Solución modelo 3](figuras_modelo3/01_Solucion.png)|Perfiles estacionarios de \hat{C}_A y \hat{W}.|

## MODELO 4

Genera en figuras_modelo4/:

| Archivo / Visualización | Contenido |
|---|---|
| ![Mapa C](figuras_modelo4/01_mapa_C_3x1.png) | Mapas espacio-temporales de $C_A$ (3 casos) |
| ![Mapa W](figuras_modelo4/01_mapa_W_3x1.png) | Mapas espacio-temporales de W (3 casos) |
| ![Perfil espacial C caso base](figuras_modelo4/02_perfiles_espaciales_C_base.png) | Perfiles espaciales de $C_A$ (caso base) |
| ![Perfil espacial W caso base](figuras_modelo4/02_perfiles_espaciales_W_base.png) | Perfiles espaciales de W (caso base) |
| ![Serie temporal C caso base](figuras_modelo4/03_series_temporales_C_base.png) | Series temporales de $C_A$ (caso base) |
| ![Serie temporal W caso base](figuras_modelo4/03_series_temporales_W_base.png) | Series temporales de W (caso base) |
| ![Convergencia al modelo 3](figuras_modelo4/04_convergencia_modelo4_hacia_modelo3.png) | Error relativo de convergencia al estacionario |
| ![Error espacial final respecto modelo 3](figuras_modelo4/05_error_espacial_final_base.png) | Error espacial final respecto al Modelo 3 |


## MODELO FINAL

Genera en figuras_modelo_final/:

| Archivo | Contenido |
|---|---|
| ![Mapa C](figuras_modelo_final/01_mapa_C_3x1.png) | Mapas espacio-temporales de $C_A$ (3 casos) |
| ![Mapa W](figuras_modelo_final/01_mapa_W_3x1.png) | Mapas espacio-temporales de W (3 casos) |
| ![Perfil espacial C caso base](figuras_modelo_final/02_perfiles_espaciales_C_caso_base.png) | Perfiles espaciales de $C_A$ (caso base) |
| ![Perfil espacial W caso base](figuras_modelo_final/02_perfiles_espaciales_W_caso_base.png) | Perfiles espaciales de W (caso base) |
| ![Serie temporal C caso base](figuras_modelo_final/03_series_temporales_C_caso_base.png) | Series temporales de $C_A$ (caso base) |
| ![Serie temporal C caso base](figuras_modelo_final/03_series_temporales_W_caso_base.png) | Series temporales de W (caso base) |
| ![Barrido hemoglobina](figuras_modelo_final/04_barrido_Hemoglobina.png) | Barrido de la hemoglobina |
| ![Validación Z0 cero](figuras_modelo_final/05_validacion_Z0_cero.png) | Validación frente a Modelo 4 |

# PARÁMETROS POR DEFECTO

| Parámetro | Valor | Descripción |
|---|---|---|
| `--N` | 160 | Nodos espaciales |
| `--T` | 4.0 | Tiempo final adimensional |
| `--dt` | 0.0025 | Paso temporal adimensional |
| `PeA` | 20 | Número de Péclet alveolar |
| `PeC` | 20 | Número de Péclet capilar |
| `DaA` | 3 | Número de Damköhler alveolar |
| `DaC` | 3 | Número de Damköhler capilar |
| `Z0` | 55.6 | Concentración adimensional de hemoglobina (fisiológico) |

# TRAZABILIDAD 

**Commit asociado a la memoria:** `1bf9ca2`

**DOI del Repositorio:** [![DOI](https://zenodo.org/badge/1235796656.svg)](https://doi.org/10.5281/zenodo.20691104)
