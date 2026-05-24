
Repositorio que acompaña al Trabajo de Fin de Grado:
   ### Modelo 1D de Reacción-Advección-Difusión para Intercambio Gaseoso en Biochips Pulmonares: Análisis y Extensiones Fisiológicas

En el desarrollo del TFG se pretende simular en Python el sistema acoplado:

$$
\frac{\partial \hat{C}_{A}}{\partial \hat{t}} + r_{v}\ \frac{ \partial \hat{C}_A}{\partial \hat{x}}
= \frac{1}{\mathrm{Pe}_{A}}\ \frac{\partial^2 \hat{C}_{A}}{\partial \hat{x}^2} - \mathrm{Da}_A\ (\hat{C}_A - \hat{W})
$$

$$ 
 \hat{\alpha}(\hat{W})\left[\frac{\partial\hat{W}}{\partial {\hat{t}}} + \frac{\partial\hat{W}}{\partial \hat{x}} \right] = \frac{1}{\mathrm{Pe}_C}\frac{\partial}{\partial \hat{x}} \left(\hat{D}(\hat{W})\frac{\partial \hat{W}}{\partial \hat{x}} \right) + Da_C(\hat{C}_A -\hat{W})
$$

La estrategia es relajar hipótesis progresivamente en cuatro etapas, validando cada modelo antes de pasar al siguiente.

    Nota: Todas las variables son adimensionales; se omite el sombrero en lo sucesivo.

# Estructura del repositorio:
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


