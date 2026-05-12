# Simulaciones-TFG
En el desarrollo de mi TFG "Modelo 1D de Reacción-Advección-Difusión para Intercambio Gaseoso en Biochips Pulmonares: Análisis y Extensiones Fisiológicas" se pretende desarrollar una simulación en Python. Se crea este repositorio para guardar los scripts de cada etapa hasta conseguir la simulación final.

La simulación que se pretende hacer es la del sistema:

$$
\frac{\partial \hat{C}_{A}}{\partial \hat{t}} + r_{v}\ \frac{ \partial \hat{C}_A}{\partial \hat{x}}
= \frac{1}{\mathrm{Pe}_{A}}\ \frac{\partial^2 \hat{C}_{A}}{\partial \hat{x}^2} - \mathrm{Da}_A\ (\hat{C}_A - \hat{W})
$$

$$ 
 \hat{\alpha}(\hat{W})\left[\frac{\partial\hat{W}}{\partial {\hat{t}}} + \frac{\partial\hat{W}}{\partial \hat{x}} \right] = \frac{1}{\mathrm{Pe}_C}\frac{\partial}{\partial \hat{x}} \left(\hat{D}(\hat{W})\frac{\partial \hat{W}}{\partial \hat{x}} \right) + Da_C(\hat{C}_A -\hat{W})
$$

Para lograrlo se hace en 4 etapas: se imponen 4 hipótesis al principio y en cada etapa se van relajando hipótesis hasta alcanzar las ecuaciones del modelo que pretendemos simular. De esta forma vamos validando poco a poco las simulaciones.

Nota: Todas las variables son adimensionales; se omite el sombrero en lo sucesivo.

- MODELO 1 : Primera etapa para llegar a la simulación final. 
  
    Se imponen 4 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$
     3) Concentración alveolar fija: $C_A(x) = C_A^0$
     4) No hay difusión axial: $Pe_C \to \inf $

    Se resuelve por tanto la ecuación:

     $$\frac{d W}{dx} = \mathrm{Da}_C (C_A^0 − W)$$
  
     $$W(0) = W_{in}$$

- MODELO 2 : Segunda etapa para llegar a la simulación final. 
  
    Se imponen 4 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$
     3) Concentración alveolar fija: $C_A(x) = C_A^0$
     4) Hay difusión axial: $Pe_C$ finito

    Se resuelve por tanto la ecuación:

     $$\frac{d W}{dx} = \frac{1}{\mathrm{Pe}_{C}} \frac{d^2 W}{dx^2} + \mathrm{Da}_C (C_A^0 − W)$$
  
     En $x=0$: $W(0) = W_{in}$  (Dirichlet en la entrada)
  
     En $x=1$: $\frac{d W}{dx}(1) = 0$ (Neumann en la salida)
 
- MODELO 3 : Tercera etapa para llegar a la simulación final. 
  
    Se imponen 2 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$

    Se resuelve por tanto el sistema:

     $$r_v \frac{d C}{d x}= \frac{1}{\mathrm{Pe}_{A}}\frac{d^2 C}{\partial x^2} - \mathrm{Da}_A(C_A − W)$$
     
     $$\frac{d W}{dx}= \frac{1}{\mathrm{Pe}_{C}} \frac{d^2 W}{dx^2} + \mathrm{Da}_C (C_A^0 − W)$$

     En $x=0$: $\hspace{1cm} C_A (0) = C_{A,in}$, $W(0) = W_{in}$  (Dirichlet en la entrada)
  
     En $x=1$: $\hspace{1cm} \frac{d C_A}{dx}(1) = 0$, $\frac{d W}{dx}(1) = 0$ (Neumann en la salida)

- MODELO 4 : Cuarta etapa para llegar a la simulación final. 
  
    Se imponen 1 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$

    Se resuelve por tanto el sistema:

     $$\frac{\partial C}{\partial t} + r_v\frac{\partial C}{\partial x} = \frac{1}{\mathrm{Pe}_{A}}\frac{\partial^2 C}{\partial x^2} - \mathrm{Da}_A(C_A − W)$$

     $$\frac{\partial W}{\partial t} + \frac{\partial W}{\partial x} = \frac{1}{\mathrm{Pe}_{C}} \frac{\partial^2 W}{\partial x^2} + \mathrm{Da}_C (C_A^0 − W)$$

     Condición inicial : $C_A (x,0) = C_{A,0}(x)$
  
     En $x=0$: $\hspace{1cm} C_A (0,t) = C_{A,in}$, $W(0,t) = W_{in}$ (Dirichlet en la entrada)
  
     En $x=1$: $\hspace{1cm} \frac{dC_A}{dx}(1,t) = 0, \frac{dW}{dx}(1,t) = 0$ (Neumann en la salida)
  
 - MODELO FINAL : simulación final 

    Se añaden las funciones no lineales
   
     $$ \alpha(W) = 1 + 4 Z_0f'(W)}$$

     $$
       D(W) = 1 + 4 \delta Z_0 f'(W)
     $$
   
    Y se tiene el sistema original que pretendíamos simular.
