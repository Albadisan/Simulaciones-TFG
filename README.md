# Simulaciones-TFG
En el desarrollo de mi TFG "Modelo 1D de Reacción-Advección-Difusión para Intercambio Gaseoso en Biochips Pulmonares: Análisis y Extensiones Fisiológicas" se pretende desarrollar una simulación en Python. Se crea este repositorio para guardar los scripts de cada etapa hasta conseguir la simulación final.

La simulación que se pretende hacer es la del sistema:

$$ \frac{\partial\hat{C}_A}{\partial {\hat{t}}} + r_v\frac{\partial\hat{C}_A}{\partial \hat{x}} = \frac{1}{Pe_A}\frac{\partial^2\hat{C}_A}{\partial x^2} - Da_A(\hat{C}_A -\hat{W})$$

$$ \hat{\alpha}(\hat{W})\left[\frac{\partial\hat{W}}{\partial {\hat{t}}} + \frac{\partial\hat{W}}{\partial \hat{x}} \right] = \frac{1}{Pe_C}\frac{\partial}{\partial \hat{x}} \left(\hat{D}(\hat{W})\frac{\partial \hat{W}}{\partial \hat{x}} \right) + Da_C(\hat{C}_A -\hat{W})$$

Para lograrlo se hace en 4 etapas: se imponen 4 hipótesis al principio y en cada etapa se van relajando hipótesis hasta alcanzar las ecuaciones del modelo que pretendemos simular. De esta forma vamos validando poco a poco las simulaciones.

- MODELO 1 : Primera etapa para llegar a la simulación final. 
  
    Se imponen 4 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$
     3) Concentración alveolar fija: $\hat{C}_A(x) = hat{C}_A^0$
     4) No hay difusión axial: $Pe_C \to \inf $

    Se resuelve por tanto la ecuación:

     $$\frac{d \hat{W}}{dx} = DaC (\hat{C}_A^0 − \hat{W})$$
  
     $$\hat{W}(0) = \hat{W}_{in}$$

- MODELO 2 : Segunda etapa para llegar a la simulación final. 
  
    Se imponen 4 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$
     3) Concentración alveolar fija: $\hat{C}_A(x) = \hat{C}_A^0$
     4) Hay difusión axial: $Pe_C$ finito

    Se resuelve por tanto la ecuación:

     $$\frac{d \hat{W}}{dx}= \frac{1}{Pe_C} \frac{d^2 \hat{W}}{dx^2} + DaC (\hat{C}_A^0 − \hat{W})$$
  
     $$\hat{W}(0) = \hat{W}_{in}$$
 
- MODELO 3 : Tercera etapa para llegar a la simulación final. 
  
    Se imponen 2 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$

    Se resuelve por tanto el sistema:

     $$r_v \frac{d \hat{C}}{dx}= \frac{1}{Pe_A}\frac{d^2 \hat{C}}{dx^2} - DaA (\hat{C}_A − \hat{W})$$

     $$\frac{d \hat{W}}{dx}= \frac{1}{Pe_C} \frac{d^2 \hat{W}}{dx^2} + DaC (\hat{C}_A − \hat{W})$$

     En $x=0$: $\hspace{1cm} \hat{C}_A (0) = \hat{C}_{A,in}$, $\hat{W}(0) = \hat{W}_{in}$
  
     En $x=1$: $\hspace{1cm} \frac{d \hat{C}_A}{dx} = 0$, $\frac{d \hat{W}}{dx} = 0$ (Neumann homogénea)

- MODELO 4 : Cuarta etapa para llegar a la simulación final. 
  
    Se imponen 1 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$

    Se resuelve por tanto el sistema:

     $$\frac{\partial \hat{C}}{\partial t} + r_v \frac{\partial \hat{C}}{\partial x}= \frac{1}{Pe_A}\frac{\partial^2 \hat{C}}{\partial x^2} - DaA (\hat{C}_A − \hat{W})$$

     $$\frac{\partial \hat{W}}{\partial t} + \frac{\partial \hat{W}}{\partial x}= \frac{1}{Pe_C} \frac{\partial^2 \hat{W}}{\partial x^2} + DaC (\hat{C}_A − \hat{W})$$

     En $x=0$: $\hspace{1cm} \hat{C}_A (0) = \hat{C}_{A,in}$, $\hat{W}(0) = \hat{W}_{in}$
  
     En $x=1$: $\hspace{1cm} \frac{d\hat{C}_A}{dx} = 0, \frac{d\hat{W}}{dx} = 0$ (Neumann homogénea)
  
 - MODELO FINAL : simulación final 

    Se añaden las funciones no lineales y se resuelve por fin el sistema:

   $$ \frac{\partial\hat{C}_A}{\partial {\hat{t}}} + r_v\frac{\partial\hat{C}_A}{\partial \hat{x}} = \frac{1}{Pe_A}\frac{\partial^2\hat{C}_A}{\partial x^2} - Da_A(\hat{C}_A -\hat{W})4$$

$$ \hat{\alpha}(\hat{W})\left[\frac{\partial\hat{W}}{\partial {\hat{t}}} + \frac{\partial\hat{W}}{\partial \hat{x}} \right] = \frac{1}{Pe_C}\frac{\partial}{\partial \hat{x}} \left(\hat{D}(\hat{W})\frac{\partial \hat{W}}{\partial \hat{x}} \right) + Da_C(\hat{C}_A -\hat{W})$$

