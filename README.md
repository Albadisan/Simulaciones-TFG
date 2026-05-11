# Simulaciones-TFG
En el desarrollo de mi TFG "Modelo 1D de Reacción-Advección-Difusión para Intercambio Gaseoso en Biochips Pulmonares: Análisis y Extensiones Fisiológicas" se pretende desarrollar una simulación en Python. Se crea este repositorio para guardar los scripts de cada etapa hasta conseguir la simulación final.

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
     4) No hay difusión axial: $Pe_C$ finito

    Se resuelve por tanto la ecuación:

     $$\frac{d \hat{W}}{dx}= \frac{1}{Pe_C} \frac{d^2 \hat{W}}{dx^2} + DaC (\hat{C}_A^0 − \hat{W})$$
  
     $$\hat{W}(0) = \hat{W}_{in}$$
 
- MODELO 3 : Tercra etapa para llegar a la simulación final. 
  
    Se imponen 2 hipótesis:

     1) Sin hemoglobina: $Z_0 = 0$
     2) Estado estacionario: $d/dt = 0$

    Se resuelve por tanto el sistema:

     $$r_v \frac{d \hat{C}}{dx}= (1/Pe_A)\frac{d^2 \hat{W}}{dx^2} - DaA (\hat{C}_A^0 − \hat{W})$$  

     $$\frac{d \hat{W}}{dx}= \frac{1}{Pe_C} \frac{d^2 \hat{W}}{dx^2} + DaC (\hat{C}_A^0 − \hat{W})$$
  
     $$\hat{W}(0) = \hat{W}_{in}$$
