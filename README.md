# Simulaciones-TFG
En el desarrollo de mi TFG "Modelo 1D de Reacción-Advección-Difusión para Intercambio Gaseoso en Biochips Pulmonares: Análisis y Extensiones Fisiológicas" se pretende desarrollar una simulación en Python. Se crea este repositorio para guardar los scripts de cada etapa hasta conseguir la simulación final.

- MODELO 1 : Primera etapa para llegar a la simulación final. 
  
    Se imponen 3 hipótesis:

        1) Sin hemoglobina: Z_0 = 0
        2) Estado estacionario: d/dt = 0
        3) Concentración alveolar fija: C_A(x)=C_A^0

    Se resuelve por tanto la ecuación:

      dW/dx= DaC ( ˆC_A^0 − W),
      W(0) = ˆWin.
