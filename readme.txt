1.-Generar la configuración inicial
    compilar el programa conf.py:
    python3 conf.py

    Seleccionar la opcion 1 (FCC)
    Seleccionar el número de réplicas 5 Número de átomos = 500
    Introducir densidad numérica adimensional, dens=0.8442
    Seleccionar la opción 2 (Foto de la configuración inicial)
    Seleccionar la opción 0 (salir)
    Sugerencia de nombre del archivo de salida fort.dat

2.-Simulación Monte Carlo
    compilar el programa MC_python.py:
    python3 MC_python.py

    Lennard-Jones       Título de la corrida
    5000                Número de ciclos
    500                 Número de ciclos entre salida
    100                 Frecuencia rdf
    5                   Intervalo para actualizar el desplazamiento max.
    fort.dat            Nombre del archivo de la configuración inicial
    0.8442              Densidad
    0.71                Temperatura
    2.5                 Distancia de corte
