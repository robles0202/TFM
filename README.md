# TFM · Análisis de complejidad genómica en muestras de la ISS

Trabajo de Fin de Máster (Bioinformática) que analiza si el vuelo espacial
deja alguna huella detectable en la complejidad del ADN. Se comparan
muestras conservadas en la Estación Espacial Internacional (condición
**Flight**) frente a sus controles en tierra (condición **Ground**),
usando los datos públicos **NASA OSD-84 / GLDS-84** (secuenciador
portátil MinION).

## Qué calcula

A partir del perfil de k-mers de cada secuencia (FASTQ/FASTA), el script
calcula cuatro métricas de complejidad genómica y las compara siempre
frente a un genoma "de control" generado al azar con la misma composición
de bases, para separar la señal biológica real del simple efecto del
tamaño de la muestra:

- **Entropía de Shannon** (global, por k-mer, acumulada y normalizada)
- **K-mers ausentes** (recuento empírico y modelo teórico del
  coleccionista de cupones)
- **IEPWRMkmer** (distancias intra/inter-grupo basadas en posición y
  rareza de cada k-mer)
- **Biobit**
- **Genome Signature (GS)**

## Instalación

```bash
pip install -r requirements.txt
```

Requiere Python 3.9+.

## Uso

`TFM.py` no tiene interfaz de línea de comandos: las rutas de las
carpetas de datos (Ground/Flight) y la semilla de aleatoriedad se editan
directamente en el bloque `if __name__ == "__main__":` al final del
fichero. Después, simplemente:

```bash
python TFM.py
```

Esto ejecuta las cuatro métricas sobre los dos conjuntos de datos del
trabajo (metagenómico y *E. coli* consenso) y genera las figuras y
tablas CSV correspondientes en el directorio actual.

## Estructura del código

El fichero está dividido en bloques (marcados con comentarios `#===`):

1. Lectura de ficheros
2. Conteo de nucleótidos y frecuencias
3. Utilidades de k-mers
4. Entropía de Shannon
5. K-mers ausentes
6. IEPWRMkmer
7. Biobit
8. Genome Signature
9. Ejecución del análisis completo

## Resultado

Con estos datos, las cuatro métricas no detectan una diferencia clara
entre Ground y Flight: las diferencias observadas quedan siempre dentro
de la variabilidad normal entre muestras del mismo grupo.

## Autor

Rubén Ballester Robles

## Licencia

Este proyecto está bajo licencia MIT (ver [LICENSE](LICENSE)).
