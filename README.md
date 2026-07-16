# Análisis comparativo de la complejidad genómica en condiciones de espacio y Tierra

Código y materiales del Trabajo de Fin de Máster del **Máster Universitario en
Bioinformática** de la Universitat de València.

- **Autor:** Rubén Ballester Robles
- **Tutores:** Juli Peretó Magraner (Universitat de València) y Javier Buceta Fernández (i2SysBio)

## Descripción

Este repositorio contiene el código empleado para evaluar si el almacenamiento y
la secuenciación a bordo de la Estación Espacial Internacional (ISS) dejan una
huella detectable en la **complejidad genómica** del ADN. El análisis compara
muestras procesadas en la ISS (*Flight*) y en tierra (*Ground*) del estudio
[OSD-84](https://osdr.nasa.gov/bio/repo/data/studies/OSD-84) de la NASA, en dos
niveles: las lecturas metagenómicas completas y los genomas consenso de
*Escherichia coli* K-12 MG1655.

Se calculan cuatro medidas de complejidad basadas en *k*-mers —entropía de
Shannon, firma genómica (*Genome Signature*), Biobit e IEPWRMkmer— junto con
**controles aleatorios pareados por longitud y composición**, y se generan las
figuras y tablas del trabajo.

## Contenido

- `TFM.py` — script principal. Incluye la lectura de FASTQ/FASTA, el
  conteo y las distribuciones de *k*-mers, las cuatro métricas de complejidad, la
  generación de los controles aleatorios y la exportación de figuras y tablas. Las
  rutas de los datos y la semilla de aleatoriedad se configuran editando el bloque
  `if __name__ == "__main__":` situado al final del fichero.
- `requirements.txt` — dependencias necesarias para ejecutar los análisis.
- `LICENSE` — licencia del repositorio.

El preprocesamiento previo —alineamiento con **minimap2**, detección de variantes
y generación de los genomas consenso con **bcftools**, y ensamblaje con **Flye**—
se realizó con herramientas externas. Este repositorio cubre el análisis de
complejidad y la visualización.

## Disponibilidad de datos

Siguiendo la práctica habitual, el repositorio contiene **el código**, no los
datos voluminosos. GitHub no está pensado para datos de gran tamaño y bloquea
ficheros de más de 100 MB. Los datos se organizan en tres niveles:

**1. Datos brutos de secuenciación (públicos).** Las lecturas de nanoporos
(Oxford Nanopore MinION) del estudio OSD-84 no se incluyen aquí por su tamaño y
por estar ya disponibles públicamente:

- Open Science Data Repository (OSDR) de la NASA, estudio OSD-84 / GLDS-84:
  <https://osdr.nasa.gov/bio/repo/data/studies/OSD-84>

**2. Genomas de referencia (públicos).** Descargables del NCBI:

- *Escherichia coli* K-12 MG1655 (NC_000913.3):
  <https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3>

**3. Datos derivados generados en este trabajo.** Los seis genomas consenso de
*Escherichia coli* K-12 MG1655 constituyen la entrada directa del análisis de
complejidad a nivel de genoma individual. A partir de estos ficheros FASTA puede
reproducirse el análisis de complejidad ejecutando `TFM.py`, sin
necesidad de repetir el alineamiento, el *variant calling* ni la generación de los
consensos.

Los genomas consenso se han archivado de forma permanente y citable en Zenodo:

- Zenodo: <https://doi.org/10.5281/zenodo.21242575>
- DOI: `10.5281/zenodo.21242575`

## Obtención de los genomas consenso de *E. coli*

Los genomas consenso de *Escherichia coli* K-12 MG1655 se obtuvieron a partir de
las muestras metagenómicas seleccionadas mediante un flujo de trabajo externo al
script `TFM.py`. Este procedimiento se describe en la sección de
Materiales y Métodos del TFM y consta de cuatro etapas principales: alineamiento
contra el genoma de referencia, indexado de los archivos BAM, detección de
variantes y generación de los genomas consenso.

### 1. Alineamiento contra el genoma de referencia

Cada archivo FASTQ se alineó de forma independiente contra el genoma de referencia
de *Escherichia coli* K-12 MG1655 (`NC_000913.3`) mediante `minimap2`, usando el
preset específico para lecturas Oxford Nanopore (`map-ont`). La salida SAM se
ordenó y convirtió directamente a BAM con `samtools sort`.

```bash
minimap2 -ax map-ont ecoli_K12_MG1655.fasta [muestra].fastq \
  | samtools sort -o [muestra]_ecoli.bam
```

El parámetro `-ax map-ont` indica que se trata de lecturas largas de Oxford
Nanopore y genera una salida en formato SAM adecuada para su procesamiento con
`samtools`.

### 2. Indexado de los archivos BAM

Los archivos BAM ordenados se indexaron con `samtools index` para permitir el
acceso rápido por coordenadas durante los pasos posteriores.

```bash
samtools index [muestra]_ecoli.bam
```

### 3. Detección de variantes (*variant calling*)

La detección de variantes se realizó con `bcftools`. En cada comparación se
analizaron conjuntamente la muestra `Ground` y la muestra `Flight`
correspondientes a la misma fecha experimental. Primero, `bcftools mpileup`
evaluó las lecturas alineadas sobre cada posición del genoma de referencia;
después, `bcftools call` identificó SNPs e INDELs y generó un archivo VCF
comprimido.

```bash
bcftools mpileup -f ecoli_K12_MG1655.fasta ground_ecoli.bam flight_ecoli.bam \
  | bcftools call -mv -Oz -o variants.vcf.gz
```

El archivo resultante `variants.vcf.gz` contiene las variantes detectadas respecto
al genoma de referencia y permite distinguir en qué muestra está presente cada
variante.

### 4. Generación de los genomas consenso

Finalmente, se aplicaron las variantes detectadas al genoma de referencia mediante
`bcftools consensus`, generando un archivo FASTA consenso para cada muestra.

```bash
bcftools consensus -f ecoli_K12_MG1655.fasta \
  -s [muestra] variants_[fecha].vcf.gz > consenso_[muestra].fasta
```

El archivo `consenso_[muestra].fasta` representa el genoma consenso de cada
muestra, es decir, una secuencia basada en la referencia de *E. coli* K-12 MG1655
pero modificada con los SNPs e INDELs identificados en esa muestra concreta.

De este modo se obtuvieron seis genomas consenso: tres correspondientes a la
condición `Ground` y tres correspondientes a la condición `Flight`. Estos seis
FASTA son los datos derivados depositados en Zenodo y utilizados como entrada para
los análisis de Genome Signature, Biobit, IEPWRMkmer y entropía de Shannon.

## Requisitos

Python 3.11 o superior.

Instalación de dependencias:

```bash
pip install -r requirements.txt
```

Principales dependencias:

- `numpy`
- `pandas`
- `matplotlib`
- `biopython`
- `scipy`

Además, para reproducir desde cero la obtención de los genomas consenso se
requieren herramientas externas:

- `minimap2`
- `samtools`
- `bcftools`
- `flye`

Estas herramientas no son necesarias si se parte directamente de los genomas
consenso depositados en Zenodo.

Para la reproducibilidad del código, hice uso de la semilla `356`.

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

## Cómo citar

Si utilizas este código, cita el Trabajo de Fin de Máster:

> Ballester Robles, R. (2026). *Análisis comparativo de la complejidad genómica
> en condiciones de espacio y Tierra*. Trabajo de Fin de Máster, Máster
> Universitario en Bioinformática, Universitat de València.

Si utilizas los genomas consenso de *E. coli* generados en este trabajo, cita
también el depósito de Zenodo:

> Ballester Robles, R. (2026). *Genomas consenso de Escherichia coli K-12 MG1655
> derivados del estudio OSD-84*. Zenodo.
> <https://doi.org/10.5281/zenodo.21242575>

## Licencia

Distribuido bajo la licencia MIT. Consulta el fichero [`LICENSE`](LICENSE).

