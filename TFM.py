#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TFM - Análisis de complejidad genómica de las muestras de ADN de la ISS
(datos NASA OSD-84, secuenciador MinION), comparando la condición Ground
(tierra) con la de Flight (vuelo).

Se leen los FASTQ/FASTA de las dos carpetas y se calculan sobre el perfil
de k-mers las 4 métricas del trabajo: entropía de Shannon, k-mers ausentes,
IEPWRMkmer, Biobit y Genome Signature. Para separar la señal biológica del
ruido de fondo, cada métrica se compara además con un genoma "de control"
generado al azar con la misma composición de bases.

Bloques del fichero (busca los comentarios con "====="):
  1. Lectura de ficheros
  2. Conteo de nucleótidos y frecuencias
  3. Utilidades de k-mers
  4. Entropía de Shannon
  5. K-mers ausentes
  6. IEPWRMkmer
  7. Biobit
  8. Genome Signature
  9. Ejecución del análisis completo

Necesita: biopython, numpy, pandas, matplotlib, scipy.

Rubén Ballester Robles - TFM Bioinformática
"""

from collections import defaultdict
import csv
import itertools
import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d import Axes3D  # sin esto no funciona projection="3d"
from Bio import SeqIO
from scipy.spatial.distance import cityblock

BASES = "ACGT"
CONJUNTO_BASES = set(BASES)


# ============================================================================
# 1. LECTURA DE FICHEROS
# ============================================================================

def FASTQ(filepath):
    """
    Lee un fichero FASTQ y devuelve todas sus lecturas concatenadas en una
    única cadena de texto.

    Parameters
    ----------
    filepath : str
        Ruta al fichero FASTQ que se desea leer.

    Returns
    -------
    str
        Cadena con todas las secuencias del fichero concatenadas, sin saltos
        de línea ni espacios.

    Raises
    ------
    FileNotFoundError
        Si la ruta indicada no existe.
    """

    secuencias = []

    with open(filepath, "r") as f:
        for numero_linea, linea in enumerate(f, start=1):
            if numero_linea % 4 == 2:
                secuencias.append(linea.strip())

    return "".join(secuencias)


def FASTA(filepath):
    """
    Lee un fichero FASTA y devuelve todas sus secuencias concatenadas en una
    única cadena de texto.

    Parameters
    ----------
    filepath : str
        Ruta al fichero FASTA que se desea leer.

    Returns
    -------
    str
        Cadena con todas las secuencias del fichero concatenadas y pasadas a
        mayúsculas.

    Raises
    ------
    FileNotFoundError
        Si la ruta indicada no existe.
    """

    secuencias = []

    for record in SeqIO.parse(filepath, "fasta"):
        secuencias.append(str(record.seq).upper())

    return "".join(secuencias)


def cargar_carpeta(carpeta):
    """
    Lee todos los ficheros FASTQ y FASTA de una carpeta y devuelve la lista de 
    secuencias y la lista de nombres de fichero en el mismo orden.

    Parameters
    ----------
    carpeta : str
        Ruta a la carpeta que contiene los ficheros de secuencias.

    Returns
    -------
    tuple of (list of str, list of str)
        - secuencias : lista con la secuencia de ADN de cada fichero.
        - nombres    : lista con el nombre de cada fichero.

    Raises
    ------
    FileNotFoundError
        Si la carpeta no existe.
    """

    secuencias = []
    nombres = []

    for archivo in os.listdir(carpeta):
        ruta = os.path.join(carpeta, archivo)
        if archivo.endswith(".fastq"):
            secuencias.append(FASTQ(ruta))
            nombres.append(archivo)
        elif archivo.endswith(".fasta"):
            secuencias.append(FASTA(ruta))
            nombres.append(archivo)
            
    return secuencias, nombres


# ============================================================================
# 2. CONTEO DE NUCLEÓTIDOS, FRECUENCIAS Y SECUENCIAS ALEATORIAS
# ============================================================================

def conteo(seq):
    """
    Cuenta cuántas veces aparece cada mono-, di-, tri-, tetra- y pentanucleótido
    en una secuencia al recorrer la secuencia con ventanas solapantes.

    Parameters
    ----------
    seq : str
        Secuencia de ADN.

    Returns
    -------
    tuple of dict
        Cinco diccionarios {subcadena: número de apariciones}, uno para cada
        tamaño de ventana (mono, di, tri, tetra y penta nucleótidos),
        en ese orden.
    """

    mono = defaultdict(int)
    di = defaultdict(int)
    tri = defaultdict(int)
    tetra = defaultdict(int)
    penta = defaultdict(int)

    for base in seq:
        mono[base] += 1
    for i in range(len(seq) - 1):
        di[seq[i:i + 2]] += 1
    for i in range(len(seq) - 2):
        tri[seq[i:i + 3]] += 1
    for i in range(len(seq) - 3):
        tetra[seq[i:i + 4]] += 1
    for i in range(len(seq) - 4):
        penta[seq[i:i + 5]] += 1

    return mono, di, tri, tetra, penta


def frecuencia(mono, seq):
    """
    Calcula la frecuencia relativa de cada nucleótido en una secuencia (número
    de apariciones dividido por la longitud de la secuencia).

    Parameters
    ----------
    mono : dict
        Diccionario {nucleótido: número de apariciones}.
    seq : str
        Secuencia de ADN de la que proceden los conteos.

    Returns
    -------
    dict
        Diccionario {nucleótido: frecuencia relativa}, con valores entre 0 y 1.

    Raises
    ------
    ZeroDivisionError
        Si seq está vacía.
    """

    longitud_seq = len(seq)
    return {base: apariciones / longitud_seq for base, apariciones in mono.items()}


def aleatorio(freq, longitud):
    """
    Genera una secuencia de ADN al azar que tiene la misma composición de bases
    proporcionada por el usuario. Primero coloca el número exacto de cada base 
    que corresponde a su frecuencia, completa los huecos por redondeo eligiendo 
    bases al azar según esas mismas frecuencias y baraja todo el conjunto.

    Parameters
    ----------
    freq : dict
        Diccionario {nucleótido: frecuencia relativa}.
    longitud : int
        Número de nucleótidos que tendrá la secuencia generada.

    Returns
    -------
    str
        Secuencia de ADN aleatoria de la longitud solicitada. Para que el
        resultado sea reproducible hay que fijar antes una semilla con
        random.seed(...).
    """

    longitud = int(longitud)
    secuencia = []

    total = 0
    for nucleotido, f in freq.items():
        cantidad = int(f * longitud)
        secuencia.extend([nucleotido] * cantidad)
        total += cantidad

    # El redondeo a entero puede dejar bases por asignar; se completan al azar.
    faltantes = longitud - total
    if faltantes > 0:
        nucleotidos = list(freq.keys())
        pesos = [freq[n] for n in nucleotidos]
        secuencia.extend(random.choices(nucleotidos, weights=pesos, k=faltantes))

    random.shuffle(secuencia)
    return "".join(secuencia)


def calcular_frecuencias_carpeta(carpeta):
    """
    Calcula la frecuencia de nucleótidos de cada fichero de una carpeta y la
    frecuencia media del conjunto.

    Parameters
    ----------
    carpeta : str
        Ruta a la carpeta con los ficheros de secuencias.

    Returns
    -------
    tuple of (list of dict, dict)
        - frecuencias_archivos : lista con el diccionario de frecuencias de
          cada fichero.
        - media_freq : diccionario {nucleótido: frecuencia media} del grupo,
          con las cuatro bases A, C, G y T.

    Raises
    ------
    FileNotFoundError
        Si la carpeta no existe.
    """

    secuencias, _ = cargar_carpeta(carpeta)

    frecuencias_archivos = []
    for seq in secuencias:
        mono, _, _, _, _ = conteo(seq)
        frecuencias_archivos.append(frecuencia(mono, seq))

    media_freq = {}
    for base in BASES:
        valores = [f.get(base, 0) for f in frecuencias_archivos]
        media_freq[base] = np.mean(valores)

    return frecuencias_archivos, media_freq


def generar_genomas_aleatorios(carpeta_ground, carpeta_flight, n_genomas=3,
                               muestra=False, semilla=None):
    """
    Genera genomas aleatorios de control a partir de la composición media de
    las muestras Ground y Flight (media entre ambos grupos).

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    n_genomas : int, optional
        Número de réplicas a generar cuando muestra=False (por defecto 3).
    muestra : bool, optional
        Si es True, genera un genoma aleatorio por cada muestra real, con su
        misma longitud. Si es False (por defecto), genera `n_genomas`
        réplicas con la longitud de la primera muestra de Ground.
    semilla : int or None, optional
        Semilla para la generación reproducible. Por defecto None.

    Returns
    -------
    tuple
        Si muestra=True: (genomas_ground, genomas_flight, nombres_ground,
        nombres_flight, freq_aleatorio). Si muestra=False: (lista_genomas,
        freq_aleatorio).

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    IndexError
        Si carpeta_ground no contiene ninguna secuencia (modo muestra=False).
    """

    if semilla is not None:
        random.seed(semilla)

    _, freq_ground = calcular_frecuencias_carpeta(carpeta_ground)
    _, freq_flight = calcular_frecuencias_carpeta(carpeta_flight)

    freq_aleatorio = {base: (freq_ground[base] + freq_flight[base]) / 2
                      for base in BASES}

    secuencias_ground, nombres_ground = cargar_carpeta(carpeta_ground)
    secuencias_flight, nombres_flight = cargar_carpeta(carpeta_flight)

    if muestra:
        print("Generando genomas aleatorios (uno por muestra real):")
        print("  Ground:")
        genomas_ground = []
        for nombre, seq in zip(nombres_ground, secuencias_ground):
            print(f"    {nombre}: {len(seq):,} bp")
            genomas_ground.append(aleatorio(freq_aleatorio, len(seq)))

        print("  Flight:")
        genomas_flight = []
        for nombre, seq in zip(nombres_flight, secuencias_flight):
            print(f"    {nombre}: {len(seq):,} bp")
            genomas_flight.append(aleatorio(freq_aleatorio, len(seq)))

        return (genomas_ground, genomas_flight,
                nombres_ground, nombres_flight, freq_aleatorio)

    longitud_ref = len(secuencias_ground[0])
    print(f"Generando {n_genomas} genomas aleatorios "
          f"(longitud de referencia: {longitud_ref:,} bp):")
    lista_genomas = []
    for i in range(n_genomas):
        lista_genomas.append(aleatorio(freq_aleatorio, longitud_ref))
        print(f"  Aleatorio {i + 1}: {longitud_ref:,} bp")

    return lista_genomas, freq_aleatorio


def aleatorio2(carpeta_ground, carpeta_flight, semilla=None):
    """
    Genera un genoma aleatorio por condición usando la frecuencia media de
    nucleótidos de cada carpeta.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    semilla : int or None, optional
        Semilla para la generación reproducible. Por defecto None.

    Returns
    -------
    tuple of str
        (genoma_ground, genoma_flight), un genoma aleatorio por condición.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    if semilla is not None:
        random.seed(semilla)

    def procesar_carpeta(carpeta):
        # saca la frecuencia media de cada nucleótido y la longitud de referencia
        secuencias, nombres_archivos = cargar_carpeta(carpeta)
        frecuencias_archivos = []
        for nombre, seq in zip(nombres_archivos, secuencias):
            mono, _, _, _, _ = conteo(seq)
            freq = frecuencia(mono, seq)
            frecuencias_archivos.append(freq)
            print(f"\nFrecuencias archivo {nombre}: {dict(freq)}")

        media_freq = {}
        for nuc in BASES:
            media_freq[nuc] = np.mean([f.get(nuc, 0) for f in frecuencias_archivos])

        print(f"\nFrecuencias medias ({os.path.basename(carpeta.strip('/'))}): {media_freq}")
        return media_freq, len(secuencias[0])

    media_ground, longitud_ground = procesar_carpeta(carpeta_ground)
    media_flight, longitud_flight = procesar_carpeta(carpeta_flight)

    genoma_ground = aleatorio(media_ground, longitud_ground)
    genoma_flight = aleatorio(media_flight, longitud_flight)

    print("\nGenoma aleatorio Ground (primeros 200 nucleótidos):")
    print(genoma_ground[:200])
    print("\nGenoma aleatorio Flight (primeros 200 nucleótidos):")
    print(genoma_flight[:200])

    return genoma_ground, genoma_flight


def longitud(carpeta):
    """
    Devuelve la longitud de la primera secuencia de una carpeta.

    Parameters
    ----------
    carpeta : str
        Ruta a la carpeta con los ficheros de secuencias (.fastq o .fasta).

    Returns
    -------
    int
        Número de bases de la primera secuencia cargada de la carpeta.

    Raises
    ------
    IndexError
        Si la carpeta no contiene ninguna secuencia.
    """

    secuencias, _ = cargar_carpeta(carpeta)
    return len(secuencias[0])


# ============================================================================
# 3. UTILIDADES GENÉRICAS DE K-MERS
# ============================================================================

def generar_todos_kmers(k):
    """
    Genera la lista de todos los k-mers posibles de longitud k, ordenados
    alfabéticamente. Ojo, para k grandes esto crece muy rápido (4**k), con
    k=12 ya son más de 16 millones.

    Parameters
    ----------
    k : int
        Longitud de los k-mers.

    Returns
    -------
    tuple of (list of str, numpy.ndarray)
        - todos_kmers : lista con los 4**k k-mers posibles ordenada
          alfabéticamente.
        - n : array de índices.
    """

    todos_kmers = ["".join(p) for p in sorted(itertools.product(BASES, repeat=k))]
    n = np.arange(1, len(todos_kmers) + 1)
    return todos_kmers, n


def kmers(seq, k):
    """
    Cuenta cuántas veces aparece cada k-mer en una secuencia.

    Desliza una ventana de longitud k a lo largo de la secuencia, de un
    nucleótido en un nucleótido.

    Parameters
    ----------
    seq : str
        Secuencia de ADN a analizar.
    k : int
        Longitud de los k-mers.

    Returns
    -------
    collections.defaultdict
        Diccionario {k-mer: número de apariciones}.
    """

    conteos = defaultdict(int)
    for i in range(len(seq) - k + 1):
        conteos[seq[i:i + k]] += 1
    return conteos


def prob_kmer(seq, k):
    """
    Calcula la probabilidad (frecuencia relativa) de aparición de cada k-mer en
    una secuencia.

    Parameters
    ----------
    seq : str
        Secuencia de ADN a analizar.
    k : int
        Longitud de los k-mers.

    Returns
    -------
    dict
        Diccionario {k-mer: probabilidad de aparición}. Las probabilidades
        suman 1.

    Raises
    ------
    ZeroDivisionError
        Si seq es más corta que k (no se encuentra ningún k-mer).
    """

    conteos = kmers(seq, k)
    total = sum(conteos.values())
    return {kmer: count / total for kmer, count in conteos.items()}


def calcular_probabilidades_teoricas(todos_kmers):
    """
    Calcula la probabilidad teórica de cada k-mer bajo la hipótesis de azar
    puro (bases equiprobables e independientes, frecuencia 0.25 cada una).

    Parameters
    ----------
    todos_kmers : list of str
        Lista de k-mers (todos de la misma longitud), normalmente la salida de
        generar_todos_kmers().

    Returns
    -------
    numpy.ndarray
        Array con la probabilidad teórica de cada k-mer, igual a 0.25 elevado
        a la longitud del k-mer.
    """

    freq = {base: 0.25 for base in BASES}
    probs_teoricos = []
    for kmer in todos_kmers:
        p = 1.0
        for base in kmer:
            p *= freq[base]
        probs_teoricos.append(p)
    return np.array(probs_teoricos)


def calcular_medias_stds_kmers(lista_probs, todos_kmers):
    """
    Calcula la media y la desviación estándar de la probabilidad de cada k-mer
    a través de varias muestras (los k-mers ausentes cuentan como 0).

    Parameters
    ----------
    lista_probs : list of dict
        Lista de diccionarios {k-mer: probabilidad}, uno por muestra.
    todos_kmers : list of str
        Lista de k-mers para los que se calculan las estadísticas (define el
        orden del resultado).

    Returns
    -------
    tuple of (numpy.ndarray, numpy.ndarray)
        - medias : array con la probabilidad media de cada k-mer.
        - stds : array con la desviación estándar de cada k-mer.
    """

    medias = []
    stds = []
    for kmer in todos_kmers:
        muestras = [d.get(kmer, 0) for d in lista_probs]
        medias.append(np.mean(muestras))
        stds.append(np.std(muestras))
    return np.array(medias), np.array(stds)


def analizar_kmers_carpeta(carpeta, k):
    """
    Realiza de una sola vez los cálculos de k-mers comunes a una carpeta: carga
    las secuencias, sus perfiles de probabilidad, la lista de todos los k-mers
    posibles y las medias/desviaciones por k-mer.

    Parameters
    ----------
    carpeta : str
        Ruta a la carpeta con los ficheros de secuencias.
    k : int
        Longitud de los k-mers.

    Returns
    -------
    tuple
        (secuencias, lista_probs, todos_kmers, medias, stds, n, nombres_archivos)

    Raises
    ------
    FileNotFoundError
        Si la carpeta no existe.
    """

    secuencias, nombres_archivos = cargar_carpeta(carpeta)
    lista_probs = [prob_kmer(seq, k) for seq in secuencias]
    todos_kmers, n = generar_todos_kmers(k)
    medias, stds = calcular_medias_stds_kmers(lista_probs, todos_kmers)
    return secuencias, lista_probs, todos_kmers, medias, stds, n, nombres_archivos


def probabilidades_kmer_por_grupo(carpeta_ground, carpeta_flight, k, semilla=None):
    """
    Prepara los perfiles de k-mers de los tres grupos (Ground, Flight y
    control Aleatorio) para un valor de k. Se usa desde varias de las
    funciones de comparación para no repetir siempre el mismo bloque.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra.
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo.
    k : int
        Longitud de los k-mers.
    semilla : int or None, optional
        Semilla para la generación reproducible de los genomas aleatorios. Por
        defecto None.

    Returns
    -------
    tuple
        (lista_probs_ground, lista_probs_flight, lista_probs_random,
        todos_kmers, n)

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    secuencias_ground, _ = cargar_carpeta(carpeta_ground)
    secuencias_flight, _ = cargar_carpeta(carpeta_flight)

    lista_probs_ground = [prob_kmer(seq, k) for seq in secuencias_ground]
    lista_probs_flight = [prob_kmer(seq, k) for seq in secuencias_flight]

    genomas_g, genomas_f, _, _, _ = generar_genomas_aleatorios(
        carpeta_ground, carpeta_flight, muestra=True, semilla=semilla)
    lista_probs_random = [prob_kmer(g, k) for g in genomas_g + genomas_f]

    todos_kmers, n = generar_todos_kmers(k)
    return lista_probs_ground, lista_probs_flight, lista_probs_random, todos_kmers, n


def exportar_tabla_kmers(k, todos_kmers, datos_dict, nombre_base):
    """
    Exporta a un fichero CSV una tabla con la información de todos los k-mers.

    Parameters
    ----------
    k : int
        Longitud de los k-mers (se usa en el nombre del fichero).
    todos_kmers : list of str
        Lista de todos los k-mers, ordenada.
    datos_dict : dict
        {'NombreGrupo': {'media': array, 'std': array}, ...}.
    nombre_base : str
        Prefijo del nombre del fichero CSV.

    Returns
    -------
    None
        El resultado es el fichero "{nombre_base}_k{k}.csv" en el directorio
        actual.

    Raises
    ------
    OSError
        Si el fichero no se puede escribir.
    """

    data = []
    for idx, kmer in enumerate(todos_kmers, start=1):
        fila = {"Índice": idx, "k-mer": kmer}
        for nombre, datos in datos_dict.items():
            if "std" in datos:
                fila[f"{nombre}_media"] = datos["media"][idx - 1]
                fila[f"{nombre}_std"] = datos["std"][idx - 1]
            else:
                fila[nombre] = datos["media"][idx - 1]
        data.append(fila)

    nombre_archivo = f"{nombre_base}_k{k}.csv"
    pd.DataFrame(data).to_csv(nombre_archivo, index=False)
    print(f"\nTabla exportada correctamente como '{nombre_archivo}'.")


def exportar_tabla_kmers_filtrada(k, todos_kmers, datos_dict, nombre_base):
    """
    Exporta a CSV solo los k-mers ausentes (probabilidad 0) en al menos uno de
    los grupos Ground, Flight o Aleatorio.

    Parameters
    ----------
    k : int
        Longitud de los k-mers (se usa en el nombre del fichero).
    todos_kmers : list of str
        Lista de todos los k-mers, ordenada.
    datos_dict : dict
        {'Ground': {'media': array, 'std': array}, 'Flight': {...}, 'Aleatorio': {...}, ...}.
    nombre_base : str
        Prefijo del nombre del fichero CSV.

    Returns
    -------
    None
        El resultado es el fichero "{nombre_base}_k{k}_filtrado.csv" y un
        resumen del recuento por categorías impreso por pantalla. Se escribe
        fila a fila para no llenar la memoria cuando k es grande.

    Raises
    ------
    OSError
        Si el fichero no se puede escribir.
    """

    nombre_archivo = f"{nombre_base}_k{k}_filtrado.csv"

    contador_ground_cero = 0
    contador_flight_cero = 0
    contador_aleatorio_cero = 0
    contador_ambos_cero = 0
    total_kmers = len(todos_kmers)

    print(f"\nProcesando {total_kmers:,} k-mers...")

    with open(nombre_archivo, "w", newline="") as csvfile:
        fieldnames = ["Índice", "k-mer"]
        for nombre, datos in datos_dict.items():
            if "std" in datos:
                fieldnames.append(f"{nombre}_media")
                fieldnames.append(f"{nombre}_std")
            else:
                fieldnames.append(nombre)

        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for idx in range(total_kmers):
            ground_val = datos_dict["Ground"]["media"][idx]
            flight_val = datos_dict["Flight"]["media"][idx]
            aleatorio_val = datos_dict["Aleatorio"]["media"][idx]

            if ground_val == 0 or flight_val == 0 or aleatorio_val == 0:
                fila = {"Índice": idx + 1, "k-mer": todos_kmers[idx]}
                for nombre, datos in datos_dict.items():
                    if "std" in datos:
                        fila[f"{nombre}_media"] = datos["media"][idx]
                        fila[f"{nombre}_std"] = datos["std"][idx]
                    else:
                        fila[nombre] = datos["media"][idx]
                writer.writerow(fila)

                if ground_val == 0 and flight_val == 0:
                    contador_ambos_cero += 1
                elif ground_val == 0:
                    contador_ground_cero += 1
                elif flight_val == 0:
                    contador_flight_cero += 1
                if aleatorio_val == 0:
                    contador_aleatorio_cero += 1

            if (idx + 1) % 100_000_000 == 0:
                print(f"  Procesados {idx + 1:,} / {total_kmers:,} k-mers...")

    total_exportados = (contador_ground_cero + contador_flight_cero +
                        contador_ambos_cero + contador_aleatorio_cero)
    if total_exportados > 0:
        print("\nTABLA EXPORTADA")
        print(f"Archivo: '{nombre_archivo}'")
        print(f"  K-mers solo en Flight (Ground=0):     {contador_ground_cero:,}")
        print(f"  K-mers solo en Ground (Flight=0):     {contador_flight_cero:,}")
        print(f"  K-mers en ninguno (ambos=0):          {contador_ambos_cero:,}")
        print(f"  K-mers con Aleatorio=0:               {contador_aleatorio_cero:,}")
        print(f"  Total exportados:                     {total_exportados:,}")


def consultar_kmer_interactivo(k, todos_kmers, datos_dict):
    """
    Consulta de forma interactiva la información asociada a k-mers concretos.

    Pide un k-mer por teclado hasta que el usuario pulsa Enter sin escribir
    nada.

    Parameters
    ----------
    k : int
        Longitud de los k-mers.
    todos_kmers : list of str
        Lista de todos los k-mers posibles (define el índice de cada k-mer).
    datos_dict : dict
        {nombre_grupo: {'media': array, 'std': array}}.

    Returns
    -------
    None
        Función interactiva: imprime la información por pantalla.
    """

    while True:
        kmer_usuario = input(
            f"\nIntroduce un k-mer de longitud {k} para ver su información "
            f"(o pulsa Enter para saltar): ").strip().upper()

        if kmer_usuario == "":
            break

        if kmer_usuario in todos_kmers:
            idx = todos_kmers.index(kmer_usuario)
            print(f"\nInformación para el k-mer '{kmer_usuario}':")
            print(f"  Índice: {idx + 1}")
            for nombre, datos in datos_dict.items():
                if "std" in datos:
                    print(f"  {nombre}: media={datos['media'][idx]:.6f}, "
                          f"std={datos['std'][idx]:.6f}")
                else:
                    print(f"  {nombre}: {datos['media'][idx]:.6f}")
        else:
            print(f"El k-mer '{kmer_usuario}' no es válido para k={k}.")


def graficar_kmer(carpeta, k, etiqueta):
    """
    Representa la probabilidad media de cada k-mer de una carpeta, con barras
    de error.

    Parameters
    ----------
    carpeta : str
        Ruta a la carpeta con los ficheros de secuencias (.fastq o .fasta).
    k : int
        Longitud de los k-mers.
    etiqueta : str
        Texto identificativo de la condición ("Flight" o "Ground"), usado en
        el título y en la leyenda.

    Returns
    -------
    None
        Muestra la figura por pantalla con plt.show().
    """

    _, _, _, medias, stds, n, _ = analizar_kmers_carpeta(carpeta, k)

    plt.figure(figsize=(14, 5))
    plt.errorbar(n, medias, yerr=stds, fmt="o", markersize=3, capsize=2,
                 ecolor="gray", label=etiqueta)
    plt.xlabel("n")
    plt.ylabel("pi (probabilidad de aparición)")
    plt.title(f"Probabilidad de aparición de k-mers - {etiqueta} (k={k})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


def graficar_kmer_individual(carpeta, k, etiqueta):
    """
    Representa una línea por muestra con la probabilidad de cada k-mer (a
    diferencia de graficar_kmer(), que promedia todas las muestras).

    Parameters
    ----------
    carpeta : str
        Ruta a la carpeta con los ficheros de secuencias (.fastq o .fasta).
    k : int
        Longitud de los k-mers.
    etiqueta : str
        Texto identificativo de la condición ("Flight" o "Ground"), usado en
        el título.

    Returns
    -------
    None
        Muestra la figura por pantalla con plt.show().
    """

    _, lista_probs, todos_kmers, _, _, n, nombres_archivos = \
        analizar_kmers_carpeta(carpeta, k)

    plt.figure(figsize=(14, 5))
    for i, d in enumerate(lista_probs):
        valores = [d.get(kmer, 0) for kmer in todos_kmers]
        plt.plot(n, valores, label=nombres_archivos[i], alpha=0.6)

    plt.xlabel("n")
    plt.ylabel("pi (probabilidad de aparición)")
    plt.title(f"Probabilidad de aparición de k-mers - {etiqueta} (k={k})")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# ============================================================================
# 4. ENTROPÍA DE SHANNON
# ============================================================================


def shannon_entropy(probabilidades):
    """
    Calcula la entropía de Shannon global de un conjunto de probabilidades, en
    bits: -suma(p_i * log2(p_i)).

    Parameters
    ----------
    probabilidades : numpy.ndarray
        Vector de probabilidades p_i. Solo se tienen en cuenta los valores
        estrictamente positivos (log2(0) no está definido).

    Returns
    -------
    float
        Entropía de Shannon en bits.
    """

    probs = probabilidades[probabilidades > 0]  # Evita el logaritmo de 0.
    return -np.sum(probs * np.log2(probs))


def calcular_entropias_por_kmer(lista_probs, todos_kmers):
    """
    Calcula la contribución entrópica individual -p*log2(p) de cada k-mer,
    promediada entre varias muestras.

    Parameters
    ----------
    lista_probs : list of dict
        Lista de perfiles {k-mer: probabilidad}, uno por muestra. Los k-mers
        ausentes se tratan como probabilidad 0.
    todos_kmers : list of str
        Lista ordenada de todos los k-mers posibles. Fija el orden de las
        columnas del resultado.

    Returns
    -------
    tuple of numpy.ndarray
        (medias, stds): contribución entrópica media de cada k-mer y su
        desviación estándar entre muestras.
    """

    lista_entropias = []
    for probs in lista_probs:
        valores = np.array([probs.get(kmer, 0) for kmer in todos_kmers], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            entropias = np.where(valores > 0, -valores * np.log2(valores), 0.0)
        lista_entropias.append(entropias)

    lista_entropias = np.array(lista_entropias)
    medias = np.mean(lista_entropias, axis=0)
    stds = np.std(lista_entropias, axis=0)
    return medias, stds


def calcular_entropias_acumuladas(lista_probs, todos_kmers):
    """
    Calcula la curva de entropía de Shannon acumulada S(n) de cada muestra: la
    suma de las contribuciones entrópicas de los n primeros k-mers en el orden
    alfabético de todos_kmers.

    Parameters
    ----------
    lista_probs : list of dict
        Lista de perfiles {k-mer: probabilidad}, uno por muestra.
    todos_kmers : list of str
        Lista ordenada de todos los k-mers posibles.

    Returns
    -------
    tuple of numpy.ndarray
        (medias, stds): la curva de entropía acumulada media y su desviación
        estándar entre muestras (un valor por posición n). Se calcula con
        cumsum en vez de recalcular la entropía en cada punto, que sería
        mucho más lento.
    """

    curvas = []
    for probs in lista_probs:
        valores = np.array([probs.get(kmer, 0) for kmer in todos_kmers], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            contribuciones = np.where(valores > 0, -valores * np.log2(valores), 0.0)
        curvas.append(np.cumsum(contribuciones))

    curvas = np.array(curvas)
    medias = np.mean(curvas, axis=0)
    stds = np.std(curvas, axis=0)
    return medias, stds


def comparacion(carpeta_ground, carpeta_flight, k, semilla=None):
    """
    Compara el perfil de probabilidad de aparición de los k-mers entre
    Ground, Flight y el control Aleatorio, junto con la curva teórica, y
    exporta la tabla de k-mers ausentes.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    k : int
        Longitud de los k-mers.
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.

    Returns
    -------
    None
        Genera la figura Comparacion_k{k}_coli.png, la tabla CSV de
        k-mers ausentes y las entropías impresas por pantalla.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    (lista_probs_ground, lista_probs_flight,
     lista_probs_random, todos_kmers, n) = probabilidades_kmer_por_grupo(
        carpeta_ground, carpeta_flight, k, semilla)

    medias_ground, stds_ground = calcular_medias_stds_kmers(lista_probs_ground, todos_kmers)
    medias_flight, stds_flight = calcular_medias_stds_kmers(lista_probs_flight, todos_kmers)
    medias_random, stds_random = calcular_medias_stds_kmers(lista_probs_random, todos_kmers)

    entropia_ground = shannon_entropy(medias_ground)
    entropia_flight = shannon_entropy(medias_flight)
    entropia_aleatorio = shannon_entropy(medias_random)
    print(f"Entropía de Shannon Ground (k={k}): {entropia_ground:.4f} bits")
    print(f"Entropía de Shannon Flight (k={k}): {entropia_flight:.4f} bits")
    print(f"Entropía de Shannon aleatorio (k={k}): {entropia_aleatorio:.4f} bits")

    probs_teoricos = calcular_probabilidades_teoricas(todos_kmers)

    fig, ax = plt.subplots(figsize=(14, 5))

    # Zonas coloreadas por primer nucleótido de cada k-mer (orden alfabético).
    n_zona = 4 ** (k - 1)
    col_zona = {"A": "#d4f1c4", "C": "#c4d4f1", "G": "#f4f4c4", "T": "#f1c4c4"}
    for i, b in enumerate("ACGT"):
        ax.axvspan(i * n_zona + 0.5, (i + 1) * n_zona + 0.5,
                   alpha=0.35, color=col_zona[b], zorder=0)
        ax.axvline(i * n_zona + 0.5, color="gray", linestyle="--", alpha=0.4,
                   linewidth=0.7, zorder=1)
    parches = [Patch(facecolor=col_zona[b], alpha=0.5, label=f"{b}···") for b in "ACGT"]

    ax.plot(n, medias_ground, "-o", label="Ground", color="green", markersize=4, zorder=2)
    ax.plot(n, medias_flight, "-o", label="Flight", color="blue", markersize=4, zorder=2)
    ax.plot(n, medias_random, "-o", label="Aleatorio", color="orange", markersize=4, zorder=2)
    ax.plot(n, probs_teoricos, "--", label="Teórico", color="red", linewidth=2, zorder=2)
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("pi (probabilidad de aparición)")
    ax.set_title(f"Comparación de k-mers entre Ground, Flight y Aleatorio (k={k})")
    ax.grid(True, which="both", alpha=0.3)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h + parches, l + [p.get_label() for p in parches], fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(f"Comparacion_k{k}_coli.png", dpi=300, bbox_inches="tight")
    plt.close()

    datos_consulta = {
        "Ground": {"media": medias_ground, "std": stds_ground},
        "Flight": {"media": medias_flight, "std": stds_flight},
        "Aleatorio": {"media": medias_random, "std": stds_random},
        "Teórico": {"media": probs_teoricos},
    }
    exportar_tabla_kmers_filtrada(k, todos_kmers, datos_consulta, "tabla")


def comparacion_shannon(carpeta_ground, carpeta_flight, k, semilla=None):
    """
    Compara la contribución entrópica de Shannon (-p*log2(p)) de cada k-mer
    entre Ground, Flight y el control Aleatorio, junto con la curva teórica.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    k : int
        Longitud de los k-mers.
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.

    Returns
    -------
    None
        Genera la figura Comparacion_shannon_k{k}_shannon.png y la tabla
        CSV tabla_shannon_k{k}.csv.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    (lista_probs_ground, lista_probs_flight,
     lista_probs_random, todos_kmers, n) = probabilidades_kmer_por_grupo(
        carpeta_ground, carpeta_flight, k, semilla)

    medias_ground, stds_ground = calcular_entropias_por_kmer(lista_probs_ground, todos_kmers)
    medias_flight, stds_flight = calcular_entropias_por_kmer(lista_probs_flight, todos_kmers)
    medias_genoma, stds_genoma = calcular_entropias_por_kmer(lista_probs_random, todos_kmers)

    probs_teoricos = calcular_probabilidades_teoricas(todos_kmers)
    dict_teorico = {kmer: prob for kmer, prob in zip(todos_kmers, probs_teoricos)}
    entropias_teoricas, _ = calcular_entropias_por_kmer([dict_teorico], todos_kmers)

    fig, ax = plt.subplots(figsize=(14, 5))

    n_zona = 4 ** (k - 1)
    col_zona = {"A": "#d4f1c4", "C": "#c4d4f1", "G": "#f4f4c4", "T": "#f1c4c4"}
    for i, b in enumerate("ACGT"):
        ax.axvspan(i * n_zona + 0.5, (i + 1) * n_zona + 0.5,
                   alpha=0.35, color=col_zona[b], zorder=0)
        ax.axvline(i * n_zona + 0.5, color="gray", linestyle="--", alpha=0.4,
                   linewidth=0.7, zorder=1)
    parches = [Patch(facecolor=col_zona[b], alpha=0.5, label=f"{b}···") for b in "ACGT"]

    ax.plot(n, medias_ground, "-o", label="Ground", color="green", markersize=4, zorder=2)
    ax.plot(n, medias_flight, "-o", label="Flight", color="blue", markersize=4, zorder=2)
    ax.plot(n, medias_genoma, "-o", label="Aleatorio", color="orange", markersize=4, zorder=2)
    ax.plot(n, entropias_teoricas, "--", label="Teórico", color="red", linewidth=2, zorder=2)
    ax.set_yscale("log")
    ax.set_xlabel("n")
    ax.set_ylabel("Entropía de Shannon (bits)")
    ax.set_title(f"Comparación de entropía de Shannon entre Ground, Flight y Aleatorio (k={k})")
    ax.grid(True, which="both", alpha=0.3)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h + parches, l + [p.get_label() for p in parches], fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(f"Comparacion_shannon_k{k}_shannon.png", dpi=300, bbox_inches="tight")
    plt.close()

    datos_consulta = {
        "Ground": {"media": medias_ground, "std": stds_ground},
        "Flight": {"media": medias_flight, "std": stds_flight},
        "Aleatorio": {"media": medias_genoma, "std": stds_genoma},
        "Teórico": {"media": entropias_teoricas},
    }
    exportar_tabla_kmers(k, todos_kmers, datos_consulta, "tabla_shannon")


def comparacion_shannon_acumulada(carpeta_ground, carpeta_flight, k, semilla=None):
    """
    Compara la curva de entropía de Shannon acumulada entre Ground, Flight y
    el control Aleatorio, junto con la curva teórica de referencia.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    k : int
        Longitud de los k-mers.
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.

    Returns
    -------
    None
        Genera la figura shannon_acumulada_k{k}_6random.png.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    (lista_probs_ground, lista_probs_flight,
     lista_probs_random, todos_kmers, n) = probabilidades_kmer_por_grupo(
        carpeta_ground, carpeta_flight, k, semilla)

    media_ground, std_ground = calcular_entropias_acumuladas(lista_probs_ground, todos_kmers)
    media_flight, std_flight = calcular_entropias_acumuladas(lista_probs_flight, todos_kmers)
    media_random, std_random = calcular_entropias_acumuladas(lista_probs_random, todos_kmers)

    probs_teoricos = calcular_probabilidades_teoricas(todos_kmers)
    dict_teorico = {kmer: prob for kmer, prob in zip(todos_kmers, probs_teoricos)}
    entropia_acum_teorica, _ = calcular_entropias_acumuladas([dict_teorico], todos_kmers)

    fig, ax = plt.subplots(figsize=(14, 5))

    n_zona = 4 ** (k - 1)
    col_zona = {"A": "#d4f1c4", "C": "#c4d4f1", "G": "#f4f4c4", "T": "#f1c4c4"}
    for i, b in enumerate("ACGT"):
        ax.axvspan(i * n_zona + 0.5, (i + 1) * n_zona + 0.5,
                   alpha=0.35, color=col_zona[b], zorder=0)
        ax.axvline(i * n_zona + 0.5, color="gray", linestyle="--", alpha=0.4,
                   linewidth=0.7, zorder=1)
    parches = [Patch(facecolor=col_zona[b], alpha=0.5, label=f"{b}···") for b in "ACGT"]

    ax.plot(n, media_ground, "-o", color="green", label="Ground", markersize=4, zorder=2)
    ax.plot(n, media_flight, "-o", color="blue", label="Flight", markersize=4, zorder=2)
    ax.plot(n, media_random, "-o", color="orange", label="Aleatorio", markersize=4, zorder=2)
    ax.plot(n, entropia_acum_teorica, "--", label="Teórico", color="red", linewidth=2, zorder=2)
    ax.set_xlabel("n")
    ax.set_ylabel("Entropía acumulada de Shannon (bits)")
    ax.set_title(f"Comparación de entropía acumulada de Shannon por k-mer (k={k})")
    ax.grid(True, which="both", alpha=0.3)
    h, l = ax.get_legend_handles_labels()
    ax.legend(h + parches, l + [p.get_label() for p in parches], fontsize=8, ncol=2)
    plt.tight_layout()
    plt.savefig(f"shannon_acumulada_k{k}_6random.png", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Figura guardada: shannon_acumulada_k{k}_6random.png")


def entropia_normalizada_multik(carpeta_ground, carpeta_flight, lista_k, semilla=None, correccion=True):
    """
    Representa la entropía de Shannon normalizada, (H - H_teorica) /
    H_teorica, para varios valores de k, comparando Ground, Flight y el
    control Aleatorio. Si correccion=True, resta la media del control
    Aleatorio a cada grupo (el control queda en 0).

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    lista_k : list of int
        Valores de k a analizar, por ejemplo [3, 6, 9, 12].
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.
    correccion : bool, optional
        Si es True (por defecto), resta la media del control Aleatorio a cada
        grupo antes de graficar.

    Returns
    -------
    None
        Genera la figura entropia_normalizad_multik_coli(sin_correccion).png.
        Al final pregunta por terminal si se quiere exportar también la
        tabla de resultados a CSV.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    resultados_ground = []
    resultados_flight = []
    resultados_aleatorio = []
    stds_ground = []
    stds_flight = []
    stds_aleatorio = []

    for k in lista_k:
        print(f"\nProcesando k={k}...")

        (lista_probs_ground, lista_probs_flight,
         lista_probs_random, todos_kmers, _) = probabilidades_kmer_por_grupo(
            carpeta_ground, carpeta_flight, k, semilla)

        probs_teoricos = calcular_probabilidades_teoricas(todos_kmers)
        entropia_teorica = shannon_entropy(probs_teoricos)

        def entropias_normalizadas(lista_probs):
            valores = []
            for probs in lista_probs:
                medias_ind, _ = calcular_medias_stds_kmers([probs], todos_kmers)
                H = shannon_entropy(medias_ind)
                valores.append((H - entropia_teorica) / entropia_teorica)
            return valores

        entropias_norm_ground = entropias_normalizadas(lista_probs_ground)
        entropias_norm_flight = entropias_normalizadas(lista_probs_flight)
        entropias_norm_random = entropias_normalizadas(lista_probs_random)

        media_ground = np.mean(entropias_norm_ground)
        std_ground = np.std(entropias_norm_ground)
        media_flight = np.mean(entropias_norm_flight)
        std_flight = np.std(entropias_norm_flight)
        media_random = np.mean(entropias_norm_random)
        std_random = np.std(entropias_norm_random)

        if correccion:
            media_ground -= media_random
            media_flight -= media_random
            media_random -= media_random

        resultados_ground.append(media_ground)
        resultados_flight.append(media_flight)
        resultados_aleatorio.append(media_random)
        stds_ground.append(std_ground)
        stds_flight.append(std_flight)
        stds_aleatorio.append(std_random)

        print(f"  Ground: H_norm_media={media_ground:.6f}, std={std_ground:.6f}")
        print(f"  Flight: H_norm_media={media_flight:.6f}, std={std_flight:.6f}")
        print(f"  Aleatorio: H_norm_media={media_random:.6f}, std={std_random:.6f}")

    plt.figure(figsize=(12, 6))
    k_array = np.array(lista_k)
    stds_ground_array = np.array(stds_ground)
    stds_flight_array = np.array(stds_flight)
    stds_aleatorio_array = np.array(stds_aleatorio)

    plt.plot(k_array, resultados_ground, "-o", color="green", label="Ground", markersize=8)
    plt.fill_between(k_array, np.array(resultados_ground) - stds_ground_array,
                     np.array(resultados_ground) + stds_ground_array,
                     color="green", alpha=0.2)

    plt.plot(k_array, resultados_flight, "-s", color="blue", label="Flight", markersize=8)
    plt.fill_between(k_array, np.array(resultados_flight) - stds_flight_array,
                     np.array(resultados_flight) + stds_flight_array,
                     color="blue", alpha=0.2)

    plt.plot(k_array, resultados_aleatorio, "-^", color="orange", label="Aleatorio", markersize=8)
    plt.fill_between(k_array, np.array(resultados_aleatorio) - stds_aleatorio_array,
                     np.array(resultados_aleatorio) + stds_aleatorio_array,
                     color="orange", alpha=0.2)

    plt.axhline(y=0, color="red", linestyle="--", linewidth=1, alpha=0.5, label="Teórico")

    plt.xlabel("k", fontsize=12)
    ylabel = (
        "((H - Hteorica) / Hteorica) - Haleatorio"
        if correccion
        else "(H - Hteorica) / Hteorica"
    )
    plt.ylabel(ylabel, fontsize=12)
    plt.title("Entropía de Shannon normalizada para diferentes valores de k", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=10)
    plt.xticks(lista_k)
    plt.tight_layout()
    plt.savefig("entropia_normalizad_multik_coli(sin_correccion).png", dpi=300, bbox_inches="tight")
    plt.show()

    exportar = input("\n¿Quieres exportar una tabla con los resultados? (s/n): ").strip().lower()
    if exportar == "s":
        data = []
        for i, k in enumerate(lista_k):
            data.append({
                "k": k,
                "Ground_normalizada": resultados_ground[i],
                "Flight_normalizada": resultados_flight[i],
                "Aleatorio_normalizada": resultados_aleatorio[i],
            })
        df = pd.DataFrame(data)
        nombre_archivo = "entropia_normalizada_multik.csv"
        df.to_csv(nombre_archivo, index=False)
        print(f"\nTabla exportada correctamente como '{nombre_archivo}' en el directorio actual.")


# ============================================================================
# 5. K-MERS AUSENTES (MISSING K-MERS)
# ============================================================================


def calcular_missing_kmers_teorico(carpeta_ground, carpeta_flight, lista_k, semilla=None):
    """
    Estima el número esperado de k-mers ausentes en cada secuencia mediante la
    fórmula exacta del problema del coleccionista de cupones:

        E[M(n, k)] = 4**k * exp(-(n - k + 1) / 4**k)

    donde n es la longitud de la secuencia y 4**k el número total de k-mers
    posibles. El cálculo se hace por secuencia y luego se promedia, ya que la
    fórmula no es lineal en n. Los valores devueltos son teóricos; los
    empíricos se obtienen con graficar_kmers_cero_por_genoma().

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    lista_k : list of int
        Valores de k a analizar, por ejemplo [3, 6, 9, 12].
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.

    Returns
    -------
    dict
        Diccionario indexado por k. Para cada k contiene las claves 'ground',
        'flight' y 'aleatorio' (cada una con 'individual', 'media' y 'std') y
        la clave 'espacio_total' (el valor 4**k). Al final pregunta por
        terminal si se quiere exportar a missing_kmers_teorico.csv.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    print("Cargando secuencias Ground...")
    secuencias_ground, nombres_ground = cargar_carpeta(carpeta_ground)
    longitudes_ground = [len(seq) for seq in secuencias_ground]
    for nombre, n in zip(nombres_ground, longitudes_ground):
        print(f"  {nombre}: {n:,} bp")

    print("\nCargando secuencias Flight...")
    secuencias_flight, nombres_flight = cargar_carpeta(carpeta_flight)
    longitudes_flight = [len(seq) for seq in secuencias_flight]
    for nombre, n in zip(nombres_flight, longitudes_flight):
        print(f"  {nombre}: {n:,} bp")

    print("\nGenerando genomas aleatorios...")
    genomas_g_alea, genomas_f_alea, _, _, _ = generar_genomas_aleatorios(
        carpeta_ground, carpeta_flight, muestra=True, semilla=semilla)
    longitudes_aleatorio = [len(g) for g in genomas_g_alea + genomas_f_alea]
    for i, n in enumerate(longitudes_aleatorio, 1):
        print(f"  Aleatorio {i}: {n:,} bp")

    resultados = {}

    header = ("k      | 4^k            | Ground E[M] ind.               | G media      "
              "| Flight E[M] ind.               | F media      | Alea media")
    print("\n" + header)
    print("-" * 140)

    for k in lista_k:
        espacio = 4 ** k

        def esperados(longitudes):
            # fórmula del coleccionista de cupones aplicada a cada longitud
            return [espacio * np.exp(-(n - k + 1) / espacio) for n in longitudes]

        missing_ground = esperados(longitudes_ground)
        missing_flight = esperados(longitudes_flight)
        missing_aleatorio = esperados(longitudes_aleatorio)

        media_g, std_g = np.mean(missing_ground), np.std(missing_ground)
        media_f, std_f = np.mean(missing_flight), np.std(missing_flight)
        media_a, std_a = np.mean(missing_aleatorio), np.std(missing_aleatorio)

        resultados[k] = {
            "ground": {"individual": missing_ground, "media": media_g, "std": std_g},
            "flight": {"individual": missing_flight, "media": media_f, "std": std_f},
            "aleatorio": {"individual": missing_aleatorio, "media": media_a, "std": std_a},
            "espacio_total": espacio,
        }

        mg_str = [f"{m:.1f}" for m in missing_ground]
        mf_str = [f"{m:.1f}" for m in missing_flight]
        print(f"{k:<6} | {espacio:<14,} | {str(mg_str):<30} | {media_g:<12,.1f} | "
              f"{str(mf_str):<30} | {media_f:<12,.1f} | {media_a:<10,.1f}")

    exportar = input("\n¿Quieres exportar los resultados a CSV? (s/n): ").strip().lower()
    if exportar == "s":
        filas = []
        for k in lista_k:
            res = resultados[k]
            for i, (mg, mf, ma) in enumerate(zip(
                    res["ground"]["individual"],
                    res["flight"]["individual"],
                    res["aleatorio"]["individual"]), 1):
                filas.append({
                    "k": k,
                    "espacio_total_4k": res["espacio_total"],
                    "muestra": i,
                    "longitud_ground": longitudes_ground[i - 1],
                    "E_missing_ground": mg,
                    "longitud_flight": longitudes_flight[i - 1],
                    "E_missing_flight": mf,
                    "longitud_aleatorio": longitudes_aleatorio[i - 1],
                    "E_missing_aleatorio": ma,
                })
            filas.append({
                "k": k,
                "espacio_total_4k": res["espacio_total"],
                "muestra": "media",
                "longitud_ground": np.mean(longitudes_ground),
                "E_missing_ground": res["ground"]["media"],
                "longitud_flight": np.mean(longitudes_flight),
                "E_missing_flight": res["flight"]["media"],
                "longitud_aleatorio": np.mean(longitudes_aleatorio),
                "E_missing_aleatorio": res["aleatorio"]["media"],
            })

        df = pd.DataFrame(filas)
        nombre_archivo = "missing_kmers_teorico.csv"
        df.to_csv(nombre_archivo, index=False)
        print(f"\nResultados exportados como '{nombre_archivo}'")

    return resultados


def graficar_kmers_cero(lista_k, nombre_base="tabla", directorio="."):
    """
    Dibuja, en función de k, el número de k-mers ausentes (probabilidad 0) en
    Ground, Flight y Aleatorio, contando sobre la unión de las muestras de
    cada grupo a partir de los CSV ya generados por
    exportar_tabla_kmers_filtrada().

    Parameters
    ----------
    lista_k : list of int
        Valores de k para los que existe un CSV filtrado.
    nombre_base : str, optional
        Prefijo de los ficheros CSV. Por defecto "tabla".
    directorio : str, optional
        Carpeta donde están los CSV. Por defecto el directorio actual.

    Returns
    -------
    None
        Genera la figura kmers_cero(segundo_intento).png y, opcionalmente,
        exporta kmers_probabilidad_cero.csv. Si un CSV está vacío, corrupto
        o pesa demasiado, se salta ese k y sigue con el siguiente.
    """

    resultados = {
        "k": [],
        "ground_ceros": [],
        "flight_ceros": [],
        "aleatorio_ceros": [],
    }

    for k in lista_k:
        nombre_archivo = os.path.join(directorio, f"{nombre_base}_k{k}_filtrado.csv")

        if not os.path.isfile(nombre_archivo):
            print(f"Archivo no encontrado: '{nombre_archivo}'. Se omite k={k}.")
            continue

        print(f"Procesando k={k}... ")

        try:
            df = pd.read_csv(nombre_archivo)
        except pd.errors.ParserError:
            print(" Error con el motor C, intentando con el motor Python...")
            try:
                df = pd.read_csv(nombre_archivo, engine="python")
            except Exception as e:
                print(f"Error al leer el archivo: {e}")
                print("El archivo puede estar corrupto o incompleto.")
                print(f"Saltando k={k}...\n")
                continue
        except MemoryError:
            print("Error de memoria: el archivo es demasiado grande.")
            print("Intenta con un k más pequeño o con más memoria RAM.")
            print(f"Saltando k={k}...\n")
            continue
        except Exception as e:
            print(f"Error inesperado: {e}")
            print(f"Saltando k={k}...\n")
            continue

        if df.empty:
            print("Archivo vacío (sin datos).")
            ground_ceros = flight_ceros = aleatorio_ceros = 0
        else:
            columnas_necesarias = ["Ground_media", "Flight_media", "Aleatorio_media"]
            columnas_faltantes = [c for c in columnas_necesarias if c not in df.columns]
            if columnas_faltantes:
                print(f"Columnas faltantes: {columnas_faltantes}")
                print(f"Columnas disponibles: {list(df.columns)}")
                print(f"Saltando k={k}...\n")
                continue
            ground_ceros = (df["Ground_media"] == 0).sum()
            flight_ceros = (df["Flight_media"] == 0).sum()
            aleatorio_ceros = (df["Aleatorio_media"] == 0).sum()

        resultados["k"].append(k)
        resultados["ground_ceros"].append(ground_ceros)
        resultados["flight_ceros"].append(flight_ceros)
        resultados["aleatorio_ceros"].append(aleatorio_ceros)

        print(f"  K-mers con Ground=0:    {ground_ceros:,}")
        print(f"  K-mers con Flight=0:    {flight_ceros:,}")
        print(f"  K-mers con Aleatorio=0: {aleatorio_ceros:,}")

    fig, ax = plt.subplots(figsize=(12, 7))
    k_array = np.array(resultados["k"])
    ax.plot(k_array, np.array(resultados["ground_ceros"]), "o-", color="green",
            label="Ground (p=0)", linewidth=2, markersize=8)
    ax.plot(k_array, np.array(resultados["flight_ceros"]), "s-", color="blue",
            label="Flight (p=0)", linewidth=2, markersize=8)
    ax.plot(k_array, np.array(resultados["aleatorio_ceros"]), "^-", color="orange",
            label="Aleatorio (p=0)", linewidth=2, markersize=8)
    ax.set_xlabel("k", fontsize=12)
    ax.set_ylabel("Número de k-mers con p = 0", fontsize=12)
    ax.set_title("K-mers con p(k)=0", fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(fontsize=11, loc="best")
    ax.set_xticks(resultados["k"])
    ax.ticklabel_format(style="plain", axis="y")
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x):,}"))
    plt.tight_layout()
    plt.savefig("kmers_cero(segundo_intento).png", dpi=300, bbox_inches="tight")
    plt.show()

    exportar = input("\n¿Quieres exportar los datos de esta gráfica a CSV? (s/n): ").strip().lower()
    if exportar == "s":
        df_resultados = pd.DataFrame({
            "k": resultados["k"],
            "Ground_p0": resultados["ground_ceros"],
            "Flight_p0": resultados["flight_ceros"],
            "Aleatorio_p0": resultados["aleatorio_ceros"],
        })
        nombre_salida = "kmers_probabilidad_cero.csv"
        df_resultados.to_csv(nombre_salida, index=False)
        print(f"\nDatos exportados a: '{nombre_salida}'\n")


def graficar_kmers_cero_por_genoma(carpeta_ground, carpeta_flight, lista_k,
                                   guardar=None, mostrar=True, semilla=None):
    """
    Versión por genoma de graficar_kmers_cero(): cuenta los k-mers
    ausentes en cada genoma individual y promedia dentro de cada condición, en
    lugar de contar sobre la unión de todos los genomas. Es el criterio
    correcto, porque la unión subestima las ausencias del control Aleatorio
    cuando el genoma es más corto que el espacio de k-mers.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con las muestras de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con las muestras de vuelo (Flight).
    lista_k : list of int
        Valores de k a analizar.
    guardar : str or None, optional
        Ruta de la figura de salida. Si es None (por defecto) no se guarda.
    mostrar : bool, optional
        Si es True (por defecto) muestra la figura por pantalla con plt.show().
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.

    Returns
    -------
    dict
        Diccionario con las listas 'k', 'ground_ceros', 'flight_ceros' y
        'aleatorio_ceros' (número medio de k-mers ausentes por genoma en cada
        condición y para cada k).
    """

    def ausentes_medios(lista_genomas, k):
        espacio = 4 ** k
        ausentes = []
        for genoma in lista_genomas:
            probs = prob_kmer(genoma, k)
            if set(genoma) <= CONJUNTO_BASES:
                presentes = len(probs)  # Sin bases ambiguas (N): cuenta directa.
            else:
                presentes = sum(1 for km in probs if set(km) <= CONJUNTO_BASES)
            ausentes.append(espacio - presentes)
        return float(np.mean(ausentes))

    secuencias_ground, _ = cargar_carpeta(carpeta_ground)
    secuencias_flight, _ = cargar_carpeta(carpeta_flight)

    genomas_g, genomas_f, _, _, _ = generar_genomas_aleatorios(
        carpeta_ground, carpeta_flight, muestra=True, semilla=semilla)
    genomas_aleatorios = genomas_g + genomas_f

    resultados = {"k": [], "ground_ceros": [], "flight_ceros": [], "aleatorio_ceros": []}

    for k in lista_k:
        print(f"Procesando k={k}...")
        g = ausentes_medios(secuencias_ground, k)
        f = ausentes_medios(secuencias_flight, k)
        a = ausentes_medios(genomas_aleatorios, k)
        resultados["k"].append(k)
        resultados["ground_ceros"].append(g)
        resultados["flight_ceros"].append(f)
        resultados["aleatorio_ceros"].append(a)
        print(f"  Ground (media por genoma):    {g:,.0f}")
        print(f"  Flight (media por genoma):    {f:,.0f}")
        print(f"  Aleatorio (media por genoma): {a:,.0f}")

    fig, ax = plt.subplots(figsize=(12, 7))
    k_array = np.array(resultados["k"])
    ax.plot(k_array, resultados["ground_ceros"], "o-", color="green",
            label="Ground (p=0)", linewidth=2, markersize=8)
    ax.plot(k_array, resultados["flight_ceros"], "s-", color="blue",
            label="Flight (p=0)", linewidth=2, markersize=8)
    ax.plot(k_array, resultados["aleatorio_ceros"], "^-", color="orange",
            label="Aleatorio (p=0)", linewidth=2, markersize=8)
    ax.set_xlabel("k", fontsize=12)
    ax.set_ylabel("Número de k-mers con p = 0", fontsize=12)
    ax.set_title("K-mers con p(k)=0", fontsize=14, fontweight="bold")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.legend(fontsize=11, loc="best")
    ax.set_xticks(resultados["k"])
    ax.ticklabel_format(style="plain", axis="y")
    ax.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, p: f"{int(x):,}"))
    plt.tight_layout()
    if guardar:
        plt.savefig(guardar, dpi=300, bbox_inches="tight")
        print(f"Figura guardada: {guardar}")
    if mostrar:
        plt.show()

    return resultados


# ============================================================================
# 6. MÉTODO IEPWRMkmer
# ============================================================================


def D_kmer(seq, k):
    """
    Calcula la medida ponderada por posición D(k-mer) de una secuencia:

        D(k-mer) = ( suma de (posición_i / L) ) / (L - k + 1)

    donde las posiciones se cuentan desde 1, L es la longitud de la secuencia
    y (L - k + 1) es el número de ventanas de tamaño k. Un k-mer que aparece
    sobre todo al final de la secuencia tiene un valor mayor que uno que
    aparece al principio.

    Parameters
    ----------
    seq : str
        Secuencia de ADN.
    k : int
        Longitud de los k-mers.

    Returns
    -------
    dict
        Diccionario {k-mer: D(k-mer)} con un valor por cada k-mer presente en
        la secuencia.

    Raises
    ------
    ValueError
        Si la secuencia es más corta que k.
    """

    L = len(seq)
    if L < k:
        raise ValueError(f"La secuencia es más corta (L={L}) que k={k}")

    posiciones_kmers = defaultdict(list)
    for i in range(L - k + 1):
        kmer = seq[i:i + k]
        posiciones_kmers[kmer].append(i + 1)

    D_kmers = {}
    denominador = L - k + 1
    for kmer, posiciones in posiciones_kmers.items():
        sum_pos = sum(pos / L for pos in posiciones)
        D_kmers[kmer] = sum_pos / denominador

    return D_kmers


def H_kmer(lista_D_kmers, todos_kmers):
    """
    Calcula la entropía discriminativa H(k-mer): la rareza de un k-mer según
    en cuántos genomas aparece.

        H(k-mer) = -( F*log2(F) + (1 - F)*log2(1 - F) )

    con F = m / N (m genomas que contienen el k-mer, de un total de N). Vale 0
    cuando el k-mer está en todos los genomas o en ninguno, y es máximo cuando
    aparece en aproximadamente la mitad.

    Parameters
    ----------
    lista_D_kmers : list of dict
        Lista de diccionarios D(k-mer), uno por genoma (la salida de
        D_kmer()). Solo se usa para saber qué k-mers están presentes en
        cada genoma (D > 0).
    todos_kmers : list of str
        Lista de todos los k-mers posibles.

    Returns
    -------
    dict
        Diccionario {k-mer: H(k-mer)}.

    Raises
    ------
    ValueError
        Si la lista de genomas está vacía.
    """

    N = len(lista_D_kmers)
    if N == 0:
        raise ValueError("lista_D_kmers está vacía")

    H_kmers = {}
    for kmer in todos_kmers:
        m = sum(1 for D_dict in lista_D_kmers if D_dict.get(kmer, 0) > 0)
        F = m / N
        if F == 0 or F == 1:
            H_kmers[kmer] = 0.0
        else:
            H_kmers[kmer] = -(F * np.log2(F) + (1 - F) * np.log2(1 - F))

    return H_kmers


def E_kmer(D_kmer, H_kmer, todos_kmers):
    """
    Combina D(k-mer) y H(k-mer) en el vector de características E(k-mer) de
    un genoma: E(k-mer) = D(k-mer) * H(k-mer).

    Parameters
    ----------
    D_kmer : dict
        Valores D(k-mer) de un genoma concreto (0 para los k-mers ausentes).
    H_kmer : dict
        Valores H(k-mer) comunes a todo el conjunto de genomas.
    todos_kmers : list of str
        Lista de todos los k-mers posibles. Fija el orden del vector.

    Returns
    -------
    numpy.ndarray
        Vector de características E del genoma (un valor por cada k-mer
        posible).
    """

    E_vector = []
    for kmer in todos_kmers:
        D_val = D_kmer.get(kmer, 0.0)
        H_val = H_kmer[kmer]
        E_vector.append(D_val * H_val)
    return np.array(E_vector)


def manhattan(E_A, E_B):
    """
    Calcula la distancia de Manhattan (suma de diferencias absolutas) entre los
    vectores de características de dos genomas.

    Parameters
    ----------
    E_A : numpy.ndarray
        Vector de características E del genoma A.
    E_B : numpy.ndarray
        Vector de características E del genoma B.

    Returns
    -------
    float
        Distancia de Manhattan entre A y B, es decir, suma de |E_A - E_B|.

    Raises
    ------
    ValueError
        Si E_A y E_B tienen longitudes distintas.
    """

    return cityblock(E_A, E_B)


def analisis_IEPWRMkmer(carpeta_ground, carpeta_flight, k):
    """
    Ejecuta el método IEPWRMkmer completo para comparar Ground frente a Flight
    a un valor de k: carga las secuencias, calcula D(k-mer) y H(k-mer),
    combina ambos en el vector E(k-mer) de cada genoma y mide las distancias
    de Manhattan intra-grupo e inter-grupo.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con los genomas de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con los genomas de vuelo (Flight).
    k : int
        Longitud de los k-mers.

    Returns
    -------
    dict
        Diccionario con todos los resultados intermedios y finales del
        análisis (vectores E por grupo y promedio, valores H, listas de
        distancias intra e inter-grupo, estadísticas de H, etc.). Lo consume
        exportar_resultados_IEPWRMkmer().

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    """

    print(f"\n{'='*70}")
    print(f"ANÁLISIS IEPWRMkmer - k={k}")
    print(f"{'='*70}\n")

    print("Paso 1: Cargando secuencias...")
    seq_ground, nombres_ground = cargar_carpeta(carpeta_ground)
    seq_flight, nombres_flight = cargar_carpeta(carpeta_flight)
    N_ground = len(seq_ground)
    N_flight = len(seq_flight)
    N_total = N_ground + N_flight
    print(f"  - Ground: {N_ground} genomas")
    print(f"  - Flight: {N_flight} genomas")
    print(f"  - Total: {N_total} genomas\n")

    todos_kmers, _ = generar_todos_kmers(k)
    n_kmers_posibles = len(todos_kmers)
    print(f"  - K-mers posibles (4^{k}): {n_kmers_posibles:,}\n")

    print("Paso 2: Calculando D(k-mer) para cada genoma...")
    lista_D_ground = []
    for i, seq in enumerate(seq_ground, 1):
        lista_D_ground.append(D_kmer(seq, k))
        if i % 5 == 0 or i == N_ground:
            print(f"  - Ground: {i}/{N_ground} procesados")

    lista_D_flight = []
    for i, seq in enumerate(seq_flight, 1):
        lista_D_flight.append(D_kmer(seq, k))
        if i % 5 == 0 or i == N_flight:
            print(f"  - Flight: {i}/{N_flight} procesados")

    lista_D_todos = lista_D_ground + lista_D_flight
    print()

    print("Paso 3: Calculando H(k-mer) para el conjunto completo...")
    H_kmers = H_kmer(lista_D_todos, todos_kmers)
    H_valores = list(H_kmers.values())
    H_nonzero = [h for h in H_valores if h > 0]
    H_media = np.mean(H_valores)
    H_mediana = np.median(H_valores)
    H_max = np.max(H_valores)
    n_kmers_informativos = len(H_nonzero)
    print(f"  - H(k-mer) media: {H_media:.4f}")
    print(f"  - H(k-mer) mediana: {H_mediana:.4f}")
    print(f"  - H(k-mer) max: {H_max:.4f}")
    print(f"  - K-mers con H>0: {n_kmers_informativos:,} "
          f"({100 * n_kmers_informativos / n_kmers_posibles:.2f}%)\n")

    print("Paso 4: Calculando E(k-mer) para cada genoma...")
    E_ground = []
    for i, D_dict in enumerate(lista_D_ground, 1):
        E_ground.append(E_kmer(D_dict, H_kmers, todos_kmers))
        if i % 5 == 0 or i == N_ground:
            print(f"  - Ground: {i}/{N_ground} vectores E calculados")

    E_flight = []
    for i, D_dict in enumerate(lista_D_flight, 1):
        E_flight.append(E_kmer(D_dict, H_kmers, todos_kmers))
        if i % 5 == 0 or i == N_flight:
            print(f"  - Flight: {i}/{N_flight} vectores E calculados")
    print()

    print("Paso 5: Calculando distancias...")
    E_ground_mean = np.mean(E_ground, axis=0)
    E_flight_mean = np.mean(E_flight, axis=0)
    dist_mean = manhattan(E_ground_mean, E_flight_mean)

    dist_intra_ground = [manhattan(E_ground[i], E_ground[j])
                         for i in range(N_ground) for j in range(i + 1, N_ground)]
    dist_intra_flight = [manhattan(E_flight[i], E_flight[j])
                         for i in range(N_flight) for j in range(i + 1, N_flight)]
    dist_inter = [manhattan(E_ground[i], E_flight[j])
                  for i in range(N_ground) for j in range(N_flight)]

    print(f"\n{'='*70}")
    print("RESULTADOS")
    print(f"{'='*70}\n")
    print(f"Distancia entre vectores E promedio: {dist_mean:.6f}\n")
    print("Distancias intra-grupo:")
    print(f"  - Ground: {np.mean(dist_intra_ground):.6f} - {np.std(dist_intra_ground):.6f}")
    print(f"  - Flight: {np.mean(dist_intra_flight):.6f} - {np.std(dist_intra_flight):.6f}\n")
    print(f"Distancia inter-grupo: {np.mean(dist_inter):.6f} - {np.std(dist_inter):.6f}\n")

    return {
        "E_ground": E_ground,
        "E_flight": E_flight,
        "H_kmers": H_kmers,
        "todos_kmers": todos_kmers,
        "nombres_ground": nombres_ground,
        "nombres_flight": nombres_flight,
        "E_ground_mean": E_ground_mean,
        "E_flight_mean": E_flight_mean,
        "dist_intra_ground": dist_intra_ground,
        "dist_intra_flight": dist_intra_flight,
        "dist_inter": dist_inter,
        "dist_mean": dist_mean,
        "n_kmers_posibles": n_kmers_posibles,
        "H_media": H_media,
        "H_mediana": H_mediana,
        "H_max": H_max,
        "n_kmers_informativos": n_kmers_informativos,
    }


def exportar_resultados_IEPWRMkmer(resultados, k):
    """
    Exporta a CSV los resultados de un análisis IEPWRMkmer: la matriz de
    distancias entre todos los genomas, los valores de H(k-mer), los vectores
    E promedio de cada condición con su diferencia, los 20 k-mers más
    discriminativos y un resumen estadístico de las distancias.

    Parameters
    ----------
    resultados : dict
        Diccionario devuelto por analisis_IEPWRMkmer().
    k : int
        Valor de k del análisis (se usa en los nombres de los ficheros).

    Returns
    -------
    None
        Escribe cinco ficheros CSV en el directorio actual.

    Raises
    ------
    OSError
        Si alguno de los ficheros no se puede escribir.
    """

    print(f"\n{'='*70}")
    print("EXPORTANDO RESULTADOS")
    print(f"{'='*70}\n")

    # 1. Matriz de distancias completa.
    E_todos = resultados["E_ground"] + resultados["E_flight"]
    nombres_todos = resultados["nombres_ground"] + resultados["nombres_flight"]
    N = len(E_todos)
    matriz = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            dist = manhattan(E_todos[i], E_todos[j])
            matriz[i, j] = dist
            matriz[j, i] = dist
    df_matriz = pd.DataFrame(matriz, index=nombres_todos, columns=nombres_todos)
    archivo_matriz = f"matriz_distancias_IEPWRMkmer_k{k}.csv"
    df_matriz.to_csv(archivo_matriz)
    print(f" Matriz de distancias: {archivo_matriz}")

    # 2. Valores H(k-mer).
    df_H = pd.DataFrame({
        "k-mer": resultados["todos_kmers"],
        "H(k-mer)": [resultados["H_kmers"][kmer] for kmer in resultados["todos_kmers"]],
    })
    archivo_H = f"H_kmers_k{k}.csv"
    df_H.to_csv(archivo_H, index=False)
    print(f" Valores H(k-mer): {archivo_H}")

    # 3. Vectores E promedio y su diferencia.
    df_E = pd.DataFrame({
        "k-mer": resultados["todos_kmers"],
        "E_ground_mean": resultados["E_ground_mean"],
        "E_flight_mean": resultados["E_flight_mean"],
        "Diferencia": np.abs(resultados["E_ground_mean"] - resultados["E_flight_mean"]),
    })
    archivo_E = f"E_kmers_promedio_k{k}.csv"
    df_E.to_csv(archivo_E, index=False)
    print(f" Vectores E promedio: {archivo_E}")

    # 4. Los 20 k-mers más discriminativos.
    diferencias = np.abs(resultados["E_ground_mean"] - resultados["E_flight_mean"])
    indices_top = np.argsort(diferencias)[-20:][::-1]
    data_top = []
    for idx in indices_top:
        kmer = resultados["todos_kmers"][idx]
        data_top.append({
            "k-mer": kmer,
            "E_ground": resultados["E_ground_mean"][idx],
            "E_flight": resultados["E_flight_mean"][idx],
            "Diferencia": diferencias[idx],
            "H(k-mer)": resultados["H_kmers"][kmer],
        })
    df_top = pd.DataFrame(data_top)
    archivo_top = f"top_20_kmers_discriminativos_k{k}.csv"
    df_top.to_csv(archivo_top, index=False)
    print(f" Top 20 k-mers discriminativos: {archivo_top}")

    # 5. Resumen estadístico.
    df_resumen = pd.DataFrame({
        "Métrica": [
            "Distancia entre promedios",
            "Distancia intra-ground (media)",
            "Distancia intra-ground (std)",
            "Distancia intra-flight (media)",
            "Distancia intra-flight (std)",
            "Distancia inter-grupo (media)",
            "Distancia inter-grupo (std)",
        ],
        "Valor": [
            resultados["dist_mean"],
            np.mean(resultados["dist_intra_ground"]),
            np.std(resultados["dist_intra_ground"]),
            np.mean(resultados["dist_intra_flight"]),
            np.std(resultados["dist_intra_flight"]),
            np.mean(resultados["dist_inter"]),
            np.std(resultados["dist_inter"]),
        ],
    })
    archivo_resumen = f"resumen_estadistico_IEPWRMkmer_k{k}.csv"
    df_resumen.to_csv(archivo_resumen, index=False)
    print(f" Resumen estadístico: {archivo_resumen}")
    print(f"\n{'='*70}\n")


def frecuencia_nucleotidica_kmers_significativos(resultados, k, n_top=20):
    """
    Calcula la frecuencia nucleotídica de los k-mers más discriminativos de
    IEPWRMkmer (mayor |E_ground_mean - E_flight_mean|) y la compara con la
    frecuencia global de todos los k-mers analizados.

    Parameters
    ----------
    resultados : dict
        Diccionario devuelto por analisis_IEPWRMkmer().
    k : int
        Valor de k del análisis (para nombrar los ficheros exportados).
    n_top : int, optional
        Número de k-mers discriminativos a seleccionar (por defecto 20).

    Returns
    -------
    dict
        Con la lista de k-mers seleccionados ('kmers_significativos'), la
        frecuencia de cada base en todo el conjunto ('freq_global') y en los
        seleccionados ('freq_significativos'), y la tabla comparativa
        ('tabla', un DataFrame).

    Raises
    ------
    OSError
        Si el CSV de salida no se puede escribir.
    """

    todos_kmers = resultados["todos_kmers"]
    E_ground_mean = resultados["E_ground_mean"]
    E_flight_mean = resultados["E_flight_mean"]

    diferencias = np.abs(E_ground_mean - E_flight_mean)
    indices_top = np.argsort(diferencias)[-n_top:][::-1]
    kmers_sig = [todos_kmers[i] for i in indices_top]

    conteo_global = {b: 0 for b in BASES}
    for kmer in todos_kmers:
        for base in kmer:
            if base in conteo_global:
                conteo_global[base] += 1
    total_global = sum(conteo_global.values())
    freq_global = {b: conteo_global[b] / total_global for b in BASES}

    conteo_sig = {b: 0 for b in BASES}
    for kmer in kmers_sig:
        for base in kmer:
            if base in conteo_sig:
                conteo_sig[base] += 1
    total_sig = sum(conteo_sig.values())
    freq_sig = {b: conteo_sig[b] / total_sig if total_sig > 0 else 0.0 for b in BASES}

    filas = []
    for base in BASES:
        filas.append({
            "Nucleótido": base,
            "Freq_global": round(freq_global[base], 6),
            "Freq_significativos": round(freq_sig[base], 6),
            "Diferencia": round(freq_sig[base] - freq_global[base], 6),
        })
    df = pd.DataFrame(filas)

    print(f"\n{'='*70}")
    print(f"FRECUENCIA NUCLEOTÍDICA — Top {n_top} k-mers discriminativos (k={k})")
    print(f"{'='*70}")
    print(f"\nK-mers seleccionados ({n_top}): {kmers_sig}")
    print(f"\nComparación de frecuencias nucleotídicas:")
    print(df.to_string(index=False))

    nombre_csv = f"freq_nucleotidica_kmers_significativos_k{k}.csv"
    df.to_csv(nombre_csv, index=False)
    print(f"\n Tabla exportada: {nombre_csv}\n")

    return {
        "kmers_significativos": kmers_sig,
        "freq_global": freq_global,
        "freq_significativos": freq_sig,
        "tabla": df,
    }


def analisis_multik_IEPWRMkmer(carpeta_ground, carpeta_flight, lista_k):
    """
    Ejecuta el análisis IEPWRMkmer para varios valores de k, dibuja la
    evolución de la distancia inter-grupo frente a k y exporta un resumen en
    CSV.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con los genomas de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con los genomas de vuelo (Flight).
    lista_k : list of int
        Valores de k a analizar.

    Returns
    -------
    dict
        Diccionario indexado por k con el resultado completo de
        analisis_IEPWRMkmer() para cada valor. Genera también las figuras
        IEPWRMKmer_coli.png e IEPWRMkmer_3D_Hkmer.png y el csv
        resumen_multik_IEPWRMkmer.csv.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    OSError
        Si alguno de los ficheros de salida no se puede escribir.
    """

    print(f"\n{'='*70}")
    print("ANÁLISIS MULTI-K CON IEPWRMkmer")
    print(f"Valores de k: {lista_k}")
    print(f"{'='*70}\n")

    resultados_multi = {}
    kmers_posibles = []
    h_medias = []
    h_medianas = []
    h_maximos = []
    kmers_informativos = []
    porcentajes_informativos = []
    dist_intra_ground_media = []
    dist_intra_ground_std = []
    dist_intra_flight_media = []
    dist_intra_flight_std = []
    dist_inter_media = []
    dist_inter_std = []
    dist_entre_promedios = []

    for k in lista_k:
        print(f"\n{'-'*70}")
        print(f"Analizando k = {k}")
        print(f"{'-'*70}\n")

        resultados = analisis_IEPWRMkmer(carpeta_ground, carpeta_flight, k)
        resultados_multi[k] = resultados

        kmers_posibles.append(resultados["n_kmers_posibles"])
        h_medias.append(resultados["H_media"])
        h_medianas.append(resultados["H_mediana"])
        h_maximos.append(resultados["H_max"])
        kmers_informativos.append(resultados["n_kmers_informativos"])
        porcentajes_informativos.append(
            100 * resultados["n_kmers_informativos"] / resultados["n_kmers_posibles"])

        intra_g = resultados["dist_intra_ground"]
        intra_f = resultados["dist_intra_flight"]
        dist_intra_ground_media.append(np.mean(intra_g) if intra_g else 0)
        dist_intra_ground_std.append(np.std(intra_g) if intra_g else 0)
        dist_intra_flight_media.append(np.mean(intra_f) if intra_f else 0)
        dist_intra_flight_std.append(np.std(intra_f) if intra_f else 0)
        dist_inter_media.append(np.mean(resultados["dist_inter"]))
        dist_inter_std.append(np.std(resultados["dist_inter"]))
        dist_entre_promedios.append(resultados["dist_mean"])

    plt.figure(figsize=(10, 6))
    plt.plot(lista_k, dist_inter_media, "o-", linewidth=2, markersize=10, color="steelblue")
    plt.xlabel("k", fontsize=12)
    plt.ylabel("Distancia inter-grupo promedio", fontsize=12)
    plt.title("Evolución de la distancia Ground-Flight (IEPWRMkmer)", fontsize=14, fontweight="bold")
    plt.grid(True, alpha=0.3)
    plt.xticks(lista_k)
    plt.tight_layout()
    plt.savefig("IEPWRMKmer_coli.png", dpi=300, bbox_inches="tight")
    plt.show()

    # Plot 3D: distribución de H(k-mer) por valor de k.
    # Eje X = H(k-mer), Eje Y = Frecuencia (conteo), Eje Z = k
    fig_3d = plt.figure(figsize=(13, 8))
    ax_3d = fig_3d.add_subplot(111, projection="3d")
    colores_3d = plt.cm.viridis(np.linspace(0, 0.85, len(lista_k)))
    for idx_k, (k_val, color) in enumerate(zip(lista_k, colores_3d)):
        H_vals = list(resultados_multi[k_val]["H_kmers"].values())
        counts, edges = np.histogram(H_vals, bins=40)
        centers = (edges[:-1] + edges[1:]) / 2
        width = edges[1] - edges[0]
        log_counts = np.log10(counts + 1)  # +1 para evitar log(0)
        ax_3d.bar(centers, log_counts, zs=k_val, zdir="y",
                  alpha=0.7, width=width, color=color, edgecolor="none")
    ax_3d.set_xlabel("H(k-mer)", fontsize=11, labelpad=8)
    ax_3d.set_ylabel("k", fontsize=11, labelpad=8)
    ax_3d.set_zlabel("log₁₀(Frecuencia + 1)", fontsize=11, labelpad=8)
    ax_3d.set_yticks(lista_k)
    ax_3d.set_title("Distribución de H(k-mer) por valor de k (IEPWRMkmer)",
                    fontsize=13, fontweight="bold")
    plt.tight_layout()
    plt.savefig("IEPWRMkmer_3D_Hkmer.png", dpi=300, bbox_inches="tight")
    plt.show()

    df_resumen = pd.DataFrame({
        "k": lista_k,
        "kmers_posibles": kmers_posibles,
        "kmers_informativos": kmers_informativos,
        "porcentaje_informativos": porcentajes_informativos,
        "H_media": h_medias,
        "H_mediana": h_medianas,
        "H_max": h_maximos,
        "dist_intra_ground_media": dist_intra_ground_media,
        "dist_intra_ground_std": dist_intra_ground_std,
        "dist_intra_flight_media": dist_intra_flight_media,
        "dist_intra_flight_std": dist_intra_flight_std,
        "dist_inter_grupo_media": dist_inter_media,
        "dist_inter_grupo_std": dist_inter_std,
        "dist_entre_promedios": dist_entre_promedios,
    })
    df_resumen.to_csv("resumen_multik_IEPWRMkmer.csv", index=False)
    print("\n Resumen exportado: resumen_multik_IEPWRMkmer.csv\n")

    return resultados_multi


# ============================================================================
# 7. MÉTRICA BIOBIT (BB)
# ============================================================================


def encontrar_k_optimo(seq, k_min=3, k_max=12):
    """
    Busca el valor de k entre k_min y k_max que produce la mayor entropía de
    Shannon en una secuencia. Este k óptimo es la base de la métrica Biobit.

    Parameters
    ----------
    seq : str
        Secuencia de ADN.
    k_min : int, optional
        Valor mínimo de k a probar (por defecto 3).
    k_max : int, optional
        Valor máximo de k a probar (por defecto 12).

    Returns
    -------
    tuple
        (k_optimo, entropia_maxima), el valor de k que maximiza la entropía y
        la entropía alcanzada (en bits).
    """

    max_entropia = 0
    k_optimo = k_min

    print(f"  Buscando k óptimo (k={k_min} a k={k_max})...")
    for k in range(k_min, k_max + 1):
        probs_dict = prob_kmer(seq, k)
        probs_array = np.array(list(probs_dict.values()))
        entropia = shannon_entropy(probs_array)
        if entropia > max_entropia:
            max_entropia = entropia
            k_optimo = k
        if k % 3 == 0:
            print(f"    k={k}: H={entropia:.4f}")

    return k_optimo, max_entropia


def calcular_biobit(seq, k_range=(3, 12), verbose=False):
    """
    Calcula la métrica Biobit (BB) de una secuencia: el equilibrio entre su
    componente entrópica E(G) = E2L(G) - L(G) y antientrópica
    A(G) = 2*L(G) - E2L(G), donde L(G) = log4(G) y E2L(G) es la entropía
    máxima (la del k óptimo).

        BB(G) = sqrt(L(G)) * sqrt(A(G)/L(G)) * (1 - 2*A(G)/L(G))**3

    Parameters
    ----------
    seq : str
        Secuencia de ADN.
    k_range : tuple of int, optional
        Rango (k_min, k_max) en el que buscar el k óptimo (por defecto
        (3, 12)).
    verbose : bool, optional
        Si es True, imprime el detalle de todos los componentes. Por defecto
        False.

    Returns
    -------
    dict
        Diccionario con la longitud, el k óptimo, L(G), E2L(G), las
        componentes E(G) y A(G), el cociente A(G)/L(G) y el valor final 'BB'.
    """

    G = len(seq)
    L_G = np.log(G) / np.log(4)

    k_optimo, E2L_G = encontrar_k_optimo(seq, k_min=k_range[0], k_max=k_range[1])

    E_G = E2L_G - L_G       # Componente entrópica (desorden).
    A_G = 2 * L_G - E2L_G   # Componente antientrópica (orden).
    ratio_A_L = A_G / L_G if L_G > 0 else 0

    if ratio_A_L >= 0 and L_G > 0:
        BB = np.sqrt(L_G) * np.sqrt(ratio_A_L) * (1 - 2 * ratio_A_L) ** 3
    else:
        BB = 0.0

    resultados = {
        "longitud": G,
        "k_optimo": k_optimo,
        "L_G": L_G,
        "E2L_G": E2L_G,
        "E_G": E_G,
        "A_G": A_G,
        "ratio_A_L": ratio_A_L,
        "BB": BB,
    }

    if verbose:
        print("\n")
        print(f"  Longitud (G): {G:,} bp")
        print(f"  K óptimo: {k_optimo}")
        print(f"  L(G) = log4(G): {L_G:.4f}")
        print(f"  E2L(G) (entropía real): {E2L_G:.4f}")
        print(f"  E(G) (entrópico): {E_G:.4f}")
        print(f"  A(G) (antientrópico): {A_G:.4f}")
        print(f"  A(G)/L(G): {ratio_A_L:.4f}")
        print(f"  Biobit (BB): {BB:.6f}")

    return resultados


def biobit_carpetas(carpeta_ground, carpeta_flight, k_range=(3, 12)):
    """
    Calcula y compara el Biobit de las muestras Ground y Flight de dos
    carpetas.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con los FASTA de tierra (Ground).
    carpeta_flight : str
        Ruta a la carpeta con los FASTA de vuelo (Flight).
    k_range : tuple of int, optional
        Rango (k_min, k_max) para buscar el k óptimo de cada secuencia (por
        defecto (3, 12)).

    Returns
    -------
    dict
        Diccionario con los Biobit individuales de cada grupo ('bb_ground',
        'bb_flight'), sus medias y desviaciones, y la diferencia Flight -
        Ground. Solo lee ficheros .fasta y genera biobit_resultados.csv.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    OSError
        Si el CSV de resultados no se puede escribir.
    """

    def cargar_fastas(carpeta):
        # carga solo los .fasta de la carpeta (ordenados por nombre)
        seqs, nombres = [], []
        for archivo in sorted(os.listdir(carpeta)):
            if archivo.endswith(".fasta"):
                ruta = os.path.join(carpeta, archivo)
                print(f"  {archivo}: ", end="")
                seq = FASTA(ruta)
                seqs.append(seq)
                nombres.append(archivo)
                print(f"{len(seq):,} bp")
        return seqs, nombres

    print("Cargando secuencias Ground...")
    seq_ground, nombres_ground = cargar_fastas(carpeta_ground)
    print("\nCargando secuencias Flight...")
    seq_flight, nombres_flight = cargar_fastas(carpeta_flight)
    print(f"\nTotal: {len(seq_ground)} Ground, {len(seq_flight)} Flight\n")

    print("=" * 70)
    print("CALCULANDO BIOBIT PARA GROUND")
    print("=" * 70)
    bb_ground = []
    for i, (seq, nombre) in enumerate(zip(seq_ground, nombres_ground), 1):
        print(f"\n[{i}/{len(seq_ground)}] {nombre}")
        res = calcular_biobit(seq, k_range=k_range, verbose=True)
        bb_ground.append(res["BB"])

    print(f"\n{'='*70}")
    print("CALCULANDO BIOBIT PARA FLIGHT")
    print("=" * 70)
    bb_flight = []
    for i, (seq, nombre) in enumerate(zip(seq_flight, nombres_flight), 1):
        print(f"\n[{i}/{len(seq_flight)}] {nombre}")
        res = calcular_biobit(seq, k_range=k_range, verbose=True)
        bb_flight.append(res["BB"])

    print(f"\n{'='*70}")
    print("ESTADÍSTICAS COMPARATIVAS")
    print("=" * 70)
    media_ground = np.mean(bb_ground)
    std_ground = np.std(bb_ground, ddof=1) if len(bb_ground) > 1 else 0
    media_flight = np.mean(bb_flight)
    std_flight = np.std(bb_flight, ddof=1) if len(bb_flight) > 1 else 0
    print(f"\nBiobit Ground: {media_ground:.6f} - {std_ground:.6f}")
    print(f"  Valores: {[f'{x:.6f}' for x in bb_ground]}")
    print(f"\nBiobit Flight: {media_flight:.6f} - {std_flight:.6f}")
    print(f"  Valores: {[f'{x:.6f}' for x in bb_flight]}")
    diferencia = media_flight - media_ground
    print(f"\nDiferencia (Flight - Ground): {diferencia:+.6f}")

    df = pd.DataFrame({
        "Muestra": nombres_ground + nombres_flight,
        "Grupo": ["Ground"] * len(bb_ground) + ["Flight"] * len(bb_flight),
        "BB": bb_ground + bb_flight,
    })
    df.to_csv("biobit_resultados.csv", index=False)
    print("\n Resultados exportados: biobit_resultados.csv\n")

    return {
        "bb_ground": bb_ground,
        "bb_flight": bb_flight,
        "media_ground": media_ground,
        "std_ground": std_ground,
        "media_flight": media_flight,
        "std_flight": std_flight,
        "diferencia": diferencia,
    }


# ============================================================================
# 8. FIRMA GENÓMICA (GENOME SIGNATURE, GS)
# ============================================================================
#
# La firma genómica resume, en un único número, cuánto se aleja la composición
# de k-mers de un genoma de una distribución completamente uniforme (aquella en
# la que todos los k-mers posibles aparecerían el mismo número de veces). Un
# valor alto indica que el genoma tiene una "huella" composicional marcada
# (algunos k-mers muy sobrerrepresentados y otros casi ausentes); un valor bajo
# indica una composición más homogénea. Aquí se emplea como métrica comparativa
# entre las condiciones Ground y Flight.
#
# Referencia: Román-Escrivá et al. (2025) Biology 14(4):338.
# ============================================================================


def calcular_GS_k(seq, k):
    """
    Calcula la firma genómica de una secuencia para un único valor de k: la
    desviación absoluta media del reparto observado de k-mers respecto al
    reparto uniforme.

        GS_k = (1 / N_total) * suma_i |n_i - EV|

    donde N_total es el número total de k-mers encontrados, EV = N_total /
    4**k es el valor esperado bajo distribución uniforme y n_i es el número de
    apariciones del k-mer i-ésimo (los k-mers ausentes contribuyen con
    |0 - EV| = EV cada uno).

    Parameters
    ----------
    seq : str
        Secuencia de ADN.
    k : int
        Longitud de los k-mers.

    Returns
    -------
    float
        Valor de GS para ese k. Devuelve 0.0 si la secuencia no contiene
        ningún k-mer (por ejemplo, si es más corta que k).
    """

    conteos = kmers(seq, k)

    N_total = sum(conteos.values())
    n_posibles = 4 ** k

    if N_total == 0:
        return 0.0

    EV = N_total / n_posibles

    suma_desviaciones = 0.0
    for kmer_str, ni in conteos.items():
        suma_desviaciones += abs(ni - EV)

    n_ausentes = n_posibles - len(conteos)
    suma_desviaciones += n_ausentes * EV

    return suma_desviaciones / N_total


def calcular_GS(seq, k_min=2, k_max=16, verbose=False):
    """
    Calcula la firma genómica completa de una secuencia recorriendo varios k,
    corrigiendo el efecto de la composición de bases: genera un genoma
    aleatorio de la misma longitud y composición, calcula GS_k del genoma real
    (GS_g) y del aleatorio (GS_r) para cada k, y se queda con la diferencia
    GS_p = GS_g - GS_r. La firma final es el máximo de GS_p y el k donde se
    alcanza.

    Parameters
    ----------
    seq : str
        Secuencia de ADN.
    k_min : int, optional
        Valor mínimo de k a evaluar (por defecto 2).
    k_max : int, optional
        Valor máximo de k a evaluar (por defecto 16).
    verbose : bool, optional
        Si es True, imprime por pantalla una tabla con GS_g, GS_r y GS_p para
        cada k (por defecto False).

    Returns
    -------
    dict
        La firma final ('GS') y el k donde se alcanza ('k_optimo'), más los
        perfiles completos por k: perfil_k, perfil_GSg (genoma real),
        perfil_GSr (genoma aleatorio) y perfil_GSp (la diferencia).

    Para que el genoma aleatorio de control salga siempre igual hay que fijar
    una semilla con random.seed(...) antes de llamar (gs_carpetas ya lo hace).
    """

    mono, _, _, _, _ = conteo(seq)
    freq = frecuencia(mono, seq)
    seq_aleatoria = aleatorio(freq, len(seq))

    perfil_k = []
    perfil_GSg = []
    perfil_GSr = []
    perfil_GSp = []

    if verbose:
        print(f"    {'k':>3}  {'GS_g':>10}  {'GS_r':>10}  {'GS_p':>10}")
        print(f"    {'---':>3}  {'----------':>10}  {'----------':>10}  {'----------':>10}")

    for k in range(k_min, k_max + 1):
        GS_g = calcular_GS_k(seq, k)
        GS_r = calcular_GS_k(seq_aleatoria, k)
        GS_p = GS_g - GS_r

        perfil_k.append(k)
        perfil_GSg.append(GS_g)
        perfil_GSr.append(GS_r)
        perfil_GSp.append(GS_p)

        if verbose:
            print(f"    {k:>3}  {GS_g:>10.6f}  {GS_r:>10.6f}  {GS_p:>10.6f}")

    idx_max = np.argmax(perfil_GSp)
    GS_final = perfil_GSp[idx_max]
    k_optimo = perfil_k[idx_max]

    if verbose:
        print(f"\n    GS = {GS_final:.6f} (k óptimo = {k_optimo})")

    return {
        'GS': GS_final,
        'k_optimo': k_optimo,
        'perfil_k': perfil_k,
        'perfil_GSg': perfil_GSg,
        'perfil_GSr': perfil_GSr,
        'perfil_GSp': perfil_GSp,
    }


def gs_carpetas(carpeta_ground, carpeta_flight, k_min=2, k_max=12, semilla=None):
    """
    Compara la firma genómica entre las condiciones Ground y Flight: calcula
    la GS de cada genoma de las dos carpetas y de un control Aleatorio, resume
    los valores por condición, construye los gráficos comparativos y exporta
    los resultados a CSV.

    Parameters
    ----------
    carpeta_ground : str
        Ruta a la carpeta con los genomas de la condición Ground (.fasta o
        .fastq).
    carpeta_flight : str
        Ruta a la carpeta con los genomas de la condición Flight (.fasta o
        .fastq).
    k_min : int, optional
        Valor mínimo de k a evaluar (por defecto 2).
    k_max : int, optional
        Valor máximo de k a evaluar (por defecto 12).
    semilla : int or None, optional
        Semilla para la generación reproducible del control Aleatorio. Por
        defecto None.

    Returns
    -------
    dict
        Las listas de GS por genoma (gs_ground, gs_flight), sus medias y
        desviaciones (media_ground, std_ground, media_flight, std_flight),
        la diferencia Flight - Ground, y los resultados detallados de
        calcular_GS() de cada genoma (resultados_ground, resultados_flight).

    Genera las figuras gs_perfiles_individuales.png y gs_perfiles_medios.png,
    y los csv gs_resultados.csv y gs_perfiles.csv.

    Raises
    ------
    FileNotFoundError
        Si alguna de las carpetas no existe.
    OSError
        Si alguno de los ficheros de salida no se puede escribir.
    """

    if semilla is not None:
        random.seed(semilla)

    print(f"\n{'='*70}")
    print(f"ANÁLISIS FIRMA GENÓMICA (GS)")
    print(f"Rango de k: {k_min} a {k_max}")
    print(f"{'='*70}\n")

    print("Cargando secuencias...")
    seq_ground, nombres_ground = cargar_carpeta(carpeta_ground)
    seq_flight, nombres_flight = cargar_carpeta(carpeta_flight)

    print(f"  Ground: {len(seq_ground)} genomas")
    print(f"  Flight: {len(seq_flight)} genomas\n")

    mono_todos, _, _, _, _ = conteo("".join(seq_ground))
    freq_media = frecuencia(mono_todos, "".join(seq_ground))
    seq_aleatorio = [aleatorio(freq_media, len(s)) for s in seq_ground]
    nombres_aleatorio = [f"aleatorio_{i + 1}" for i in range(len(seq_aleatorio))]

    print("=" * 70)
    print("CALCULANDO GS PARA GROUND")
    print("=" * 70)

    resultados_ground = []
    gs_ground = []
    for i, (seq, nombre) in enumerate(zip(seq_ground, nombres_ground), 1):
        print(f"\n  [{i}/{len(seq_ground)}] {nombre} ({len(seq):,} bp)")
        res = calcular_GS(seq, k_min=k_min, k_max=k_max, verbose=True)
        resultados_ground.append(res)
        gs_ground.append(res['GS'])

    print(f"\n{'='*70}")
    print("CALCULANDO GS PARA FLIGHT")
    print("=" * 70)

    resultados_flight = []
    gs_flight = []
    for i, (seq, nombre) in enumerate(zip(seq_flight, nombres_flight), 1):
        print(f"\n  [{i}/{len(seq_flight)}] {nombre} ({len(seq):,} bp)")
        res = calcular_GS(seq, k_min=k_min, k_max=k_max, verbose=True)
        resultados_flight.append(res)
        gs_flight.append(res['GS'])

    print(f"\n{'='*70}")
    print("CALCULANDO GS PARA ALEATORIO")
    print("=" * 70)

    resultados_aleatorio = []
    gs_aleatorio = []
    for i, (seq, nombre) in enumerate(zip(seq_aleatorio, nombres_aleatorio), 1):
        print(f"\n  [{i}/{len(seq_aleatorio)}] {nombre} ({len(seq):,} bp)")
        res = calcular_GS(seq, k_min=k_min, k_max=k_max, verbose=True)
        resultados_aleatorio.append(res)
        gs_aleatorio.append(res['GS'])

    print(f"\n{'='*70}")
    print("ESTADÍSTICAS COMPARATIVAS")
    print("=" * 70)

    media_ground = np.mean(gs_ground)
    std_ground = np.std(gs_ground, ddof=1) if len(gs_ground) > 1 else 0
    media_flight = np.mean(gs_flight)
    std_flight = np.std(gs_flight, ddof=1) if len(gs_flight) > 1 else 0
    media_aleatorio = np.mean(gs_aleatorio)
    std_aleatorio = np.std(gs_aleatorio, ddof=1) if len(gs_aleatorio) > 1 else 0

    print(f"\nGS Ground:    {media_ground:.6f} ± {std_ground:.6f}")
    print(f"  Valores: {[f'{x:.6f}' for x in gs_ground]}")
    print(f"  k óptimos: {[r['k_optimo'] for r in resultados_ground]}")

    print(f"\nGS Flight:    {media_flight:.6f} ± {std_flight:.6f}")
    print(f"  Valores: {[f'{x:.6f}' for x in gs_flight]}")
    print(f"  k óptimos: {[r['k_optimo'] for r in resultados_flight]}")

    print(f"\nGS Aleatorio: {media_aleatorio:.6f} ± {std_aleatorio:.6f}")
    print(f"  Valores: {[f'{x:.6f}' for x in gs_aleatorio]}")

    diferencia = media_flight - media_ground
    print(f"\nDiferencia (Flight - Ground): {diferencia:+.6f}")

    # --- GRÁFICO 1: perfiles GS_p(k) individuales (sin boxplot) ---
    fig1, ax1 = plt.subplots(figsize=(11, 6))

    for i, (res, nombre) in enumerate(zip(resultados_ground, nombres_ground)):
        ax1.plot(res['perfil_k'], res['perfil_GSp'], 'o-',
                 color='green', alpha=0.4 + 0.2 * i, label=f'G: {nombre}')
    for i, (res, nombre) in enumerate(zip(resultados_flight, nombres_flight)):
        ax1.plot(res['perfil_k'], res['perfil_GSp'], 's--',
                 color='blue', alpha=0.4 + 0.2 * i, label=f'F: {nombre}')
    for i, (res, nombre) in enumerate(zip(resultados_aleatorio, nombres_aleatorio)):
        ax1.plot(res['perfil_k'], res['perfil_GSp'], '^:',
                 color='gold', alpha=0.4 + 0.2 * i, label=f'A: {nombre}')

    ax1.set_xlabel('k', fontsize=12)
    ax1.set_ylabel('GS_p (GS_g − GS_r)', fontsize=12)
    ax1.set_title('Perfiles de firma genómica (GS_p)', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=8, loc='best')
    ax1.grid(True, alpha=0.3, which='both')
    ax1.set_xticks(range(k_min, k_max + 1))
    plt.tight_layout()
    plt.savefig('gs_perfiles_individuales.png', dpi=300, bbox_inches='tight')
    plt.show()

    # --- GRÁFICO 2: perfiles medios por condición con bandas de error ---
    fig2, ax2 = plt.subplots(figsize=(11, 6))

    ks = resultados_ground[0]['perfil_k']

    def perfil_medio(resultados):
        mat = np.array([r['perfil_GSp'] for r in resultados])
        media = np.mean(mat, axis=0)
        std = np.std(mat, axis=0, ddof=1) if mat.shape[0] > 1 else np.zeros(mat.shape[1])
        return media, std

    mg, sg = perfil_medio(resultados_ground)
    mf, sf = perfil_medio(resultados_flight)
    ma, sa = perfil_medio(resultados_aleatorio)

    ax2.plot(ks, mg, 'o-', color='green', label='Ground (media)', linewidth=2)
    ax2.fill_between(ks, mg - sg, mg + sg, alpha=0.2, color='green')
    ax2.plot(ks, mf, 's--', color='blue', label='Flight (media)', linewidth=2)
    ax2.fill_between(ks, mf - sf, mf + sf, alpha=0.2, color='blue')
    ax2.plot(ks, ma, '^:', color='gold', label='Aleatorio (media)', linewidth=2)
    ax2.fill_between(ks, ma - sa, ma + sa, alpha=0.2, color='gold')

    ax2.set_xlabel('k', fontsize=12)
    ax2.set_ylabel('GS_p (GS_g − GS_r)', fontsize=12)
    ax2.set_title('Perfil medio de firma genómica: Ground vs Flight vs Aleatorio',
                  fontsize=14, fontweight='bold')
    ax2.legend(fontsize=11)
    ax2.grid(True, alpha=0.3, which='both')
    ax2.set_xticks(range(k_min, k_max + 1))
    plt.tight_layout()
    plt.savefig('gs_perfiles_medios.png', dpi=300, bbox_inches='tight')
    plt.show()

    filas = []
    for nombre, res in zip(nombres_ground, resultados_ground):
        filas.append({'Muestra': nombre, 'Grupo': 'Ground',
                      'GS': res['GS'], 'k_optimo': res['k_optimo']})
    for nombre, res in zip(nombres_flight, resultados_flight):
        filas.append({'Muestra': nombre, 'Grupo': 'Flight',
                      'GS': res['GS'], 'k_optimo': res['k_optimo']})
    for nombre, res in zip(nombres_aleatorio, resultados_aleatorio):
        filas.append({'Muestra': nombre, 'Grupo': 'Aleatorio',
                      'GS': res['GS'], 'k_optimo': res['k_optimo']})

    df = pd.DataFrame(filas)
    df.to_csv('gs_resultados.csv', index=False)
    print(f"\nResultados exportados: gs_resultados.csv")

    filas_perfil = []
    todos_nombres = nombres_ground + nombres_flight + nombres_aleatorio
    todos_res = resultados_ground + resultados_flight + resultados_aleatorio
    todos_grupos = (['Ground'] * len(nombres_ground) + ['Flight'] * len(nombres_flight)
                    + ['Aleatorio'] * len(nombres_aleatorio))
    for nombre, grupo, res in zip(todos_nombres, todos_grupos, todos_res):
        for j, k_val in enumerate(res['perfil_k']):
            filas_perfil.append({
                'Muestra': nombre, 'Grupo': grupo, 'k': k_val,
                'GS_g': res['perfil_GSg'][j],
                'GS_r': res['perfil_GSr'][j],
                'GS_p': res['perfil_GSp'][j]
            })

    df_perfil = pd.DataFrame(filas_perfil)
    df_perfil.to_csv('gs_perfiles.csv', index=False)
    print(f"Perfiles exportados: gs_perfiles.csv\n")

    return {
        'gs_ground': gs_ground,
        'gs_flight': gs_flight,
        'gs_aleatorio': gs_aleatorio,
        'media_ground': media_ground,
        'std_ground': std_ground,
        'media_flight': media_flight,
        'std_flight': std_flight,
        'media_aleatorio': media_aleatorio,
        'std_aleatorio': std_aleatorio,
        'diferencia': diferencia,
        'resultados_ground': resultados_ground,
        'resultados_flight': resultados_flight,
        'resultados_aleatorio': resultados_aleatorio,
    }


# ============================================================================
# 9. EJECUCIÓN DEL ANÁLISIS COMPLETO
# ============================================================================

if __name__ == "__main__":

    # Rutas de las carpetas de datos (cambiar según donde estén las muestras)
    GROUND_META  = r"/home/ada/guest/ruben_data/Programacion/ground/"
    FLIGHT_META  = r"/home/ada/guest/ruben_data/Programacion/flight/"
    GROUND_ECOLI = r"/home/ada/guest/ruben_data/Programacion/E.coli/genomas_consenso/ground/"
    FLIGHT_ECOLI = r"/home/ada/guest/ruben_data/Programacion/E.coli/genomas_consenso/flight/"

    SEMILLA = 42
    K_META  = [3, 6, 9, 12]
    K_ECOLI = [3, 6, 8, 9, 12]
    random.seed(SEMILLA)
    np.random.seed(SEMILLA)

    # ========================================================================
    # 1. Shannon + k-mers (datos metagenómicos)
    # ========================================================================
    for k in K_META:
        comparacion(GROUND_META, FLIGHT_META, k, semilla=SEMILLA)
        comparacion_shannon(GROUND_META, FLIGHT_META, k, semilla=SEMILLA)
        comparacion_shannon_acumulada(GROUND_META, FLIGHT_META, k, semilla=SEMILLA)

    entropia_normalizada_multik(GROUND_META, FLIGHT_META, K_META,
                                semilla=SEMILLA, correccion=False)
    entropia_normalizada_multik(GROUND_META, FLIGHT_META, K_META,
                                semilla=SEMILLA, correccion=True)

    # ========================================================================
    # 2. K-mers ausentes (datos metagenómicos)
    # ========================================================================
    graficar_kmers_cero_por_genoma(GROUND_META, FLIGHT_META, K_META,
                                   semilla=SEMILLA)
    calcular_missing_kmers_teorico(GROUND_META, FLIGHT_META, K_META,
                                   semilla=SEMILLA)

    # ========================================================================
    # 3. Genome Signature (genomas E. coli consenso)
    # ========================================================================
    gs_carpetas(GROUND_ECOLI, FLIGHT_ECOLI,
                k_min=2, k_max=12, semilla=SEMILLA)

    # ========================================================================
    # 4. Biobit (genomas E. coli consenso)
    # ========================================================================
    biobit_carpetas(GROUND_ECOLI, FLIGHT_ECOLI, k_range=(3, 12))

    # ========================================================================
    # 5. IEPWRMkmer (genomas E. coli consenso)
    #    analisis_multik genera: IEPWRMKmer_coli.png + IEPWRMkmer_3D_Hkmer.png
    #    Luego se exportan CSV y tabla nucleotídica por cada k
    # ========================================================================
    resultados_multi = analisis_multik_IEPWRMkmer(
        GROUND_ECOLI, FLIGHT_ECOLI, K_ECOLI)

    for k in K_ECOLI:
        exportar_resultados_IEPWRMkmer(resultados_multi[k], k)
        frecuencia_nucleotidica_kmers_significativos(
            resultados_multi[k], k, n_top=20)

    # ========================================================================
    # 6. Shannon + k-mers (datos E. coli)
    # ========================================================================
    for k in K_ECOLI:
        comparacion(GROUND_ECOLI, FLIGHT_ECOLI, k, semilla=SEMILLA)
        comparacion_shannon(GROUND_ECOLI, FLIGHT_ECOLI, k, semilla=SEMILLA)
        comparacion_shannon_acumulada(GROUND_ECOLI, FLIGHT_ECOLI, k, semilla=SEMILLA)

    entropia_normalizada_multik(GROUND_ECOLI, FLIGHT_ECOLI, K_ECOLI,
                                semilla=SEMILLA, correccion=False)
    entropia_normalizada_multik(GROUND_ECOLI, FLIGHT_ECOLI, K_ECOLI,
                                semilla=SEMILLA, correccion=True)

    # ========================================================================
    # 7. K-mers ausentes (datos E. coli)
    # ========================================================================
    graficar_kmers_cero_por_genoma(GROUND_ECOLI, FLIGHT_ECOLI, K_ECOLI,
                                   semilla=SEMILLA)
    calcular_missing_kmers_teorico(GROUND_ECOLI, FLIGHT_ECOLI, K_ECOLI,
                                   semilla=SEMILLA)