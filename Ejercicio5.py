import json
import xml.etree.ElementTree as ET

def cargar_transcripto_desde_archivo(ruta):
    """Carga la secuencia del transcripto desde un archivo de texto"""
    with open(ruta, 'r') as archivo:
        secuencia = archivo.read().strip()
    return secuencia

def generar_primers(secuencia, min_length=18, max_length=24):
    """Genera todos los posibles primers dentro del rango de longitud dado"""
    primers = [secuencia[i:i+length] for i in range(len(secuencia)-min_length+1) 
               for length in range(min_length, max_length+1)]
    return primers

def calcular_gc_content(primer):
    """Calcula el porcentaje de GC en el primer"""
    gc_count = primer.count('G') + primer.count('C')
    return (gc_count / len(primer)) * 100

def evitar_gc_extremos(primer):
    """Evita que el primer tenga GC en los extremos"""
    return primer[0] not in ['G', 'C'] and primer[-1] not in ['G', 'C']

def calcular_tm(primer):
    """Calcula la temperatura de melting (Tm) del primer"""
    gc_count = primer.count('G') + primer.count('C')
    at_count = len(primer) - gc_count
    return 2 * at_count + 4 * gc_count

def filtrar_primers(primers, min_gc=50, max_gc=60, max_tm=67):
    """Filtra los primers según el contenido de GC, extremos, y temperatura de melting"""
    primers_filtrados = []
    for primer in primers:
        gc_content = calcular_gc_content(primer)
        tm = calcular_tm(primer)
        if (min_gc <= gc_content <= max_gc) and evitar_gc_extremos(primer) and tm <= max_tm:
            primers_filtrados.append(primer)
    return primers_filtrados

def cargar_parametros_json(ruta):
    """Carga los parámetros desde un archivo JSON"""
    with open(ruta, 'r') as archivo:
        parametros = json.load(archivo)
    return parametros

# Ejemplo de uso del script:
ruta_transcripto = 'transcripto.txt'  # Ruta al archivo de transcripto
ruta_parametros = 'parametros.json'   # Ruta al archivo JSON con los parámetros

# Cargar la secuencia y parámetros
secuencia = cargar_transcripto_desde_archivo(ruta_transcripto)
parametros = cargar_parametros_json(ruta_parametros)

# Generar primers
primers = generar_primers(secuencia, parametros['min_length'], parametros['max_length'])

# Filtrar primers según los criterios
primers_filtrados = filtrar_primers(primers, parametros['min_gc'], parametros['max_gc'], parametros['max_tm'])

# Guardar resultados
print("Primers filtrados que cumplen los criterios:")
for primer in primers_filtrados:
    print(primer)
