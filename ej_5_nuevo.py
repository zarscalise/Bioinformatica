import json
import xml.etree.ElementTree as ET

def cargar_transcripto_desde_archivo(ruta):
    """Carga la secuencia del transcripto desde un archivo de texto"""
    with open(ruta, 'r') as archivo:
        secuencia = archivo.read().strip()
    return secuencia

def generar_primers(secuencia, min_length=18, max_length=24):
    """Genera todos los posibles primers dentro del rango de longitud dado"""
    primers = [secuencia[i:i+length] for i in range(len(secuencia) - min_length + 1) 
               for length in range(min_length, max_length + 1)]
    return primers

def calcular_gc_content(primer):
    """Calcula el porcentaje de GC en el primer"""
    if len(primer) == 0:  # Verificación para evitar división por cero
        return 0
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

def filtrar_primers(primers, min_gc=50, max_gc=60, max_tm=67, max_primers=5):
    """Filtra y selecciona hasta 5 primers que cumplan los criterios de GC, extremos y Tm"""
    primers_filtrados = []
    for primer in primers:
        if len(primer) == 0:  # Evita procesar primers vacíos
            continue
        gc_content = calcular_gc_content(primer)
        tm = calcular_tm(primer)
        if (min_gc <= gc_content <= max_gc) and evitar_gc_extremos(primer) and tm <= max_tm:
            primers_filtrados.append(primer)
            if len(primers_filtrados) == max_primers:
                break
    return primers_filtrados or []  # Asegura que siempre retorne una lista

def cargar_parametros_json(ruta):
    """Carga los parámetros desde un archivo JSON"""
    with open(ruta, 'r') as archivo:
        parametros = json.load(archivo)
    return {
        'min_length': parametros.get('min_length', 18),
        'max_length': parametros.get('max_length', 24),
        'min_gc': float(parametros.get('min_gc', 50)),
        'max_gc': float(parametros.get('max_gc', 60)),
        'max_tm': float(parametros.get('max_tm', 67))
    }

def cargar_parametros_xml(ruta):
    """Carga los parámetros desde un archivo XML"""
    tree = ET.parse(ruta)
    root = tree.getroot()
    return {
        'min_length': int(root.find('min_length').text),
        'max_length': int(root.find('max_length').text),
        'min_gc': float(root.find('min_gc').text),
        'max_gc': float(root.find('max_gc').text),
        'max_tm': float(root.find('max_tm').text)
    }


# Ejemplo de uso del script:
ruta_transcripto = "C:\\Users\\VALENTINA\\Documents\\Bioinfo\\transcript.txt"  # Ruta al archivo de transcripto
ruta_parametros = "C:\\Users\\VALENTINA\\Documents\\Bioinfo\\parametros.json"  # Ruta al archivo de configuración (JSON o XML)


# Cargar la secuencia y parámetros
secuencia = cargar_transcripto_desde_archivo(ruta_transcripto)
print(secuencia)
try:
    parametros = cargar_parametros_json(ruta_parametros)
except json.JSONDecodeError:
    parametros = cargar_parametros_xml(ruta_parametros)

# Generar primers
primers = generar_primers(secuencia, parametros['min_length'], parametros['max_length'])
#print("Primers generados:", primers)

# Filtrar primers y seleccionar los primeros cinco
primers_filtrados = filtrar_primers(primers, parametros['min_gc'], parametros['max_gc'], parametros['max_tm'])

# Guardar resultados
if primers_filtrados:
    print("Primers filtrados que cumplen los criterios:")
    for primer in primers_filtrados:
        print(primer)
else:
    print("No se encontraron primers que cumplan los criterios.")