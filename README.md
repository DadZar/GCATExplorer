# GCATExplorer

## Descripción General
Este script está diseñado para analizar secuencias genéticas, extrayendo características clave como el contenido GC, contenido AT, longitud de las secuencias e identificando valores atípicos basados en métricas estadísticas. Además, genera visualizaciones y un informe consolidado en formato PDF.

## Funcionalidades

1. **Análisis de Secuencias**: 
   - Procesa secuencias desde un archivo FASTA.
   - Extrae metadatos como números de acceso, nombres de genes, nombres de proteínas e información de organismos.

2. **Cálculo de Contenido GC y AT**:
   - Calcula el porcentaje de contenido GC y AT para cada secuencia.

3. **Detección de Valores Atípicos**:
   - Identifica valores atípicos para el contenido GC, contenido AT y longitud de las secuencias usando el método del Rango Intercuartílico (IQR).


4. **Visualizaciones**:
   - Genera gráficos de barras, diagramas de caja y mapas de calor para la distribución de secuencias y métricas.

5. **Generación de Informes en PDF**:
   - Consolida los resultados en un informe estructurado en formato PDF.
   - Incluye gráficos y datos resumidos.

## Dependencias
El script requiere las siguientes bibliotecas de Python:
- `pandas`
- `matplotlib`
- `seaborn`
- `fpdf`
- `BioPython`
- `re`

Asegúrese de que todas las dependencias estén instaladas antes de ejecutar el script. Use `pip install <paquete>` para instalar cualquier biblioteca faltante.

## Archivos de Entrada
1. **Archivo FASTA**: Contiene las secuencias genéticas y metadatos.
   - Ejemplo: `sequence_insects.txt` o `sequence_moth.txt`
2. **Archivo Complementario de Metadatos**: Contiene información sobre organismos y genomas.
   - Ejemplo: `nuccore_result.txt`

## Archivos de Salida
1. **Informe en PDF**: 
   - Incluye hallazgos resumidos y visualizaciones.
   - Ejemplo: `informe_secuencias_intento_mod26.pdf`


## Instrucciones de Uso
1. Coloque los archivos de entrada (`sequence_insects.txt`, `nuccore_result.txt`) en el directorio del script.
2. Ejecute el script utilizando:
   ```bash
   python <nombre_del_script>.py
   ```
3. Los archivos de salida se generarán en el mismo directorio.

## Funciones Clave

### Procesamiento y Preprocesamiento
- `gc_content(sequence)`: Calcula el contenido GC como porcentaje.
- `at_content(sequence)`: Calcula el contenido AT como porcentaje.
- `replace_roman_with_arabic(name)`: Convierte números romanos en nombres de genes a números arábigos.

### Análisis Estadístico
- `find_outliers_binary(group)`: Identifica valores atípicos en un conjunto de datos usando el método IQR.

### Reestructuración de Datos
- `restructure_gene(df)`: Ajusta los nombres de genes basándose en los nombres de proteínas.

### Visualización
- Gráficos de barras y diagramas de caja para métricas específicas de genes.
- Mapas de calor para la distribución del contenido GC entre genes y accesiones.

### Generación de Informes
- Utiliza la biblioteca `FPDF` para compilar datos y visualizaciones en un informe PDF.

## Estructura de Archivos
```
|-- directorio_del_script/
    |-- sequence_insects.txt        # Archivo FASTA de entrada
    |-- nuccore_result.txt          # Archivo de metadatos de entrada
    |-- informe_pdf_salida.pdf      # Informe PDF generado
```

## Notas
- Asegúrese de que los archivos de entrada estén correctamente formateados.
- Ajuste las rutas de los archivos en el script si es necesario.
