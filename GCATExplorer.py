import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF
from Bio import SeqIO
import re


# Leer archivo de secuencias
file_path = "sequence_moth.txt"
# Procesar el archivo complementario para obtener Organism
complementary_file_path = 'nuccore_result.txt'

output_pdf_path = "informe_secuencias_intento_mod2656.pdf"



# Función para identificar outliers y marcarlos como 1 o 0
def find_outliers_binary(group):
    Q1 = group.quantile(0.25)
    Q3 = group.quantile(0.75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    return group.apply(lambda x: 1 if x < lower_bound or x > upper_bound else 0)

# Función para reemplazar números romanos en genes
def replace_roman_with_arabic(name):
    roman_to_arabic_map = {
        "I": "1",
        "II": "2",
        "III": "3",
        "IV": "4",
        "V": "5",
        "VI": "6",
        "VII": "7",
        "VIII": "8",
        "IX": "9",
        "X": "10"
    }
    match = re.search(r"(I|II|III|IV|V|VI|VII|VIII|IX|X)$", name)
    if match:
        roman = match.group(0)
        if roman in roman_to_arabic_map:
            name = name.replace(roman, roman_to_arabic_map[roman])
    return name

# Función para calcular el contenido GC
def gc_content(sequence):
    if len(sequence) == 0:
        return 0
    sequence = sequence.upper()
    G_count = sequence.count("G")
    C_count = sequence.count("C")
    GC_content = (G_count + C_count) / len(sequence)
    return round(100 * GC_content, 2)


# Funciones para calcular métricas individuales
def base_count(sequence, base):
    """Cuenta la frecuencia de una base específica en la secuencia."""
    sequence = sequence.upper()
    return sequence.count(base)

def at_content(sequence):
    """Calcula el contenido AT (%) de la secuencia."""
    if len(sequence) == 0:
        return 0
    sequence = sequence.upper()
    A_count = sequence.count("A")
    T_count = sequence.count("T")
    AT_content = (A_count + T_count) / len(sequence)
    return round(100 * AT_content, 2)

data = []

for record in SeqIO.parse(file_path, "fasta"):
    description = record.description
    accession = description.split("|")[1].split("_cds")[0]
    gene = "Unknown"
    name_protein = "Unknown"

    if "[gene=" in description:
        gene = description.split("[gene=")[1].split("]")[0].upper()
        gene = replace_roman_with_arabic(gene)

    if "[protein=" in description:
        name_protein = description.split("[protein=")[1].split("]")[0].upper()

    sequence = str(record.seq)
    data.append((accession, gene, name_protein, sequence))

# Crear un DataFrame con los datos
df = pd.DataFrame(data, columns=["Accession", "Gene", "Name_Protein", "Sequence"])
df["GC_Content"] = df["Sequence"].apply(gc_content)
df["AT_Content"] = df["Sequence"].apply(at_content)

data_complementary = []

with open(complementary_file_path, "r") as file:
    lines = file.readlines()
    previous_lines = []

    for line in lines:
        line = line.strip()
        if line.startswith("NC_"):
            if len(previous_lines) >= 2:
                potential_organism_line = previous_lines[-2]
                tokens = potential_organism_line.split()
                organism = " ".join(tokens[1:3]) if len(tokens) >= 3 else "Desconocido"
                accession = line.split()[0]
                data_complementary.append({"Accession": accession, "Organism": organism})
        previous_lines.append(line)
        if len(previous_lines) > 2:
            previous_lines.pop(0)

df_complementary = pd.DataFrame(data_complementary)
df = pd.merge(df, df_complementary, on="Accession", how="left")

# Reestructurar Gene basado en Name_Protein
def restructure_gene(df):
    first_genes = df.groupby("Name_Protein")['Gene'].transform(lambda x: x.iloc[0].upper())
    df['Gene'] = df.apply(lambda row: [first_genes[row.name], row['Gene'].upper()], axis=1)
    return df

df = restructure_gene(df)

# Calcular longitud de las secuencias
df["Sequence_Length"] = df["Sequence"].apply(len)

# Calcular outliers para GC_Content
df["Outlier_GC_Content"] = (
    df.groupby(df["Gene"].apply(lambda x: x[0]))["GC_Content"]
    .transform(lambda g: find_outliers_binary(g))
)

# Calcular outliers para AT_Content
df["Outlier_AT_Content"] = (
    df.groupby(df["Gene"].apply(lambda x: x[0]))["AT_Content"]
    .transform(lambda g: find_outliers_binary(g))
)

# Calcular outliers para la longitud de las secuencias
df["Outlier_Length"] = (
    df.groupby(df["Gene"].apply(lambda x: x[0]))["Sequence_Length"]
    .transform(lambda g: find_outliers_binary(g))
)

df["A_Count"] = df["Sequence"].apply(lambda seq: base_count(seq, "A"))
df["T_Count"] = df["Sequence"].apply(lambda seq: base_count(seq, "T"))
df["G_Count"] = df["Sequence"].apply(lambda seq: base_count(seq, "G"))
df["C_Count"] = df["Sequence"].apply(lambda seq: base_count(seq, "C"))

# Consolidar toda la información en un único DataFrame
Secuencia_dataFrame = df[[
    "Gene",
    "Accession",
    "Name_Protein",
    "GC_Content",
    "Outlier_GC_Content",
    "AT_Content",
    "Outlier_AT_Content",
    "Sequence_Length",
    "Outlier_Length",
    "A_Count",
    "T_Count",
    "G_Count",
    "C_Count",
    "Organism"
]]

# Imprimir el DataFrame consolidado para inspección
print(Secuencia_dataFrame)

# Crear el informe en PDF
class PDF(FPDF):
    def header(self):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, 'Informe de Análisis de Secuencias', 0, 1, 'C')

    def chapter_title(self, title):
        self.set_font('Arial', 'B', 12)
        self.cell(0, 10, title, 0, 1, 'L')
        self.ln(5)

    def chapter_body(self, body):
        self.set_font('Arial', '', 12)
        self.multi_cell(0, 10, body)
        self.ln()

pdf = PDF()
pdf.add_page()
pdf.chapter_title("Introducción")
pdf.chapter_body(
    "Este informe analiza el contenido GC de las secuencias agrupadas por gen y organiza los resultados con base en su acceso y organismo.")

# Barplot del número de secuencias
pdf.add_page()
pdf.chapter_title("Número de Secuencias por Gen")
gen_counts = df["Gene"].apply(lambda x: x[0]).value_counts().reset_index()
gen_counts.columns = ["Gene", "Count"]

# Calcular límites del eje y
min_count = gen_counts["Count"].min()
max_count = gen_counts["Count"].max()
lower_limit = (min_count * 0.85) // 10 * 10
upper_limit = (max_count * 1.15 + 9) // 10 * 10

plt.figure(figsize=(12, 6))
sns.barplot(data=gen_counts, x="Gene", y="Count", palette="tab10")
plt.title("Número de Secuencias por Gen", fontsize=14)
plt.xlabel("Gen", fontsize=12)
plt.ylabel("Cantidad de Secuencias", fontsize=12)
plt.xticks(rotation=45)
plt.ylim(lower_limit, upper_limit)
plt.tight_layout()
barplot_path = "barplot_gen_counts.png"
plt.savefig(barplot_path)
plt.close()
pdf.image(barplot_path, x=10, y=None, w=190)

# Extraer solo el primer elemento de Gene y generar una lista única
unique_genes = Secuencia_dataFrame.copy()
unique_genes["Gene"] = unique_genes["Gene"].apply(lambda x: x[0] if isinstance(x, list) else x)

# Crear DataFrame con Genes únicos y sus proteínas asociadas
gene_protein_pairs = unique_genes[["Gene", "Name_Protein"]].drop_duplicates()

# Iterar sobre los pares únicos para agregarlos al PDF
for _, row in gene_protein_pairs.iterrows():
    pdf.set_font("Arial", size=10)
    pdf.cell(0, 10, f"Gene: {row['Gene']}, Proteina: {row['Name_Protein']}", ln=True)

# Boxplot del %GC
plt.figure(figsize=(12, 6))
sns.boxplot(data=df, x=df["Gene"].apply(lambda x: x[0]), y="GC_Content", palette="tab10", dodge=False)
plt.title("Distribución del %GC por Gen", fontsize=14)
plt.xlabel("Gen", fontsize=12)
plt.ylabel("%GC", fontsize=12)
plt.xticks(rotation=90)
plt.tight_layout()
boxplot_path = "boxplot_gc.png"
plt.savefig(boxplot_path)
plt.close()
pdf.chapter_title("Distribución del %GC por Gen")
pdf.image(boxplot_path, x=10, y=None, w=190)

# Boxplot del %AT
plt.figure(figsize=(12, 6))
sns.boxplot(data=df, x=df["Gene"].apply(lambda x: x[0]), y="AT_Content", palette="tab10", dodge=False)
plt.title("Distribución del %AT por Gen", fontsize=14)
plt.xlabel("Gen", fontsize=12)
plt.ylabel("%AT", fontsize=12)
plt.xticks(rotation=90)
plt.tight_layout()
at_boxplot_path = "boxplot_at.png"
plt.savefig(at_boxplot_path)
plt.close()
pdf.chapter_title("Distribución del %AT por Gen")
pdf.image(at_boxplot_path, x=10, y=None, w=190)

# Boxplot de la longitud de las secuencias
plt.figure(figsize=(12, 6))
sns.boxplot(data=df, x=df["Gene"].apply(lambda x: x[0]), y="Sequence_Length", palette="tab10", dodge=False)
plt.title("Distribución de la Longitud de las Secuencias por Gen", fontsize=14)
plt.xlabel("Gen", fontsize=12)
plt.ylabel("Longitud de la Secuencia (bp)", fontsize=12)
plt.xticks(rotation=90)
plt.tight_layout()
sequence_length_boxplot_path = "boxplot_sequence_length.png"
plt.savefig(sequence_length_boxplot_path)
plt.close()
pdf.chapter_title("Distribución de la Longitud de las Secuencias por Gen")
pdf.image(sequence_length_boxplot_path, x=10, y=None, w=190)

# Agregar una página en blanco
pdf.add_page()

# Filtrar Accessions con outliers en AT_Content, GC_Content y Length
outlier_accessions = df[
    (df["Outlier_AT_Content"] == 1) &
    (df["Outlier_GC_Content"] == 1) &
    (df["Outlier_Length"] == 1)
    ]["Accession"].unique()

# Gráfico de barras para la cantidad de nucleótidos por gen y por accession filtrados
# Inicializar la posición para colocar las imágenes en el PDF
x_positions = [10, 105]  # Coordenadas x para dos columnas
y_position = 30  # Coordenada y inicial (más abajo para evitar tapar el título)
image_width = 90  # Ancho de cada imagen
row_height = 95  # Altura de cada fila

for accession in outlier_accessions:
    sub_df = df[df["Accession"] == accession].copy()  # Evitar SettingWithCopyWarning
    sub_df.loc[:, "Gene"] = sub_df["Gene"].apply(
        lambda x: x[0] if isinstance(x, list) else x)  # Usar loc para modificaciones

    organism = sub_df["Organism"].iloc[0]

    melt_df = sub_df.melt(id_vars=["Gene"],
                          value_vars=["A_Count", "T_Count", "G_Count", "C_Count"],
                          var_name="Nucleotide",
                          value_name="Count")

    plt.figure(figsize=(14, 8))
    sns.barplot(data=melt_df, x="Gene", y="Count", hue="Nucleotide", palette="Set2")
    plt.title(f"Cantidad de Nucleótidos por Gen para {organism} {accession}", fontsize=12)
    plt.xlabel("Gen", fontsize=10)
    plt.ylabel("Cantidad de Nucleótidos", fontsize=10)
    plt.xticks(rotation=45)
    plt.tight_layout()

    # Configurar leyenda con menos elementos si es necesario
    handles, labels = plt.gca().get_legend_handles_labels()
    if len(labels) > 10:
        plt.legend(handles[:10], labels[:10], title="Nucleótido", loc="upper right", bbox_to_anchor=(1.05, 1))
    else:
        plt.legend(title="Nucleótido", loc="upper right", bbox_to_anchor=(1.05, 1))

    nucleotide_barplot_path = f"nucleotide_barplot_{accession}.png"
    try:
        plt.savefig(nucleotide_barplot_path)
        #print(f"Guardada correctamente: {nucleotide_barplot_path}")
    except Exception as e:
        print(f"Error al guardar la gráfica: {e}")
    finally:
        plt.close()

    # Agregar la gráfica al PDF
    pdf.image(nucleotide_barplot_path, x=x_positions[0], y=y_position, w=image_width)
    x_positions = x_positions[::-1]  # Alternar columnas
    if x_positions[0] == 10:  # Si volvemos a la primera columna
        y_position += row_height  # Ajustar coordenada y para la próxima fila
    if y_position > 250:  # Si la posición y excede el límite, agregar nueva página
        pdf.add_page()
        y_position = 30  # Reiniciar posición para la nueva página

pdf.output(output_pdf_path)
print(f"Informe generado: {output_pdf_path}")