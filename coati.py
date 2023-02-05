
from Bio import Entrez 
from Bio import SeqIO 
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Align.Applications import ClustalwCommandline
import os
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import matplotlib.pyplot as plt
from Bio import Phylo
##############Primera funcion
def fasta_downloader(filename):
    id_coati=[]
    rep= []
    coati= []
    
    """
    Función para descargar información en formato genbank del NCBI usando los identificadores de accesión almacenados 
    en el archivo 'coati.txt' y guardarlo en la variable 'coati'.
    Además, guarda el resultado en un archivo genbank 'coati.gb' en la carpeta 'data'.
    """
    from Bio import Entrez 
    from Bio import SeqIO 
    with open('C:/Users/Usuario/Documents/GitHub/2022II_GBI6_exam2/data/coati.txt') as coat:
        line = coat.readline()
        for i in coat:
            rep = i.replace('\n','')
            id_coati.append(rep)
        ##Entrar a la web
    Entrez.email = "cristhoplay@gmail.com" 
    with Entrez.efetch( db="nucleotide", rettype="gb", retmode="text", id= id_coati
                  ) as handle: 
        for seq_record in SeqIO.parse(handle, "gb"): 
            coati.append(seq_record)
    SeqIO.write(coati, "data/coati.gb", "genbank")
    return coati
####################
##########Segunda funcion de alineamiento
def alignment(coati):
    """
    Esta función realiza un alineamiento de secuencias utilizando clustalW
    a partir de una variable 'coati' que contiene secuencias en formato genbank.
    El resultado del alineamiento se guarda en los archivos 'coati.aln' y 'coati.dnd'
    en la carpeta 'data'.
    """
    from Bio.Align.Applications import ClustalwCommandline
    from Bio import AlignIO

    # extraer secuencias de la variable 'coati'
    secuencias = [registro.seq for registro in coati]
    
    # crear un archivo temporal para almacenar las secuencias
    with open("data/temp_sequences.fasta", "w") as f:
        for i, seq in enumerate(secuencias):
            f.write(f">seq{i}\n{seq}\n")
    
    # ejecutar clustalW para realizar el alineamiento
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile = "data/temp_sequences.fasta")
    assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing or not found"
    stdout, stderr = clustalw_cline()
    print(clustalw_cline)
    # leer el resultado del alineamiento y guardarlo en archivos .aln y .dnd
    with open("data/temp_sequences.aln","r") as aln: 
        alignment = AlignIO.read(aln,"clustal")
    AlignIO.write(alignment, "data/coati.aln", "fasta")
    AlignIO.write(alignment, "data/coati.dnd", "phylip-relaxed")
    print(alignment)
    return alignment
#############Tercera funcion
def tree(alignment):
    """
    Esta funcion realiza el calculo de las distancias usando un archivo aln, ademas genera un arbol filogenetico 
    y guarda el resultado en un archivo pdf
    """
    
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    import matplotlib.pyplot as plt
    from Bio import Phylo

    
    # Calculo distancias
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    
    # Construct del arbol
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(alignment)
    
    # Mostrar figura
    Phylo.draw(tree)
    plt.show()
    
    # Guardar en pdf
    plt.savefig("data/coati_phylotree.pdf")
  ##########
