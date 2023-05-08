# Import modules
import pandas as pd
import openpyxl

# Disable warning when replacing "GO:" with an empty string
pd.options.mode.chained_assignment = None  # default='warn'.

# Convert the excel file into a dictionary where the keys are the sheet names and the values are the resulting dataframes
df = pd.read_excel("EggNog_v261_proteome_HopBA1_RNAseq.xlsx",
                   sheet_name = "EggNog v261 proteome",
                   header=2,
                   usecols=["query", "GOs"])

# Exclude rows without GO terms
df = df[df.GOs != "-"]

# Delete "GO:" string, keeping only the GO term number
df['GOs'] = df['GOs'].str.replace("GO:", "")

# Iterate through the dictionary and write a new excel file containing only the gene IDs and associated GO terms for each sheet from the original spreadsheet
df.to_excel("Extracted_EggNog_v261_IDs_GOs.xlsx")
print("HERE")
input()

# Load the previously generated workbook, get its sheet names and start a sheet counter
wb = openpyxl.load_workbook("Extracted_EggNog_v261_IDs_GOs.xlsx")
ws = wb["Sheet1"]

# Iterate through the rows returning row objects
rows = ws.iter_rows(min_row = 2, min_col = 2) # Desired data starts from row and column 2

# Append the IDs to a dictionary as keys and the GO terms as values after splitting the latter by commas.
ID_GO_Dict = {}
for ID,GOs in rows:
    ID_GO_Dict[ID.value] = str(GOs.value).split(",")

# Create a BinGO annotation file for UP genes in leaves specifying the species name, type of GO terms and the BinGO curator in the first line
with open("Nb_GO_Terms_EggNog_v261_HopBA1_RNAseq.txt", "w") as f:
    print("(species=Nicotiana benthamiana)(type=Biological Process)(curator=GO)", file=f)    

ID_GO_lst = list()

# Iterate through the dictionary to write to the BinGO annotation file Nb's gene IDs and their associated GO terms
for ID,GOs in ID_GO_Dict.items():
    for GO in GOs:
        if ID and GO is not None:
            ID_GO_lst.append(ID + " = " + str(GO))

for element in ID_GO_lst:
    with open("Nb_GO_Terms_EggNog_v261_HopBA1_RNAseq.txt","a") as f:
        print(element, file=f)

print("####Done!####")
