# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: creacionmatriz.R (dated 2019-06-04 in the deposited annex),
# original encoding ASCII, CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: e3bb123609c2a34cd0182a49191e4864acd72277285a7ef5a3b28e4657bfd941
# ----------------------------------------------------------------------------
#Cargamos las librerias que vamos a usar:
library(readxl)
library(data.table)

#Creacion de la matriz de respuesta a farmacos:

##Cargamos los datos de respuesta a farmacos:
response <-read_excel("v17.3_fitted_dose_response.xlsx")
##Seleccionamos las columnas de interes: nombre del farmaco, linea celular y valor AUC.
response1 <- data.frame(DRUG_NAME=response$DRUG_NAME, CELL_LINE_NAME=response$CELL_LINE_NAME, AUC=response$AUC)

##Bucle para creacion de matriz de farmaco en el que mencionando el farmaco deseado nos 
##indica si se encuentra en nuestro datos, y si es asi, devuelve un output con los valores AUC asociados a la
##linea celular tumoral en la que actua:
farmaco = 'Erlotinib'
if(farmaco %in% response1[,1]) {
  matriz.drug.celline.auc <- response1[response1$DRUG_NAME ==farmaco,]
  matriz.celline.auc <-  matriz.drug.celline.auc[,-1]
  saveRDS(matriz.celline.auc, file= 'farmaco.rds')
}else{
  print("te lo has inventado o no existe en este database")
} 

dim(matriz.celline.auc)
sapply( matriz.celline.auc, function(x) sum(is.na(x)))

#Creacion de la matriz de expresion:
##Leemos los metadatos asociados a los datos de expresion:
metadata <- read.csv("DepMap-2019q1-celllines.csv")

##Seleccionamos las columnas de metadatos de interes y los combinamos
## con la funcion merge() con las columnas del dataframe 'response' que tengan en comun 
##para obtener una relacion comun de lenguaje de analisis entre la matriz de farmacos y la de expresion: 
response2 <-data.frame (CELL_LINE_NAME=response$CELL_LINE_NAME, COSMIC_ID=response$COSMIC_ID)
metadata1 <- data.frame(DepMap_ID=metadata$DepMap_ID, COSMIC_ID=metadata$COSMIC_ID)

##Obtenemos una asociacion comun entre la nomenclatura 'DEPMAP_ID' de los datos de expresion y 'Cell_line_name'
##de los datos de respuesta a farmacos a traves de 'COSMIC_ID'que ambos datos tienen en comun:
rm <- merge(response2,metadata1)

##Nos quedamos con las columnas que vamos a necesitar elimando la columna 'COSMIC_ID' que solo hemos necesitado para
##asociar las otras dos y que solo guarde los valores unicos quitando cualquier repeticion:
rm2 <- unique(data.frame(CELL_LINE_NAME=rm$CELL_LINE_NAME, DepMap_ID=rm$DepMap_ID))

##Cargamos los datos de expresion con la funcion 'fread()' para agilizar el proceso de lectura ya que el archivo es relativamente
##grande:
expresion <-fread("CCLE_depMap_19Q1_TPM.csv")

##Nombramos la primera columna con el mismo nombre que en los metadatos para posteriormente relacionar los metadatos
##filtrados en la variable 'c2'con los datos de expresion:
names(expresion)[1]="DepMap_ID"
rmexp <- merge(rm2, expresion, by= "DepMap_ID")

##Eliminamos la columna 'DepMap_ID'pues ya tenemos la relacion que necesitabamos entre los nombres de las lineas celulares
##de los datos de respuesta a farmacos y los datos de expresion 'DepMap_ID'y lo guardamos (por si acaso para comprobaciones):
rmexp1 <- rmexp[,-1]
saveRDS(rmexp1, file = 'expresion.rds')

##Fusionamos los datos obtenidos con los valores AUC asociados, 
##y usamos la columna de lineas celulares para nombrar las filas de la matriz:
finalmatrix <- merge(rmexp1, matriz.celline.auc, all = FALSE)
rownames(finalmatrix) <- finalmatrix$CELL_LINE_NAME
finalmatrix1 <- finalmatrix[,-1]


##Revisamos que se hayan eliminado:
sapply(finalmatrix1, function(x) sum(is.na(x)))
dim(finalmatrix1)
##Guardamos la matriz de expresion:
saveRDS(finalmatrix1, file = 'finalmatrix.rds')

