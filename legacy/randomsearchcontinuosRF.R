# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: randomsearchcontinuosRF.R (dated 2019-06-06 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: 53c68a53637c89824b8317cc9a66f95b1964b04b409db24ed5da3fae633c9924
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(randomForest)
library(caret)
library(miscTools)
library(resample)
library(plyr)
library(caTools)
library(ggplot2)

##Leemos el archivo.rds que creamos anteriormente y representa la matriz de expresion:
data <- readRDS("finalmatrix.rds")

##Creacion de subsets dependiendo de la Varianza:
##Calculamos las varianzas de cada columna del dataframe y las ordenamos en orden decreciente
##para seleccionar mas tarde los X genes de expresion con varianzas mas altas.
n <- 200
data.var <- sort(colVars(data, na.rm = TRUE), decreasing = TRUE) 
data.var.n <- data.var[1:n]


##Seleccionamos los nombres de los X genes con las varianzas mas altas:
vectornombres <- names(data.var.n)

##Filtramos el dataframe original que hemos llamado "data" con las X columnas pertenecientes
## a los genes que hemos filtrado por sus varianzas.
gen.var.n <- subset(data, select = vectornombres)

##Añadimos la columna de AUC que no habiamos seleccionado previamente pues no corresponde
##a ningun gen.
gen.var.AUC <- data.frame(gen.var.n, data$AUC)
gen.var.AUC = rename(gen.var.AUC, c(data.AUC= "class"))




#Particion 80/20:
##Tomamos un 80% de los datos como conjunto de entrenamiento y el 20% restante
##como conjunto de test:
set.seed(1234)
split <- sample.split(gen.var.AUC$class, SplitRatio = 0.80)
data_train <- subset(gen.var.AUC, split==TRUE)
data_test_total <- subset(gen.var.AUC, split==FALSE)
data_test <- data_test_total[,1:n]
data_test_class <- data_test_total$class

#Busqueda de hiperparametros con Random Search:

control <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 5)
set.seed(1)
rf_random <- train(class ~., 
                   data = data_train, 
                   method = "rf", 
                   metric = "RMSE", 
                   tuneLength = 20, trControl = control)

print(rf_random)

##Ajuste del modelo de Random Forest:
rf <- randomForest(class ~ ., data=data_train, mtry=2)
print(rf)
##Preccion del modelo ajustado con los datos de prueba o testSet:
prediction <-predict(rf, data_test)

##Calculo de la R^2:
rsq <- function (x, y) cor(x, y) ^ 2
r2 <- rsq(data_test_class, data_test_class - prediction)


##Creacion de la gráfica R^2:
p <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=data_test_class, pred=prediction))
p + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("RandomForest Regression in R r^2=", r2, sep=""))

##Guardado del modelo ajustado de random Forest:##Guardado del modelo ajustado de random Forest:
save(rf, file= "rf_ramdomsearch_RF_erlotinib_8020_continuo.RData")

#60/20/20 (3:1:1 training, validation, testing)#####
##Division de los datos en 80% y 20%:
samp <- sample(1:nrow(gen.var.AUC),.8*nrow(gen.var.AUC))
remain80 <- gen.var.AUC[samp,] ##  80% 
testing20 <- gen.var.AUC[-samp,]  ## 20%
class_testing20 <- testing20[,"class"]

##Division del 80% en train y validation (3:1)
samp2 <- sample(1:nrow(remain80),.75*nrow(remain80))
train60 <- remain80[samp2,] ## 60%
validation20 <- remain80[-samp2,] ## 20%

##Comprobacion de las particiones:
Reduce("intersect",list(train60,validation20,testing20))


#Busqueda de hiperparametros con Random Search:
set.seed(1)
rf_random602020_validation <- train(class ~., 
                   data = validation20, 
                   method = "rf", 
                   metric = "RMSE", 
                   tuneLength = 20)
print(rf_random602020_validation)

##Ajuste del modelo de Random Forest con los datos de entrenamiento:
rf_random602020 <- randomForest(class ~ ., data=remain80, mtry=2)

##Preccion del modelo ajustado con los datos de prueba o testSet:
prediction602020 <-predict(rf_random602020, testing20)

##Calculo de la R^2:
rsq <- function (x, y) cor(x, y) ^ 2
r2_602020 <- rsq(class_testing20, class_testing20 - prediction602020)
##Creacion de la gráfica R^2:
p602020 <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=class_testing20, pred=prediction602020))
p602020 + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("RandomForest Regression in R r^2=", r2_602020, sep=""))

##Guardado del modelo ajustado de random Forest:
save(rf_random602020, file= "rf_ramdomsearch_RF_erlotinib_602020_continuo.RData")
