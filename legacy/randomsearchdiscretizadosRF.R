# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: randomsearchdiscretizadosRF.R (dated 2019-06-04 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: 01c360dcc3eece4da846108df76116a99dc856653d68593482cfaf0c5dd2ca23
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(randomForest)
library(caret)
library(ROCR)
library(resample)
library(plyr)
library(caTools)

##Leemos el archivo.rds que creamos anteriormente y representa la matriz de expresion:
data <- readRDS("finalmatrix.rds")
quantile(data$AUC)

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
##Preparamos los valores AUC para indicar cuando son sensibles o resistentes
##a un farmaco las lineas celulares de estudio, es decir, discretizarlos:
tercercuantil <- Quantile(gen.var.AUC$class, 0.75)
gen.var.AUC$class[gen.var.AUC$class<tercercuantil] <- 0
gen.var.AUC$class[gen.var.AUC$class>=tercercuantil] <- 1
gen.var.AUC$class <- factor(gen.var.AUC$class, levels = c("0", "1"), labels= c("Resistente", "Sensible"))

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
                   metric = "Accuracy", 
                   tuneLength = 20, trControl = control)
print(rf_random)

##Ajuste del modelo de Random Forest:
fit_rf10 <- randomForest(class~.,
                         data = data_train,
                         mtry= 43)
##Preccion del modelo ajustado con los datos de prueba o testSet y creacion de la matriz
##de confusion y curva Roc:
prediction <-predict(fit_rf10, data_test)
confusionMatrix(prediction, data_test_class)
probs <- predict(fit_rf10, data_test, type= "prob")
pred <- prediction(probs[,2] , data_test_class)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(perf)

##Guardado del modelo ajustado de random Forest:##Guardado del modelo ajustado de random Forest:
save(fit_rf10, file = "fit_randomsearch_RF_erlotinib_8020_discretizado.RData")

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
                   metric = "Accuracy", 
                   tuneLength = 20)
print(rf_random602020_validation)

##Ajuste del modelo de Random Forest con los datos de entrenamiento:
rf_random602020 <- randomForest(class ~., 
                   data = remain80, 
                   mtry= 2)

##Preccion del modelo ajustado con los datos de prueba o testSet y creacion de la matriz
##de confusion y curva Roc:
prediction602020 <-predict(rf_random602020, testing20)
confusionMatrix(prediction602020, class_testing20)
probs1 <- predict(rf_random602020, testing20, type= "prob")
pred1 <- prediction(probs1[,2] , class_testing20)
perf1 <- performance(pred1, measure = "tpr", x.measure = "fpr")
plot(perf1)

##Guardado del modelo ajustado de random Forest:##Guardado del modelo ajustado de random Forest:
save(rf_random602020, file = "fit_randomsearch_RF_erlotinib_602020discretizado.RData")