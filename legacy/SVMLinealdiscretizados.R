# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: SVMLinealdiscretizados.R (dated 2019-06-04 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: 521c57a09f770fcbf47e827f52b0952c67ec2661690f207a18aa696db3b6ea11
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(ROCR)
library(caret)
library(e1071)
library(resample)
library(plyr)
library(caTools)

#Leemos el archivo.rds que creamos anteriormente y representa la matriz de expresion
##y calculamos los cuartiles para discretizar posteriormente los valores predictores:
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

##Para confirmar que la proporcion de resistentes y sensibles es aproximadamente
##la misma en ambos conjuntos de datos, vemos la distribucion de los tipos en cada set:
table(data_train$class)
table(data_test_class)



##Busqueda de hiperparámetros:
set.seed(1)
tune.out = tune(svm,
                class ~ .,
                data = data_train,
                kernel = "linear",
                ranges = list(cost = c(0.1, 1, 10, 100, 1000)))
summary (tune.out)

##Ajuste del model SVM lineal con los hiperparámetros
##óptimos seleccinados en la búsqueda anterior en 10-fold Cross-Validation:
svmfit.opt <- svm(class ~.,
                 data = data_train,
                 kernel = "linear",
                 family = "binomial",
                 cost = 0.1,
                 cross=10,
                 decision.values = T)


##Prediccion con el modelo ajustado de los 
##datos reservados como testSet y cracion de la matriz de confusion.
fitted2 <- predict(svmfit.opt, data_test, type= "prob")
confusionMatrix(fitted2, data_test_class)

####Prediccion con el modelo ajustado de los 
##datos reservados como testSet y creacion de la curva ROC:
rocplot <- function(pred, truth, ...) {
  predob <- prediction(pred, truth)
  perf <- performance(predob, "tpr", "fpr")
  plot(perf, ...)
}
fitted <- attributes(predict(svmfit.opt, data_test, type= "prob",
                             decision.values =T))$decision.values
rocplot(fitted,
        data_test_class)

##Guardado como archivo .RData el modelo ajustado:
save(svmfit.opt, file = "fit_SVMlineal_discretizado_erlotinib_8020.RData")


#Partición 60/20/20 (3:1:1 training, validation, testing)#####
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

##Busqueda de hiperparámetros con el set de validacion:
set.seed(1)
tune.out = tune(svm,
                class ~ .,
                data = validation20,
                kernel = "linear",
                ranges = list(cost = c(0.1, 1, 10, 100, 1000)))
summary (tune.out)

##Ajuste del model SVM lineal con los hiperparámetros
##óptimos seleccinados en la búsqueda anterior con los el set de entrenamiento:
svmfit.opt <- svm(class ~.,
                  data = remain80,
                  kernel = "linear",
                  family = "binomial",
                  cost = 0.1,
                  cross=10,
                  decision.values = T)



##Prediccion con el modelo ajustado de los 
##datos reservados como testSet y cracion de la matriz de confusion.
fitted2 <- predict(svmfit.opt, testing20, type= "prob")
confusionMatrix(fitted2, class_testing20)

####Prediccion con el modelo ajustado de los 
##datos reservados como testSet y creacion de la curva ROC:
rocplot <- function(pred, truth, ...) {
  predob <- prediction(pred, truth)
  perf <- performance(predob, "tpr", "fpr")
  plot(perf, ...)
}
fitted <- attributes(predict(svmfit.opt, testing20, type= "prob",
                             decision.values =T))$decision.values
rocplot(fitted,
        class_testing20)
save(svmfit.opt, file = "fit_SVMlineal_discretizado_erlotinib_602020.RData")