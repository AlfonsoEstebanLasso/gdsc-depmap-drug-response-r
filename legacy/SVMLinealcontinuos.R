# ----------------------------------------------------------------------------
# Legacy script from the MSc thesis annex (June 2019).
# Thesis: "Predicción de respuesta a fármacos quimioterapéuticos a partir de
# datos genómicos", MSc in Bioinformatics and Biostatistics (UOC-UB),
# deposited at https://hdl.handle.net/10609/97486
# Original file: SVMLinealcontinuos.R (dated 2019-06-06 in the deposited annex),
# original encoding ISO-8859-1 (Latin-1), CRLF line endings.
# Unchanged except for text encoding (converted to UTF-8) and this header.
# Original SHA-256: f84ca6c617801a5fbe4f3451d0e9d8dc7348ef13d3fd3a3f2473e40a9e54622c
# ----------------------------------------------------------------------------
#Preprocesamiento de los datos:
##Cargamos las librerias que vamos a usar:
library(e1071)
library(resample)
library(plyr)
library(caTools)
library(ggplot2)
library(miscTools)


##Leemos el archivo.rds que creamos anteriormente y representa la matriz de expresion:
data <- readRDS("finalmatrix.rds")

##Creacion de subsets dependiendo de la Varianza:
##Calculamos las varianzas de cada columna del dataframe y las ordenamos en orden decreciente
##para seleccionar mas tarde los X genes de expresion con varianzas mas altas.
n <- 2000
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
                  cost = 0.1,
                  cross=10,
                  decision.values = T)

##Prediccion con el modelo ajustado de los 
##datos reservados como testSet.
prediction <- predict(svmfit.opt,data_test)
prediction <- as.numeric(prediction)

##Calculo de R^2:
rsq <- function (x, y) cor(x, y) ^ 2
r2 <- rsq(data_test_class, data_test_class - prediction)

##Representacion gráfica de la R^2 obtenida:
p <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=data_test_class, pred=prediction))
p + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("SVM Regression in R r^2=", r2, sep=""))

##Guardado como archivo .RData el modelo ajustado:
save(svmfit.opt, file = "fit_SVMlineal_continuo_8020_erlotinib.RData")


#Partición 60/20/20 (3:1:1 training, validation, testing)#####
##Division de los datos en 80% y 20%:
samp <- sample(1:nrow(gen.var.AUC),.8*nrow(gen.var.AUC))
remain80 <- gen.var.AUC[samp,] ##  80% 
testing20 <- gen.var.AUC[-samp,]  ## 20%
class_testing20 <- testing20[,"class"]

##Division del 80% en train y validation (3:1):
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
svmfit.opt602020 <- svm(class ~.,
                  data = remain80,
                  kernel = "linear",
                  cost = 0.1,
                  cross=10,
                  decision.values = T)

##Prediccion con el modelo ajustado de los 
##datos reservados como testSet.
prediction602020 <- predict(svmfit.opt602020,testing20)

##Calculo de R^2:
r2_602020<- rSquared(class_testing20, class_testing20-prediction602020)

##Representacion gráfica de la R^2 obtenida:
p_602020 <- ggplot(aes(x=actual, y=pred),
            data=data.frame(actual=class_testing20, pred=prediction602020))
p_602020 + geom_point() +
  geom_abline(color="red") +
  ggtitle(paste("SVM Regression in R r^2=", r2_602020, sep=""))

##Guardado como archivo .RData el modelo ajustado:
save(svmfit.opt, file = "fit_SVMlineal_continuo_602020_erlotinib.RData")