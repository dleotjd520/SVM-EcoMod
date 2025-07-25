
# install.packages("plot3D")
# install.packages("e1071")
# install.packages("concaveman")
# install.packages("caret")
# install.packages("sf")
##############################################################
########################################  1.1. Kernel trick
library(plot3D)
# Data generation with specific pattern
set.seed(123)
n <- 300
theta <- runif(n, 0, 2*pi)
r <- c(runif(n/2, 0, 0.5), runif(n/2, 0.8, 1.2))
x1 <- r * cos(theta)
x2 <- r * sin(theta)
label <- (c(rep(1, n/2), rep(2, n/2)))
data <- data.frame(x1 = x1, x2 = x2, label = label)

### Data visualization (2D: original space)
par(mar=c(4,4,3,3))
plot(data$x1, data$x2, bg = c("blue", "tomato")[data$label], 
     las=1, pch = 21, cex=1.2, xlab = "x1", ylab = "x2", main = "2D")

# Data mapping Using RBF Kernel
gamma <- 1
data$z <- exp(-gamma * (x1^2 + x2^2)) *2

### Data visualization (3D; feature space)
par(mar=c(4,4,3,3))
scatter3D(data$x1, data$x2, data$z, bg  = c("blue", "tomato")[data$label], col="black",
          pch = 21,
          xlab = "x1", ylab = "x2", zlab = "z",
          phi = 15, theta = 45,
          colkey = F,
          # colkey = list(side = 1, length = 0.5, width = 0.5),
          main = "3D Transformation using RBF Kernel",
          bty = "g", ticktype = "detailed", d = 3)

# Generate grid for x and y sequences
grid <- expand.grid(x_seq = seq(min(data$x1), max(data$x1), length.out = 10),
                    y_seq = seq(min(data$x2), max(data$x2), length.out = 10))

# Calculate z values for the hyperplane (arbitrary hyperplane)
z_plane <- with(grid, 0.05 * x_seq - 0.05 * y_seq + 1.4)

# Reshape the data into matrices for surf3D
grid_x <- matrix(grid$x_seq, nrow = 10)
grid_y <- matrix(grid$y_seq, nrow = 10)
z_matrix <- matrix(z_plane, nrow = 10)

# Visualize of hyperplane
surf3D(x = grid_x, y = grid_y, z = z_matrix, 
       colvar = NULL, col=NA, lwd=1.5, border = "gray30", add = TRUE)

##############################################################
########################################  1.2. SVC
########################################  (1) Linear SVC
# Random data generation
set.seed(1512)
data <- data.frame(x1 = rnorm(20), x2 = rnorm(20), y = rep(c(1, 2), c(10, 10)))
data[c("x1", "x2")][data$y == 1, ] <- data[c("x1", "x2")][data$y == 1, ] + 1.5

# Check data
par(mar=c(4,4,3,3))
plot(data[1:2], col=c("blue","tomato")[data$y], pch=19, las=1, cex=1.5)

# Develope SVM model
library(e1071)
svmfit = svm(y ~ ., data = data, kernel = "linear", type="C-classification", scale = FALSE)

# Check support vectors
points(data[1:2][svmfit$index, ], cex=2)

# Calculate decision boundary (linear)
beta = as.numeric(t(svmfit$coefs)%*%as.matrix((data[-3])[svmfit$index,]))
beta0 = svmfit$rho

# Background
x1 = seq(from = range(data$x1)[1], to = range(data$x1)[2], length = 30)
x2 = seq(from = range(data$x2)[1], to = range(data$x2)[2], length = 30)
grid <- expand.grid(x1 = x1, x2 = x2)

# prediction
pred.grid = predict(svmfit, grid)

### Final figure (data points + support vectors + decision boundary + class area)
par(mar=c(4,4,3,3))
plot(grid, col = c("blue","tomato")[as.numeric(pred.grid)], pch = "+", cex = 0.4, las=1)
abline(beta0 / beta[2], -beta[1] / beta[2], lwd=2) # optimal hyperplane
abline((beta0 - 1) / beta[2], -beta[1] / beta[2], lty = 5, lwd=1.5) 
abline((beta0 + 1) / beta[2], -beta[1] / beta[2], lty = 5, lwd=1.5)
points(data[1:2][svmfit$index, ], pch=19, col="#a9def9", cex=3) # support vectors
points(data[1:2], col=c("blue","tomato")[data$y], pch=19, cex=1.5, lwd=2) # data points


########################################  (2) Non-linear SVC
# Random data generation
# Develope SVM model
library(e1071)
svmfit = svm(y ~ ., data = data, kernel = "radial", type="C-classification", scale = FALSE)

# Check support vectors
par(mar=c(4,4,3,3))
plot(data[1:2], col=c("blue","tomato")[data$y], pch=19, las=1)

# prediction
pred.grid = predict(svmfit, grid)

### Final figure (data points + support vectors + decision boundary + class area)
par(mar=c(4,4,3,3))
plot(grid, col = c("blue","tomato")[as.numeric(pred.grid)], pch = "+", cex = 0.4, las=1)
contour(x = x1, y = x2, 
        z = matrix(attributes(predict(svmfit, grid, decision.values = TRUE))$decision.values, 
                   nrow = length(x1), ncol = length(x2)), 
        nlevels = 0, lwd=2, drawlabels =F, add=T)
points(data[1:2][svmfit$index, ], pch=19, col="#a9def9", cex=3) # support vectors
points(data[1:2], col=c("blue","tomato")[data$y], pch=19, cex=1.5, lwd=2) # data points


########################################  (3) Non-linear SVC : gamma - C relationship
# Random data generation
set.seed(450)
data <- data.frame(x1 = rnorm(100), x2 = rnorm(100), y = rep(c(1, 2), c(50, 50)))
data[c("x1", "x2")][data$y == 1, ] <- data[c("x1", "x2")][data$y == 1, ] + 2
plot(data[1:2], col=c("blue","tomato")[data$y], pch=19, las=1)

# Background
x1 = seq(from = range(data$x1)[1]*1.5, to = range(data$x1)[2]*1.5, length = 500)
x2 = seq(from = range(data$x2)[1]*1.5, to = range(data$x2)[2]*1.5, length = 500)
grid <- expand.grid(x1 = x1, x2 = x2)

# Model development & prediction
r.gamma = 2^(-1:1)
r.cost = 5^(-1:1)
para.grids <- expand.grid(g=r.gamma, c=r.cost)

for( k in 1:nrow(para.grids)){ # k=3
svmfit = svm(y ~ ., data = data, kernel = "radial", type="C-classification", scale = T, 
             gamma=para.grids$g[k], cost=para.grids$c[k])
pred.grid = predict(svmfit, grid)

# evaluation
library(caret)
ACC <- confusionMatrix( table(data$y, predict(svmfit, data)) )$overall["Accuracy"]

### Final figure (data points + decision boundary + class area)
par(mar=c(4,4,3,3))
plot(data[1:2], type="n", las=1, main= paste0("C=", para.grids$c[k],", gamma=",para.grids$g[k], ", ACC=", ACC))
library(concaveman)
grid.1 <- grid[as.numeric(pred.grid)==1,]
grid.1b <- concaveman(as.matrix(grid.1), concavity=0.1) # boundary

library(sf)
poly.bg <- st_sfc(st_polygon(list(concaveman(as.matrix(grid)))))
poly.1 <- st_sfc(st_polygon(list(grid.1b)))
poly.bg <- st_difference(poly.bg, poly.1)

plot(poly.bg, col=adjustcolor("tomato", alpha.f = 0.2), add=T, lwd=2)
plot(poly.1, col=adjustcolor("blue", alpha.f = 0.2), add=T, lwd=2)
points(data[1:2], bg=c("blue","tomato")[data$y], pch=21, cex=1.5, lwd=2)
}# for k

##############################################################
########################################  1.3. SVR
# Random data generation
set.seed(123)
x <- seq(0, 10, length = 100)
y <- sin(x) + rnorm(100, sd = 0.2)
data <- data.frame(x, y)

# Check data
par(mar=c(4,4,3,3))
plot(data[c("x", "y")], cex=1, pch=21, bg="lightgray")

# Develope SVM model
library(e1071)
svmfit = svm(y ~ x, data = data, type = "eps-regression", kernel = "radial")

# Check support vectors
points(data[1:2][svmfit$index, ], cex=2)

# prediction
data$pred <- predict(svmfit, data)

# e-tube (margin)
data$upper <- data$pred + svmfit$epsilon
data$lower <- data$pred - svmfit$epsilon

### Final figure 
par(mar=c(4,4,3,3))
plot(data[c("x", "y")], pch=19, las=1, cex=0.8, las=1)
lines(data$x, data$pred, col = "tomato", lwd = 2) # regresion line
lines(data$x, data$upper, col = "darkred", lty = 2, lwd=1.5)    # e-tube upper boundary
lines(data$x, data$lower, col = "darkred", lty = 2, lwd=1.5)    # e-tube lower boundary
points(data[1:2][svmfit$index, ], cex=1.6, col="#a9def9", pch=19)
points(data[c("x", "y")], cex=1, pch=21, bg="lightgray")
