## ----setup, include=FALSE-----------------------------------------------------
library(knitr)
knitr::opts_chunk$set(
	fig.align = "center",
	fig.height = 5.5,
	fig.width = 6,
	warning = FALSE,
	collapse = TRUE,
	dev.args = list(pointsize = 10),
	out.width = "65%"
)

## ---- message = FALSE, echo = FALSE-------------------------------------------
library(smmR)
library(DiscreteWeibull)

## ----echo = FALSE, out.width = '400px', fig.cap = 'Waste treatment for a textile factory'----
knitr::include_graphics("Waste_treatment.png")

## ----echo = FALSE, out.width = '400px', fig.cap = 'Three-state discrete-time semi-Markov system modelization'----
knitr::include_graphics("Three_state_modelization.png")

## -----------------------------------------------------------------------------
states <- c("1", "2", "3") # State space

alpha <- c(1, 0, 0) # Initial distribution

p <- matrix(data = c(0, 1, 0, 
                     0.95, 0, 0.05, 
                     1, 0, 0), nrow = 3, byrow = TRUE) # Transition matrix

distr <- matrix(c(NA, "geom", NA, 
                  "dweibull", NA, "dweibull", 
                  "dweibull", NA, NA), 
                nrow = 3, ncol = 3, byrow = TRUE) # Distribution matrix

param1 <- matrix(c(NA, 0.8, NA, 
                   0.3, NA, 0.5,
                   0.6, NA, NA), 
                 nrow = 3, ncol = 3, byrow = TRUE)

param2 <- matrix(c(NA, NA, NA, 
                   0.5, NA, 0.7,
                   0.9, NA, NA), 
                 nrow = 3, ncol = 3, byrow = TRUE)

parameters <- array(c(param1, param2), c(3, 3, 2))

factory <- smmparametric(states = states, init = alpha, ptrans = p, 
                         type.sojourn = "fij", distr = distr, param = parameters)

## -----------------------------------------------------------------------------
M <- 10000
seq <- simulate(object = factory, nsim = M)

## -----------------------------------------------------------------------------
estimate <- fitsmm(sequences = seq, states = states, type.sojourn = "fij")

## -----------------------------------------------------------------------------
print(x = estimate$ptrans, digits = 2)

## -----------------------------------------------------------------------------
plot(x = estimate, i = "2", j = "3", type = "l", col = "blue")

lines(x = 1:estimate$kmax, y = ddweibull(x = 1:estimate$kmax, q = 0.5, beta = 0.7), 
       col = "red", pch = "x")

legend(x = "topright", 
       legend = c("True value", "Estimate"), 
       col = c("red", "blue"), lty = c(1, 1))

## -----------------------------------------------------------------------------
k <- 300
upstates <- c("1", "2") # Working states of the semi-Markov system

## -----------------------------------------------------------------------------
trueReliab <- reliability(x = factory, k = k, upstates = upstates)
estReliab <- reliability(x = estimate, k = k, upstates = upstates)

## -----------------------------------------------------------------------------
plot(x = 0:k, y = trueReliab[, 1], type = "l", cex = 2.5, ylim = c(0, 1), 
     col = "red", main = "Reliability", xlab = "k", ylab = "R(k)")

lines(x = estReliab[, 1], col = "blue")
lines(x = estReliab[, 3], lty = 4, col = "blue")
lines(x = estReliab[, 4], lty = 4, col = "blue")
legend(x = "topright", 
       legend = c("True value", "Estimated value", "95% confidence interval"), 
       col = c("red", "blue", "blue"), lty = c(1, 1, 4))

## -----------------------------------------------------------------------------
trueAvail <- availability(x = factory, k = k, upstates = upstates)
estAvail <- availability(x = estimate, k = k, upstates = upstates)

## -----------------------------------------------------------------------------
plot(x = 0:k, y = trueAvail[, 1], type = "l", cex = 2.5, ylim = c(0.95, 1), 
     col = "red", main = "Availability", xlab = "k", ylab = "A(k)")

lines(x = estAvail[, 1], col = "blue")
lines(x = estAvail[, 3], lty = 4, col = "blue")
lines(x = estAvail[, 4], lty = 4, col = "blue")
legend(x = "topright", 
       legend = c("True value", "Estimated value", "95% confidence interval"), 
       col = c("red", "blue", "blue"), lty = c(1, 1, 4))

## -----------------------------------------------------------------------------
trueBMP <- failureRate(x = factory, k = k, upstates = upstates)
estBMP <- failureRate(x = estimate, k = k, upstates = upstates)

## -----------------------------------------------------------------------------
plot(x = 0:k, y =  trueBMP[, 1], type = "l", cex = 2.5, ylim = c(0, 0.025), 
     col = "red", main = "BMP-failure rate", xlab = "k", ylab = bquote(lambda(k)))

lines(x = estBMP[, 1], col = "blue")
lines(x = estBMP[, 3], lty = 4, col = "blue")
lines(x = estBMP[, 4], lty = 4, col = "blue")
legend(x = "topright", 
       legend = c("True value", "Estimated value", "95% confidence interval"), 
       col = c("red", "blue", "blue"), lty = c(1, 1, 4))

## -----------------------------------------------------------------------------
trueRG <- failureRate(x = factory, k = k, upstates = upstates, failure.rate = "RG")
estRG <- failureRate(x = estimate, k = k, upstates = upstates, failure.rate = "RG")

## -----------------------------------------------------------------------------
plot(x = 0:k, y =  trueRG[, 1], type = "l", cex = 2.5, ylim = c(0, 0.03), 
     col = "red", main = "RG-failure rate", xlab = "k", ylab = "r(k)")

lines(x = estRG[, 1], col = "blue")
lines(x = estRG[, 3], lty = 4, col = "blue")
lines(x = estRG[, 4], lty = 4, col = "blue")
legend(x = "topright", 
       legend = c("True value", "Estimated value", "95% confidence interval"), 
       col = c("red", "blue", "blue"), lty = c(1, 1, 4))

## -----------------------------------------------------------------------------
trueMTTF <- mttf(x = factory, upstates = upstates)
estMTTF <- mttf(x = estimate, upstates = upstates)

## -----------------------------------------------------------------------------
print(trueMTTF)
print(estMTTF)

## -----------------------------------------------------------------------------
trueMTTR <- mttr(x = factory, upstates = upstates)
estMTTR <- mttr(x = estimate, upstates = upstates)

