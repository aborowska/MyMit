
GelmanMeng <- function(x, A = 0, B = sqrt(0.5), C1 = 1, C2 = 2, log = TRUE)
{
  if (is.vector(x))
    x <- matrix(x, nrow = 1)
  r <- -.5 * (A * x[,1]^2 * x[,2]^2 + x[,1]^2 + x[,2]^2
              - 2 * B * x[,1] * x[,2] - 2 * C1 * x[,1] - 2 * C2 * x[,2])
  if (!log)
    r <- exp(r)
  as.vector(r)
}

PlotGelmanMeng <- function(x1, x2)
{
  GelmanMeng(cbind(x1, x2), log = FALSE)
}
x1 <- x2 <- seq(from = -1.0, to = 10.0, by = 0.02)
z <- outer(x1, x2, FUN = PlotGelmanMeng)
image(x1, x2, z, las = 1, col = gray((20:0)/20),
      cex.axis = 1.1, cex.lab = 1.2,
      xlab = expression(X[1]), ylab = expression(X[2]))
box()
abline(a = 0, b = 1, lty = "dotted")


## Run AdMit (with default values)
set.seed(1234)
outAdMit <- AdMit(GelmanMeng, mu0 = c(0.0, 0.1), control = list(Ns = 1e4))
print(outAdMit)