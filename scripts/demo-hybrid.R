devtools::load_all("./")
n <- 200

x <- sort(runif(n))

Fhat <- ecdf(x)(x)


bands <- hybrid_band(Fhat, alpha = 0.05, m = n)