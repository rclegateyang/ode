# Russell Legate-Yang
# Period doubling bifurcations for logistic map

library(ggplot2)

phi = function(r, x) {
  return(r * x * (1 - x))
}

orbit = function(r, N, x0, Nm = 2) {
  X = vector(length = (Nm * N + 1), mode = 'numeric')
  X[1] = x0
  for (i in 1:(length(X)-1)) {
    X[i+1] = phi(r, X[i])
  }
  Xkeep = unique(tail(X, (N+1)))
  df = data.frame(r = rep(r, length(Xkeep)), x = Xkeep)
  return(df)
}

rorbit = function(N, x0 = 0.5, bounds = c(0, 4)) {
  step = 1/N
  R = seq(bounds[1], bounds[2], step)
  df = do.call(rbind, lapply(R, function(r) orbit(r, N, x0)))
}

N = 10 ** 3
filename = '/Users/rclegateyang/Documents/Math/ODEs/results/logistic/1000_34.png'
plot_df = rorbit(N)
p = ggplot(data = plot_df, aes(x = r, y = x)) +
  xlim(3, 4) +
  geom_point(size = 0.001, alpha = 0.01)
ggsave(filename, plot = p, width = 7, height = 5)