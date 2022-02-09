# Russell Legate-Yang
# Lorenz system
library(plotly)
library(ggplot2)

next_state = function(x, y, z, sigma, rho, beta, dt) {
  next_x = x + (sigma * (y - x)) * dt
  next_y = y + (x * (rho - z) - y) * dt
  next_z = z + (x * y - beta * z) * dt
  return(c(next_x, next_y, next_z))
}

lorenz = function(x0, y0, z0, sigma, rho, beta, dt, maxt) {
  ts = seq(from = 0, to = maxt, by = dt)
  nperiods = length(ts)
  xs = vector(length = nperiods, mode = 'numeric')
  ys = vector(length = nperiods, mode = 'numeric')
  zs = vector(length = nperiods, mode = 'numeric')
  xs[1] = x0
  ys[1] = y0
  zs[1] = z0
  
  for(t in 1:(nperiods-1)) {
    next_sys = next_state(xs[t], ys[t], zs[t], sigma, rho, beta, dt)
    xs[t+1] = next_sys[1]
    ys[t+1] = next_sys[2]
    zs[t+1] = next_sys[3]
  }
  initial_text = sprintf('(%g, %g, %g)', x0, y0, z0)
  df = data.frame(x = xs, y = ys, z = zs, t = ts, initial = rep(initial_text, nperiods))
  return(df)
}

sigma = 10
beta = 8/3
rho = 28
dt = 0.001
maxt = 50

x0s = c(1,1.0001)
y0s = c(1,1.0001)
z0s = c(1,1.0001)
df = do.call(rbind, lapply(1:length(x0s), function(i) lorenz(x0s[i], y0s[i], z0s[i], sigma, rho, beta, dt, maxt)))

plot_3d = function(df) {
  plot_ly(df, x = ~x, y = ~y, z = ~z) %>%
    add_paths(color = ~initial, colors = c('black', 'yellow'))
}
plot_3d(df)

plot_x = function(df, filename) {
  p = ggplot(data = df, aes(x = t, y = x, color = initial)) + 
    geom_line() + 
    labs(title = "Lorenz: x(t)") 
  ggsave(filename, p)
}

ext = '100x.png'
plot_x(df, paste0('/Users/rclegateyang/Documents/Math/ODEs/results/lorenz/', ext))
