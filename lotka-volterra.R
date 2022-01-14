# Russell Legate-Yang
# Simulating Lotka-Volterra

library(ggplot2)

evolve = function(x, y, k, l, a, b, h) {
    "
    Evolve the system by numerically computing the derivative with width h
    "
    xdot = k * x - a * x * y
    ydot = -l * y + b * x * y
    return(c(x + xdot * h, y + ydot * h))
}

sim = function(x_naught, y_naught, k, l, a, b, h, nsims) {
  "
  Given starting value (x_naught, y_naught) and parameters k, l, a, b > 0,
  with derivative width h > 0 , simulate the system for nsims time periods
  "
  xstates = vector(length = nsims + 1, mode = 'numeric')
  xstates[1] = x_naught
  ystates = vector(length = nsims + 1, mode = 'numeric')
  ystates[1] = y_naught

  for (i in 1:nsims) {
    next_states = evolve(xstates[i], ystates[i], k, l, a, b, h)
    xstates[i+1] = next_states[1]
    ystates[i+1] = next_states[2]
  }
  states = data.frame(rabbits = xstates, wolves = ystates, initial = rep(sprintf('(%g,%g)',x_naught, y_naught), (nsims + 1)))
  return(states)
}

gen_plot = function(x_naughts, y_naughts, k, l, a, b, h, nsims, filename) {
  "
  Simulate and generate plot given parameters
  "
  x_eqm = l / b
  y_eqm = k / a
  subtitle = sprintf('Parameters: k=%g, l=%g, a=%g, b=%g', k, l, a, b)
  curves = do.call(rbind, lapply(1:length(x_naughts), function(x) sim(x_naughts[x], y_naughts[x], k, l, a, b, h, nsims)))

  ggplot(data = curves, aes(x = rabbits, y = wolves, color = initial)) +
    geom_path(arrow=arrow(length = unit(0.2, 'cm'))) +
    geom_point(aes(x=x_eqm, y = y_eqm), color = 'black') +
    labs(title = 'Lotka-Volterra phase curves', subtitle = subtitle) 
}

x_naughts = c(0.5, 0.75, 1, 1.25)
y_naughts = c(0.5, 0.75, 1, 1.25)
k = 3
l = 1
a = 4/3
b = 4
h = 0.0001
nsims = 100000
gen_plot(x_naughts, y_naughts, k, l, a, b, h, nsims) 

