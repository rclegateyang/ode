# Russell Legate-Yang
# van der Pol model

# Let x = x and xdot = x', so the system is:
# x' = xdot, xdot' = mu(1-x)^2 * xdot + x

library(ggplot2)

update_state = function(mu, x, xdot, time_step) {
  "
  Evolve the system in time increment time_step according to 
  x' = xdot, xdot' = mu(1-x)^2 * xdot + x
  "
  new_x = x + xdot * time_step
  new_xdot = xdot + (mu * (1 - x ** 2) * xdot - x) * time_step
  return(c(new_x, new_xdot))
}

solve_ode = function(mu, time_step, max_time, initial_x, initial_xdot) {
  "
  Evolve ODE for max_time in steps of time_step given intial data and mu
  "
  ts = seq(0, max_time, by = time_step)
  numperiods = length(ts)
  xs = vector(mode = 'numeric', length = numperiods)
  xs[1] = initial_x
  xdots = vector(mode = 'numeric', length = numperiods)
  xdots[1] = initial_xdot
  
  for (p in 1:(numperiods-1)) {
    next_state = update_state(mu, xs[p], xdots[p], time_step)
    xs[p+1] = next_state[1]
    xdots[p+1] = next_state[2]
  }
  return(list(xs, xdots))
}

plot_odes = function(states_by_group, group_name, filename) {
  "
  Given a list states_by_group with names group and entries return values of 
  solve_ode, plot them all on the same (x, xdot) plane.
  The grouping variable can be mu or initial conditions. group_name identifies
  the grouping variable
  "
  numperiods = length(states_by_group[[1]][[1]])
  df = data.frame(do.call(rbind, lapply(states_by_group, function(x) matrix(unlist(x), nrow = numperiods))))
  df$temp = rep(names(states_by_group), each = numperiods)
  names(df) = c('x', 'xdot', group_name)
  
  subtitle = paste('Varying', group_name)
  p = ggplot(data = df, aes_string(x = 'x', y = 'xdot', color = group_name)) +
    geom_path(arrow=arrow(length = unit(0.2, 'cm'))) + 
    labs(title = 'van der Pol model', subtitle = subtitle)
  ggsave(filename, p, width = 7, height = 5)
}

vary_mu = function(initial_x, initial_xdot, mus, filename, max_time) {
  states_by_mu = lapply(mus, solve_ode, time_step = time_step, max_time = max_time, initial_x = initial_x, initial_xdot = initial_xdot)
  names(states_by_mu) = mus
  plot_odes(states_by_mu, 'Mu', filename)
}
# Vary mu with initial condition fixed
time_step = 0.0001
max_time = 50
initial_x = 1
initial_xdot = 1

filenamelow = '/Users/rclegateyang/Documents/Math/ODEs/results/van-der-pols/11mulow.png'
vary_mu(initial_x, initial_xdot, c(0.01, 0.1, 0.5), filenamelow, max_time)
filenamehigh = '/Users/rclegateyang/Documents/Math/ODEs/results/van-der-pols/11muhigh.png'
vary_mu(initial_x, initial_xdot, c(1, 2, 5), filenamehigh, max_time = 100)
filenamehigh = '/Users/rclegateyang/Documents/Math/ODEs/results/van-der-pols/11mulowhigh.png'
vary_mu(initial_x, initial_xdot, c(0.01, 10), filenamehigh, max_time = 100)

