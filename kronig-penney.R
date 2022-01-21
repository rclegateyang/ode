# Russell Legate-Yang
# Kronig-Penney model

library(ggplot2)

potential = function(x, H) {
  "
  Compute the potential V(x) of height H
  "
  modtwo = x %% 2
  if (0 < modtwo & modtwo < 1) {
    return(0)
  }
  else {
    return(H)
  }
}

update_state = function(E, x, y, ydot, time_step, H) {
  "
  Compute the second derviative d^2y/dx^2 = (E-V(x))y
  Then update ydot and y 
  "
  ydotdot = (E - potential(x, H)) * y
  ydot = ydot + ydotdot * time_step
  y = y + ydot * time_step
  return(c(y, ydot))
}

plot_ode = function(states, H, E, filename) {
  "
  Plot evolution of y 
  "
  subtitle = sprintf("H = %g, E = %g, (y, y') = (%g, %g)", H, E, states[[2]][1], states[[3]][1])
  p = ggplot(data = data.frame(X = states[[1]], Y = states[[2]]), aes(x = X, y = Y)) +
    geom_path() + 
    labs(title = 'Kronig-Penney model', subtitle = subtitle)
  ggsave(filename, p)
}

compute_ode = function(y, ydot, H, E, max_time, time_step) {
  "
  Given initial conditions (y, ydot) and parameters (H, E), evolve the system
  for max_time in step time_step
  "
  xs = seq(0, max_time, by = time_step)
  numperiods = length(xs)
  
  ys = vector(mode = 'numeric', length = numperiods)
  ys[1] = y 
  ydots = vector(mode = 'numeric', length = numperiods)
  ydots[1] = ydot
  
  for (p in 1:(numperiods-1)) {
    next_state = update_state(E, xs[p], ys[p], ydots[p], time_step, H)
    ys[p+1] = next_state[1]
    ydots[p+1] = next_state[2]
  }
  return(list(xs, ys, ydots))
}

check_stability = function(states, split, bound) {
  "
  Numerically evaluate whether ODE is stable. Split the values Y, taking max |y|
  in each split. Unstable if max in last split is > bound * max in first split
  "
  mod_ys = abs(states[[2]])
  split_point = floor(length(mod_ys) * split)
  first_max = max(mod_ys[1:split_point])
  last_max = max(mod_ys[(split_point+1):length(mod_ys)])
  if (last_max > first_max * bound) {
    return('Unstable')
  }
  else {
    return('Stable')
  }
}

check_stability_IVP = function(H, E, ymin, ymax, ystep, max_time, time_step, split, bound) {
  "
  Given H, E, check stability for intial values in the square
  [ymin, ymax]^2 in steps ystep
  "
  for (y in seq(ymin, ymax, ystep)) {
    for (ydot in seq(ymin, ymax, ystep)) {
      states = compute_ode(y, ydot, H, E, max_time, time_step)
      if (check_stability(states, split, bound)=='Unstable') {
        return('Unstable')
      }
    }
  }
  return('Stable')
}

check_stability_E = function(H, Estep, ymin, ymax, ystep, max_time, time_step, split, bound, filename) {
  "
  For each value of E in -5H to 5H in steps of Estep, check if ODE is stable in
  the square of intial values [ymin, ymax]^2 in steps ystep. Evolve system for
  max_time in time_step. Check stability with split and bound
  (see fn check_stability). Save results to plot in filename
  "
  minE = -5 * H
  maxE = 5 * H
  Es = seq(minE, maxE, by = Estep)
  stability_names = vapply(Es, function(E) check_stability_IVP(H, E, ymin, ymax, ystep, max_time, time_step, split, bound), character(1))
  stability = ifelse(stability_names == 'Stable', 0, 1)
  stable_Es = data.frame(E = Es, Stable = stability, Stability = stability_names)
  p = ggplot(stable_Es, aes(x=E, y=0, color = Stability)) +
    geom_point(size = 3)  +
    annotate("segment", x = minE, xend = maxE, y = 0, yend = 0, size = 1) +
    annotate("segment", x = minE, xend = minE, y = -0.1, yend = 0.1, size = 1) +
    annotate("segment", x = maxE, xend = maxE, y = -0.1, yend = 0.1, size = 1) +
    scale_x_continuous(limits = c(minE, maxE), breaks = stable_Es$E[seq(1, nrow(stable_Es), 2)]) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_color_manual(values = c(Unstable = 'red', Stable = 'blue')) +
    labs(title = 'Stability in E', subtitle = sprintf('H = %g', H)) + 
    theme_bw() + 
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_line(color = "grey",
                                            size = 0.1,
                                            linetype = 2),
          panel.grid.major.x = element_line(color = "grey",
                                            size = 0.1),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank())
  ggsave(filename, p, width = 7, height = 3)
}

# Check stability in E
H = 2
Estep = 0.5
ymin = -2.5
ymax = 2.5
ystep = 0.5
max_time = 25
time_step = 0.001
split = 0.75
bound = 1.25
filename = '/Users/rclegateyang/Documents/Math/ODEs/results/kronig-penney/stabilityE.png'

check_stability_E(H, Estep, ymin, ymax, ystep, max_time, time_step, split, bound, filename)

# Plot some stable and unstable solutions
y = 1
ydot = 1
Es_to_plot = c(-3, -1, 0, 2)
filename = '/Users/rclegateyang/Documents/Math/ODEs/results/kronig-penney/E_'

for (E in Es_to_plot) {
  states = compute_ode(y = y, ydot = ydot, E = E, H = H, max_time = max_time, time_step = time_step)
  plot_ode(states, H, E, paste0(filename, E, '.png'))
}

