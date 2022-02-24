# Russell Legate-Yang
# Methods of integration
library(kableExtra)

# Integration methods -----------------------------------------------------

euler = function(t, x, f, delta) {
  return(x + delta * f(t, x))
}

midpoint = function(t, x, f, delta) {
  ktemp = x + delta * f(t, x) / 2
  return(x + delta * f(t + delta / 2, ktemp))
}

runge_kutta = function(t, x, f, delta) {
  k1 = f(t, x)
  k2 = f(t + delta / 2, x + delta * k1 / 2)
  k3 = f(t + delta / 2, x + delta * k2 / 2)
  k4 = f(t + delta, x + delta * k3)
  return(x + delta * (k1 + 2 * k2 + 2 * k3 + k4) / 6)
}


# Evaluate method ---------------------------------------------------------
evaluate_method = function(f, method, xt, x0, stop_t) {
  "
  Given x' = f with initial condition x(0) = x0 and known solution x(t) = xt,
  evaluate method by integrating to approximate x(stop) and comparing to known
  value. Search N = 2^k for k in [10, 20] and delta = 1/N
  "
  true_val = xt(stop_t, x0)
  
  Ns = 2 ** seq(10, 20)
  errors_abs = vector(length(Ns), mode = 'numeric')
  errors_rel = vector(length(Ns), mode = 'numeric')
  
  for (i in 1:length(Ns)) {
    delta = 1 / Ns[i]
    prev_x = x0
    ts = seq(delta, stop_t, by = delta)
    for (t in ts) {
      prev_x = method(t, prev_x, f, delta)
    }
    errors_abs[i] = true_val - prev_x
    errors_rel[i] = 100 * errors_abs[i] / true_val 
  }
  
  df = data.frame(N = Ns, Error = errors_abs, Size = errors_rel)
  colnames(df) = c("N", 'Absolute error', '% error')
  return(df)
}

# Execution -------------------------------------------------------------------
prop_growth = function(t, x) {
  return(x)
}
soln = function(t, const) {
  return(const * exp(t))
}
x0 = 1
stop_t = 1

eval_euler = evaluate_method(prop_growth, euler, soln, x0, stop_t)
eval_midpoint = evaluate_method(prop_growth, midpoint, soln, x0, stop_t)
eval_runge_kutta = evaluate_method(prop_growth, runge_kutta, soln, x0, stop_t)
eval_df = do.call(cbind, list(eval_euler, eval_midpoint[c(2,3)], eval_runge_kutta[c(2,3)]))
eval_df[, !names(eval_df) == 'N'] = lapply(lapply(eval_df[, !names(eval_df) == 'N'], format, scientific = T, digits = 4), as.character)

output_file = '/Users/rclegateyang/Documents/Math/ODEs/results/numerical-integration/comparison.png'
kbl(eval_df, caption = "First-order numerical integration comparison", booktabs = T, digits = 10) %>%
  kable_styling() %>%
  add_footnote("Error = true - estimated") %>%
  add_header_above(c(' ' = 1, "Euler" = 2, "Midpoint" = 2, 'Runge-Kutta' = 2)) %>%
  save_kable(output_file, zoom = 3)
