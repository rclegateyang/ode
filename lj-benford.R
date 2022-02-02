# Russell Legate-Yang
# Lissajous curves and Benford's law
library(ggplot2)


# Lissajous curves --------------------------------------------------------
tmax = 100
tstep = 0.001

lj_paths = function(a, b, tmax, tstep) {
  ts = seq(0, tmax, by = tstep)
  xs = cos(a * ts)
  ys = sin(b * ts)
  params = rep(sprintf('a=%g, b=%g', a, b), length(ts))
  df = data.frame(x = xs, y = ys, Parameters = params)
  return(df)
}

lj_plot = function(df, filename) {
  p = ggplot(data = df, aes(x = x, y = y, color = Parameters)) +
    geom_path(arrow = arrow(length = unit(0.2, 'cm'))) + 
    labs(title = 'Lissajous curves')
    ggsave(filename, plot = p, width = 7, height = 5)
}

# a = 1, b integer 
filenamea1bint = '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/a1bint.png'
intub = 3
as = rep(1, intub)
bs = seq(intub)
dfa1bint = do.call(rbind, lapply(1:intub, function(i) lj_paths(as[i], bs[i], tmax = tmax, tstep = tstep)))
lj_plot(dfa1bint, filenamea1bint)

# a integer, b = 1  
filenameaintb1 = '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/aintb1.png'
as = seq(intub) 
bs = rep(1, intub)
dfaintb1 = do.call(rbind, lapply(1:intub, function(i) lj_paths(as[i], bs[i], tmax = tmax, tstep = tstep)))
lj_plot(dfaintb1, filenameaintb1)

# a = 1, b irrational
filenamea1birr = '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/a1birr.png'
as = rep(1, intub)
bs = c(0.5 + sqrt(5) * 0.5, 0.5 - sqrt(5) * 0.5, pi)
dfa1birr = do.call(rbind, lapply(1:intub, function(i) lj_paths(as[i], bs[i], tmax = tmax, tstep = tstep)))
lj_plot(dfa1birr, filenamea1birr)

# lambda * a, lambda * b
filenamescaled = '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/scaled.png'
lambda = 2
as = 1 * (lambda ** seq(0, 1))
bs = 2 * (lambda ** seq(0, 1))
dfscaled = do.call(rbind, lapply(1:2, function(i) lj_paths(as[i], bs[i], tmax = tmax, tstep = tstep)))
lj_plot(dfscaled, filenamescaled)

# Benford's law -----------------------------------------------------------
benford_plot = function(base, N, filename) {
  counts = table(as.numeric(substr(base ** seq(1, N), 1, 1)))
  df = data.frame(Digit = names(counts), Frequency = (as.vector(counts) / N))
  subtitle = sprintf('%g^n, n=1, ..., %g', base, N)
  p = ggplot(data = df, aes(x = Digit, y = Frequency)) +
    geom_bar(stat = "identity") + 
    labs(title = "Benford's law", subtitle = subtitle)
  ggsave(filename, plot = p, width = 7, height = 5)
}
benford_plot(2, 1000, '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/base2N1000.png')
benford_plot(3, 500, '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/base3N500.png')
benford_plot(13, 250, '/Users/rclegateyang/Documents/Math/ODEs/results/lj-benford/base13N250.png')

