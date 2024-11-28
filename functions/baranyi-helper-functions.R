#ATB
#growth rate analysis

library(zoo)

fit_baranyi = function(x, y, tries = 100){
  m1 <- NA
  #figure out start lag guess
  lag_guess = guess_half_max(x,y)
  while(is.na(m1[1]) && (tries > 0)){
    m1 <- tryCatch( nls(y ~ baranyi(x, r, lag, ymax, y0),
                        start = list(r = runif(1, 0.1, 0.4), 
                                     lag = lag_guess + sample(-20:20, 1), 
                                     ymax = max(y), 
                                     y0 = min(y))),
                    error = function(e) return(NA))
    tries = tries - 1
  }
  if (is.na(m1)){
    return(c(NA, NA,NA,NA))
  }
  #   m1 <- nls(y ~ baranyi(x, r, lag, ymax, y0),
  #             start = list(r = 0.1, lag = x[round(length(x) / 3)], ymax = max(y), y0 = min(y)))
  r = coef(m1)[1]
  lag = coef(m1)[2]
  ymax = coef(m1)[3]
  y0 = coef(m1)[4]
  return(c(r, lag, ymax, y0))
}

guess_half_max = function(x, y){
  # this function looks at a logistically-increasing time series
  # and guesses when the growth rate is at a maximum (which is also when
  # the trend is at half-max)
  
  # put in order
  y = y[sort(x, index.return = TRUE)$ix]
  x = sort(x)
  
  #find approximate time of max growth rate, using diffs, on smoothed y
  y2 = diff(y)
  y2 = zoo::rollmean(y2, 11, fill = NA, align = "center")
  half_max_idx = which(y2 == max(y2, na.rm = TRUE))[1]
  half_max = x[half_max_idx]
  return(half_max)
}

baranyi <- function(t, r, lag, logymax, logy0){
  At = t + (1 / r) * log(exp(-r * t) + exp(-r * lag) - exp(-r * (t + lag)))
  logy = logy0 + r * At - log(1 + (exp(r * At) - 1) / exp(logymax - logy0))
  return(logy)
}