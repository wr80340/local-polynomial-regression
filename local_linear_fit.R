library(MASS)
library(locpol)
epan = function(u,h){
  tmp = ifelse(abs(u/h)<1, 0.75*(1-(u/h)^2)/h, 0)
  return(tmp)
}
cosine = function(u,h){
  tmp = ifelse(abs(u/h)<1, (pi/4)*cos(pi*(u/h)/2)/h, 0)
  return(tmp) 
}
biweight = function(u,h){
  tmp = ifelse(abs(u/h)<1, (15/16)*((1-(u/h)^2)^2)/h, 0)
  return(tmp)
}
triweight = function(u,h){
  tmp = ifelse(abs(u/h)<1, (35/32)*((1-(u/h)^2)^3)/h, 0)
  return(tmp)
}
gauss = function(u,h){
  tmp = (1/sqrt(2*pi))*exp(-(0.5*(u/h)^2))/h
  return(tmp)
}
# summary.lpr <- function(lpr_object){
# }

plot.lpr = function(lpr_obj){
  par(mfrow = c(2,1), mai = c(0.8,0.8,0.5,0.1))
  ## basis plot
  plot(lpr_obj$X, lpr_obj$y, main = "kernel fit of")
  for(i in 1:length(lpr_obj$knots)){
    lines(x = lpr_obj$plot_element_1[[i]][1,], y = lpr_obj$plot_element_1[[i]][2,], col = "red", lwd = 1)
  }
  ## tot_plot
  tmp = t(lpr_obj$plot_element_1[[1]])
  for(i in 2:length(lpr_obj$plot_element_1)){
    tmp = rbind(tmp, t(lpr_obj$plot_element_1[[i]]))
  }
  d = matrix(rep(lpr_obj$plot_element_2,2), nrow = length(lpr_obj$plot_element_2), ncol = 2)
  for(i in 1:length(lpr_obj$plot_element_2)){
    d[i,2] = mean(tmp[tmp[,1] == lpr_obj$plot_element_2[i],2])
  }
  # d is a dim = (length(fitted_line), 401) matrix
  plot(lpr_obj$X, lpr_obj$y, main = "fitted line")
  lines(d[,1], d[,2], lwd = 2, col = "red")
  points(lpr_obj$X, lpr_obj$fitted_value, pch = 24, col = "black", bg = "orange")
  par(mfrow = c(1,1))
}
predict.lpr = function(lpr_obj, new_data){
  predict_val = sapply(new_data, function(x){
    ind = x-lpr_obj$h < lpr_obj$knots & lpr_obj$knots < x+lpr_obj$h
    tmp_knots = lpr_obj$knots[ind]
    # dim = (deg+1 * -)
    tmp_beta = ifelse(lpr_obj$degree == 0,
                      matrix(lpr_obj$coef[ind], nrow = 1),
                      lpr_obj$coef[, ind])
    # dim = (deg+1 * -)
    tmp_X = sapply(tmp_knots, function(t){
      return(c(1, (x-t), (x-t)^2)[1:(lpr_obj$degree+1)])
    })
    # same shape as tmp_beta
    res = ifelse(lpr_obj$degree == 0,
                 mean(tmp_beta*tmp_X),
                 mean(apply(tmp_beta*tmp_X, 2, sum)))
    return(res)
  })
  return(list(X = new_data,
              y = predict_val))
}
X = iris$Sepal.Length; y = iris$Petal.Length; h = 0.5; deg = 1; ker = "epan"; nok = 51
X_test = iris$Sepal.Length[101:150]; y_test = iris$Petal.Length[101:150]
lpr = function(X, y, h = 0.5, deg = 1, ker = "epan", nok = 51, space = "equal"){
  if(ker == "epan"){
    Ker = epan
  }else if(ker == "cosine"){
    Ker = cosine
  }else if(ker == "biweight"){
    Ker = biweight
  }else if(ker == "triweight"){
    Ker = triweight
  }else if(ker == "gaussian"){
    Ker = gauss
  }else{
    Ker = epan
    print("kernel must use (1.) cosine (2.) epan, (3.) biweight, (4.) triweight, (5.) gaussian, use epan instead")
  }
  stopifnot(length(X) == length(y), h >= 0, deg %in% c(0,1,2), nok >= 0)
  N = length(X)
  r = range(X)
  knot_line = seq(from = r[1], to = r[2], length.out = nok)
  beta = sapply(knot_line, function(t){
    # xi - x0
    xi = X - t
    # under certain t 
    x = matrix(c(rep(1, N), xi, xi^2), N, 3)
    # (N, 3) matrix
    if(deg == 0){
      x = x[,1]
    }else if(deg == 1){
      x = x[,1:2]
    }
    W = diag(Ker(u = xi, h = h))
    # (N, N) matrix
    beta = ginv(t(x)%*%W%*%x)%*%t(x)%*%W%*%y
    # length = deg vector
    return(beta)
  })
  beta = matrix(beta, nrow = deg+1, ncol = nok)
  fitted_val = sapply(X, function(x){
    tmp_knots = knot_line[x-h < knot_line & knot_line < x+h]
    if(deg == 0){
      tmp_beta = beta[x-h < knot_line & knot_line < x+h]
      tmp_X = rep(1, length(tmp_knots))
      return(mean(tmp_beta*tmp_X))
      }else if(deg == 1){
        tmp_beta = beta[, x-h < knot_line & knot_line < x+h]
        tmp_X = sapply(tmp_knots, function(t){
          return(c(1, (x-t)))
          })
        return(mean(apply(tmp_beta*tmp_X, 2, sum)))
        }else if(deg == 2){
          tmp_beta = beta[, x-h < knot_line & knot_line < x+h]
          tmp_X = sapply(tmp_knots, function(t){
            return(c(1, (x-t), (x-t)^2))
            })
          return(mean(apply(tmp_beta*tmp_X, 2, sum)))
          }
    # same shape as tmp_beta
    })
  plot_line = seq(from = r[1], to = r[2], length = 401)
  fitted_line = lapply(1:length(knot_line), function(t){
    basis_line = plot_line - knot_line[t]
    ind = which(abs(basis_line) < h)
    basis_line = basis_line[ind]
    plot_line = plot_line[ind]
    # basis_line & plot_line have same length
    # less points
    basis_line = rbind(1,basis_line, basis_line^2)
    # (3, nok) matrix
    if(deg == 0){
      basis_line = basis_line[1,]
    }else if(deg == 1){
      basis_line = basis_line[1:2,]
    }
    # basis_line changes to a (deq+1, nok) matrix
    if(length(beta[,t]) == 1){
      y_pred = matrix(beta[,t]*basis_line, nrow = 1)
    }else{
      y_pred = matrix(t(beta[,t])%*%basis_line, nrow = 1)
    }
    # y_pred ia a (1, nok) matrix
    fit_line = rbind(plot_line, y_pred)
    return(fit_line)
  })
  SSE = sum((y - fitted_val)^2)
  MSE = mean((y - fitted_val)^2)
  call = paste0("lpr(formula = ", deparse(substitute(y)),
                " ~ ", deparse(substitute(X)), 
                ", h = ", h,
                ", degree = ", deg, 
                ", number of knots = ", nok, ")")
  res = list(call = call,
             X = X,
             # from data
             y = y,
             # from data
             degree = deg,
             # from data
             h = h,
             range = r,
             # from data
             SSE = SSE,
             MSE = MSE,
             coef = beta,
             knots = knot_line,
             # knots location
             fitted_value = fitted_val,
             # with length = length of X
             plot_element_1 = fitted_line,
             plot_element_2 = plot_line)
  class(res) = "lpr"
  return(res)
}
methods(class = "lpr")
lpr.fit = lpr(X, y, h = 0.35, deg = 2, nok = 401)
plot(lpr.fit)
pred = predict(lpr.fit, X_test)

plot(X,y)
points(X_test, y_test, col = "blue", pch = 19)
points(pred$X, pred$y, col = "green", pch = 19)

lpr.fit = lpr(X,y, h = 0.7, deg = 2)
plot(lpr.fit)

lpr.fit = lpr(X,y, h = 0.7, deg = 1)
lpr.fit = lpr(X,y, h = 0.7, deg = 0)

s_other = Sys.time()
locpol(Petal.Length~Sepal.Length, data = iris, bw = 0.5)
time_other = Sys.time()-s_other

s = Sys.time()
lpr(X,y, h = 0.5, deg = 1, nok = 100)
time_us = Sys.time()-s