library(phyloseq) 
library(INLA)
######################################################## 
load("/Users/Sooran/Desktop/Sooran Kim/!! Research/04 Huilin Li/02 Enviroment-Host-Microbe Interactions/Manuscript/R/sample_data.RData")
######################################################## Data
X <- data_gen$X_li
adj_mat <- data_gen$adj_mat
micro_adj_mat <- data_gen$adj_mat_micro

physeq <- data_gen$physeq
micro_mat <- t(otu_table(physeq)) 
micro_vec <- c(t(micro_mat))
######################################################## Prior
beta_prec1 <- 0.5; beta_prec2 <- 0.001
u_prec1 <- 1; u_prec2 <- 1
u_lambda1 <- 1; u_lambda2 <- 1
v_prec1 <- 0.1; v_prec2 <- 0.1
w_prec1 <- 0.01; w_prec2 <- 0.01
w_lambda1 <- 5; w_lambda2 <- 1
######################################################## 
set.seed(123)
n <- 81; 
indiv <- c(1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 3, 1, 4, 2, 2, 1, 1, 3, 4, 1, 1, 1, 1, 5, 
           1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 3, 1, 2, 1, 1, 1, 13, 2, 3, 7, 1, 
           2, 4, 5, 38, 2, 6, 1, 1, 3, 1, 2, 5, 11, 1, 1, 1, 3, 1, 65, 17, 2, 13, 
           29, 21, 7, 4, 7, 8, 3, 7, 4, 2, 4, 1)
N <- dim(micro_mat)[1]; J <- dim(micro_mat)[2]
Xmat <- kronecker(X*rep(1, N), diag(1, J))
colnames(Xmat) <- paste("X", 1:ncol(Xmat), sep = "")

data_df <- data.frame(micro_vec = micro_vec, 
                      int=factor(rep(1:J)),
                      beta_idx = factor(rep(1:J)),
                      Xmat,
                      Xdf = rep(X, each = J),
                      region=factor(rep(1:n, times = indiv * J)),
                      sampling_rnd=factor(rep(1:N, each = J)),
                      w_rnd = factor(rep(1:J, N)))

######################################################## independent model
formula_ind <- paste("micro_vec ~ 0 +  int +",
                     paste(paste0("X", 1:J), collapse = " + "),
                     "+ f(beta_idx, Xdf, model='iid', hyper=list(prec = list(prior = 'pc.prec', param = c(", beta_prec1, ",", beta_prec2,"))))",
                     "+ f(sampling_rnd, model = 'iid', hyper=list(prec = list(prior = 'loggamma', param = c(", v_prec1, ",", v_prec2, "))))"
)
formula_ind <- as.formula(formula_ind)

result_ind <- inla(formula_ind,
       family = "zeroinflatedpoisson1", data = data_df, num.threads = 10,
       control.fixed = list(prec.intercept = 0.001),
       control.compute = list(config = TRUE))

######################################################## proposed model
regional_values <- as.factor(unique(levels(data_df$region)))
micro_values <- as.factor(unique(levels(data_df$w_rnd)))

formula <- paste("micro_vec ~ 0 + int +",
                 paste(paste0("X", 1:J), collapse = " + "),
                 "+ f(beta_idx, Xdf, model='iid', hyper=list(prec = list(prior = 'pc.prec', param = c(", beta_prec1, ",", beta_prec2,"))))",
                 "+ f(region, values = regional_values, model = 'besagproper2', graph = adj_mat,
                                    hyper = list(prec = list(prior = 'loggamma', param = c(", u_prec1, ",", u_prec2,")),
                                                 lambda = list(prior = 'logitbeta', param=c(", u_lambda1, ",", u_lambda2,"))))",
                 "+ f(sampling_rnd, model = 'iid', hyper=list(prec = list(prior = 'loggamma', param = c(", v_prec1, ",", v_prec2, "))))",
                 "+ f(w_rnd, values = micro_values, model = 'besagproper2', graph=micro_adj_mat,
                                    hyper = list(prec = list(prior = 'loggamma', param = c(", w_prec1, ",", w_prec2, ")),
                                                 lambda = list(prior = 'logitbeta', param=c(", w_lambda1, ",", w_lambda2, "))))"
)

formula <- as.formula(formula)

result <- inla(formula,
       family = "zeroinflatedpoisson1", data = data_df, num.threads = 10,
       control.fixed = list(prec.intercept = 0.001),
       control.compute = list(config = TRUE))

