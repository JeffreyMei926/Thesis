


# Install Libraries 
# =================
library("glmnet")
library("selectiveInference")
library("lars")
library("MASS")
library("leaps")
library("bestglm")
library("parcor")

get.cor = function( p, cor = 0.5 ){
  # ---------------------------------------------------------------------------
  # DESCRIPTION: Calculates the correlation matrix in Tibshirani (1996) simulations. 
  #   The correlation matrix is of the form $$ \rho^{ |i-j| } $$
  #  
  # ARGUMENTS: 
  #   p   : number of predictors  
  #   cor : rho  
  # 
  # VALUES: 
  #   cor.matrix : constructed matrix 
  # ---------------------------------------------------------------------------
  
  # Initialize Variables
  cor.matrix = matrix( nrow = p, ncol = p ) 
  
  # Fill in Matrix  
  for( i in 1:p ){
    for( j in 1:p ){
      cor.matrix[ i, j ] = cor^( abs( i - j ) )     
    }
  }
  return( cor.matrix ) 
}

pred = function( b.hat, x ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: obtains the y value such that 
  #         $$ y = x \beta $$
  #   NOTE: it does not include an intercept 
  # 
  # ARGUMENTS:
  #   b.hat : coefficients for model
  #   x     : x values for model  
  #  
  # VALUES:
  #   y.hat: predicted values
  # ----------------------------------------------------------------------------
  y.hat = x %*% b.hat  
  return( y.hat )
}

rss = function( y, y.hat ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains the RSS defined as 
  #   $$ RSS = \sum_{i=1}^{n} ( y_i - f( x_i ) )^2 $$
  # 
  # ARGUMENTS:
  #   y    : empirical y values 
  #   y.hat: predicted y values  
  #  
  # VALUES:
  #   rss:  returned RSS value
  # ----------------------------------------------------------------------------
  return( sum( y - y.hat )^2 )  
}

aic = function( n, p, rss ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains the AIC defined as 
  #   $$ AIC = n + nlog(2\pi) + nlog(RSS/n ) + 2(p+1) $$
  # 
  # ARGUMENTS:
  #   n  : number of observations  
  #   p  : number of predictors  
  #   rss: residual sum of squares  
  #  
  # VALUES:
  #   AIC:  returned AIC value
  # ----------------------------------------------------------------------------
  
  return( n + n*log( 2 / pi ) + n*log( rss / n ) + 2*( p + 1 ) )
}

bic = function( n, p, rss ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains the AIC defined as 
  #   $$ BIC = n + nlog( 2 \pi ) + nlog( rss / n ) + ( logn )( p + 1 ) $$
  # 
  # ARGUMENTS:
  #   n  : number of observations  
  #   p  : number of predictors  
  #   rss: residual sum of squares  
  #  
  # VALUES:
  #   BIC:  returned BIC value
  # ----------------------------------------------------------------------------
  return( n + n*log( 2 / pi ) + n*log( rss / n ) + ( log(n) * ( p + 1 ) ) )
}

sim.cv.lasso = function( n, p, beta, rho, sigma = 1, N = 1000, K = 5 ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains various diagnostics from a simulation using the lasso 
  #   with K-fold cross-validation
  # 
  # ARGUMENTS:
  #   n:    number of observations
  #   p:    number of predictors
  #   beta: theoretical coefficients of predictors 
  #   rho:  correlation matrix of predictors 
  #   sigma: noise associated with the error term of lasso 
  #   N:    number of simulations 
  #   K:    number of CV folds 
  #  
  # VALUES:
  #   avg.zero: average number of 0 coefficients 
  #   med.mse:  median MSE of all simulations  
  #   avg.s:    average lambda (fraction) chosen  
  # ----------------------------------------------------------------------------
  
  # Initialize Variables
  # --------------------
  return = list()   # List of Simulation Diagnoostics 
  n.zero = c()      # Number of 0 Coefficients in Each Simulation 
  mse = c()    # Vector of MSE 
  rss = c()    # Vector of RSS
  aic = c()    # Vector of AIC 
  bic = c()    # Vector of BIC  
  
  # Run Simulations 
  # ---------------
  for( i in 1:N ){
    
    # Simulate Observations 
    X = mvrnorm( n, rep( 0, p ), rho )  # Generate RVs
    X = scale( X, TRUE, TRUE )          # Standardize RVs
    
    # Calculate Y = B^T x + sigma * e
    Y = apply( X, 1, function( x ) t( beta ) %*% as.matrix( x ) )  
    Y = Y + sigma*rnorm( n )            # Add Noise Term 
    
    
    # Find Regularization Parameter (lambda) with CV 
    fit.cv = cv.lars( X, Y, type = "lasso", plot.it = F, K = K )
    min.cv = which.min( fit.cv$cv )      # Index of Smallest CV Error 
    lambda[i] = fit.cv$index[ min.cv ]   # Lambda with Smallest CV Error  
    
    
    # Get Lasso Estimate at Chosen Lambda   
    fit = lars( X, Y, type = "lasso", normalize = F, intercept = F )
    
    # Extract Coefficients of Chosen Lambda  
    b.hat = predict.lars( fit, s = lambda[i], type = "coefficients", 
                            mode = "fraction" )$coef
    
    # Obtain Diagnostics
    # ------------------
    rss[i] = rss( Y, pred( b.hat, X ) )   # Calculate RSS 
    aic[i] = aic( n, p, rss[i] )          # Calculate AIC 
    bic[i] = bic( n, p, rss[i] )          # Calculate BIC 
    n.zero[i] = sum( b.hat == 0 )         # Number of Nonzero Coefficents 
    
    # Get Mean Squared Error  
    b.error = matrix( ( b.hat - beta ), ncol = 1 )
    mse[i] = t( b.error ) %*% rho %*% b.error
  }
  
  # Save Values into List  
  # ---------------------
  return$avg.zero = sum( n.zero )/N   # Save Avg Number of 0 Coefficients 
  return$med.mse = median( mse )      # Save Median MSE 
  return$avg.s = sum( lambda )/N      # Save Avg Regularization Parameter 
  return$rss = sum( rss )/N     # Save RSS
  return$aic = sum( aic )/N     # Save AIC
  return$bic = sum( bic )/N     # Save BIC  
  return( return )
}

sim.ols = function( n, p, beta, rho, sigma = 1, N = 1000, K = 5 ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains various diagnostics from a simulation using the lasso
  # 
  # ARGUMENTS:
  #   n:    number of observations
  #   p:    number of predictors
  #   beta: theoretical coefficients of predictors 
  #   rho:  correlation matrix of predictors 
  #   sigma: noise associated with the error term of lasso 
  #   N:    number of simulations 
  #   K:    number of CV folds 
  #  
  # VALUES:
  #   avg.zero: average number of 0 coefficients 
  #   med.mse:  median MSE of all simulations  
  #   avg.s:    average lambda (fraction) chosen  
  # ----------------------------------------------------------------------------
  
  # Initialize Variables
  # --------------------
  return = list()   # List of Return Values  
  n.zero = c()      # Number of 0 Coefficients in Each Simulation 
  mse = c()         # Vector of MSE 
  rss = c()         # Vector of RSS
  aic = c()         # Vector of AIC
  bic = c()         # Vector of BIC 
  
  # Run Simulations 
  # ---------------
  for( i in 1:N ){
    
    # Simulate Data  
    # -------------
    X = mvrnorm( n, rep( 0, p ), rho )  # Generate RVs
    X = scale( X, TRUE, TRUE)           # Standardize RVs
    X = matrix( X, nrow = n, ncol = p ) # Convert to Matrix 
    
    # Calculate Y = B^T x + sigma * e
    Y = apply( X, 1, function( x ) t( beta ) %*% as.matrix( x ) )  
    Y = Y + sigma * rnorm( n ) 
    
    
    # Obtain Model
    # ------------
    fit = lm( Y ~ ., as.data.frame( X ) )
    b.hat = coefficients( fit )[-1]     # Extract Coefficients 
    
    
    # Obtain Diagnostics
    # ------------------
    rss[i] = rss( Y, pred( b.hat, X ) ) # Calculate RSS
    aic[i] = aic( n, p, rss[i] )        # Calculate AIC
    bic[i] = bic( n, p, rss[i] )        # Calculate BIC
    n.zero[i] = sum( b.hat == 0 )       # Calculate Number of Zero Coefficients 
    
    # Get Mean Squared Error  
    b.error = matrix( ( b.hat - beta ), ncol = 1 )
    mse[i] = t( b.error ) %*% cov(X) %*% b.error
  }
  
  # Save Diagnostics  
  # ----------------
  return$avg.zero = sum( n.zero )/N   # Save Avg Number of 0 Coefficients 
  return$med.mse = median( mse )      # Save Median MSE 
  return$avg.s = NA                   # OLS Doesn't Rely on Regularization  
  return$rss = sum( rss )/N   # Save RSS
  return$aic = sum( aic )/N   # Save BIC
  return$bic = sum( bic )/N   # Save BIC
  
  return( return )
}

sim.ridge = function( n, p, beta, rho, sigma = 1, N = 1000, K = 5 ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains various diagnostics from a simulation using ridge 
  #   regression with K-fold cross-validation
  # 
  # ARGUMENTS:
  #   n:    number of observations
  #   p:    number of predictors
  #   beta: theoretical coefficients of predictors 
  #   rho:  correlation matrix of predictors 
  #   sigma: noise associated with the error term of model 
  #   N:    number of simulations 
  #   K:    number of CV folds 
  #  
  # VALUES:
  #   avg.zero: average number of 0 coefficients 
  #   med.mse:  median MSE of all simulations  
  #   avg.s:    average lambda (fraction) chosen  
  # ----------------------------------------------------------------------------
  
  # Initialize Variables
  # --------------------
  return = list()   # List of Simulation Diagnoostics 
  n.zero = c()      # Number of 0 Coefficients in Each Simulation 
  mse = c()         # Vector of MSE 
  
  # Run Simulations 
  # ---------------
  for( i in 1:N ){
    
    # Simulate Observations 
    X = mvrnorm( n, rep( 0, p ), rho )  # Generate RVs
    X = scale( X, TRUE, TRUE )          # Standardize RVs
    X = matrix( X, nrow = n, ncol = p ) # Convert to Matrix 
    
    # Calculate Y = B^T x + sigma * e
    Y = apply( X, 1, function( x ) t( beta ) %*% as.matrix( x ) )  
    Y = Y + sigma*rnorm( n )    # Add Noise Term 
    
    # Find Best Regularization Parameter
    fit.cv = ridge.cv( X, Y, k = K )
    
    # Extract Coefficients ( Default Excludes Intercept ) 
    b.hat = coefficients( fit.cv )
    
    # Obtain Diagnostics 
    # ------------------
    rss[i] = rss( Y, pred( b.hat, X ) ) # Calculate RSS
    aic[i] = aic( n, p, rss[i] )        # Calculate AIC
    bic[i] = bic( n, p, rss[i] )        # Calculate BIC
    n.zero[i] = sum( b.hat == 0 )       # Number of Zero Coefficients 
    
    # Get Mean Squared Error  
    b.error = matrix( ( b.hat - beta ), ncol = 1 )
    mse[i] = t( b.error ) %*% cov(X) %*% b.error
  
  }
  
  # Save Values into List  
  # ---------------------
  return$avg.zero = sum( n.zero )/N   # Save Avg Number of 0 Coefficients 
  return$med.mse = median( mse )      # Save Median MSE 
  return$avg.s = NA                   # Save Avg Regularization Parameter 
  return$rss = sum( rss )/N     # Save RSS
  return$aic = sum( aic )/N     # Save AIC
  return$bic = sum( bic )/N     # Save BIC
  return( return )
}

sim.gcv.lasso = function( n, p, beta, rho, sigma = 1, N = 1000, K = 5 ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains various diagnostics from a simulation using the lasso
  #   with generalized cross-validation. 
  # 
  # ARGUMENTS:
  #   n:    number of observations
  #   p:    number of predictors
  #   beta: theoretical coefficients of predictors 
  #   rho:  correlation matrix of predictors 
  #   sigma: noise associated with the error term of lasso 
  #   N:    number of simulations 
  #   K:    number of CV folds 
  #  
  # VALUES:
  #   avg.zero: average number of 0 coefficients 
  #   med.mse:  median MSE of all simulations  
  #   avg.s:    average lambda (fraction) chosen  
  # ----------------------------------------------------------------------------
  
  # Initialize Variables
  # --------------------
  return = list()   # List of Simulation Diagnoostics 
  n.zero = c()      # Number of 0 Coefficients in Each Simulation 
  mse = c()         # Vector of MSE 
  gcv = c()         # Vector of GCV
  rss = c()         # Vector of RSS
  aic = c()         # Vector of AIC
  bic = c()         # Vector of BIC 
  
  # Run Simulations 
  # ---------------
  for( i in 1:N ){
    
    # Simulate Observations 
    X = mvrnorm( n, rep( 0, p ), rho )  # Generate RVs
    X = scale( X, TRUE, TRUE )          # Standardize RVs
    
    # Calculate Y = B^T x + sigma * e
    Y = apply( X, 1, function( x ) t( beta ) %*% as.matrix( x ) )  
    Y = Y + sigma*rnorm( n )    # Add Noise Term 
    
    # Find Regularization Parameter( lambda ) with GCV
    fit = lars( X, Y, type = "lasso", normalize = F, intercept = F )
    
    # Find Best GCV 
    for( j in 1:100 ){
      
      # Get Model at Various Lambda 
      fit.coef = predict( fit, s = j/100, type = "coefficients", 
                          mode = "fraction" )$coef    
      
      # Calculate GCV
      fit.rss = rss( y, pred( b.hat, x ) )    # RSS of Current Model 
      P = sum( fit.coef = 0 )                 # Number of Degrees of Freedom 
      gcv[ j ] = 1/n * fit.rss / ( 1 - P / n )^2 
      
    }
    # Get Minimum GCV to Choose Lambda 
    min.gcv = which.min( gcv )
    lambda[ i ] = min.gcv / 100 
    
    # Extract Coefficients of Chosen Lambda  
    b.hat = predict.lars( fit, s = lambda[i], type = "coefficients", 
                            mode = "fraction" )$coef
    
    # Obtain Diagnostics
    # ------------------
    rss[i] = rss( Y, pred( b.hat, X ) )   # Calculate RSS 
    aic[i] = aic( n, p, rss[i] )          # Calculate AIC 
    bic[i] = bic( n, p, rss[i] )          # Calculate BIC 
    n.zero[i] = sum( b.hat == 0 )         # Number of Zero Coefficients 
    
    # Get Mean Squared Error  
    b.error = matrix( ( b.hat - beta ), ncol = 1 )
    mse[i] = t( b.error ) %*% cov(X) %*% b.error
  
  }
  
  # Save Values into List  
  # ---------------------
  return$avg.zero = sum( n.zero )/N   # Save Avg Number of 0 Coefficients 
  return$med.mse = median( mse )      # Save Median MSE 
  return$avg.s = sum( lambda )/N      # Save Avg Regularization Parameter 
  return$rss = sum( rss )/N     # Save RSS
  return$aic = sum( aic )/N     # Save AIC 
  return$bic = sum( bic )/N     # Save BIC 
  return( return )
}

sim.subset = function( n, p, beta, rho, sigma = 1, N = 1000, K = 5 ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains various diagnostics from a simulation using the best 
  #   subset selection 
  # 
  # ARGUMENTS:
  #   n:    number of observations
  #   p:    number of predictors
  #   beta: theoretical coefficients of predictors 
  #   rho:  correlation matrix of predictors 
  #   sigma: noise associated with the error term of model  
  #   N:    number of simulations 
  #   K:    number of CV folds 
  #  
  # VALUES:
  #   avg.zero: average number of 0 coefficients 
  #   med.mse:  median MSE of all simulations  
  #   avg.s:    average lambda (fraction) chosen  
  # ----------------------------------------------------------------------------
  
  # Initialize Variables
  # --------------------
  return = list()   # List of Simulation Diagnoostics 
  n.zero = c()      # Number of 0 Coefficients in Each Simulation 
  mse = c()         # Vector of MSE 
  
  # Run Simulations 
  # ---------------
  for( i in 1:N ){
    
    # Simulate Observations 
    X = mvrnorm( n, rep( 0, p ), rho )  # Generate RVs
    X = scale( X, TRUE, TRUE )          # Standardize RVs
    X = matrix( X, nrow = n, ncol = p ) # Convert to Matrix 
    
    # Calculate Y = B^T x + sigma * e
    Y = apply( X, 1, function( x ) t( beta ) %*% as.matrix( x ) )  
    Y = Y + sigma*rnorm( n )    # Add Noise Term 
    
    # Get Best Subset 
    fit = bestglm( data.frame( Y, X ), IC = "CV" )
    
    # Extract Coefficients of Chosen Lambda  
    b.hat = predict.lars( fit, s = lambda[i], type = "coefficients", 
                            mode = "fraction" )$coef
    
    # Get Number of 0 Coefficients 
    n.zero[i] = sum( b.hat == 0 ) 
    
    # Get Mean Squared Error  
    b.error = matrix( ( b.hat - beta ), ncol = 1 )
    mse[i] = t( b.error ) %*% cov(X) %*% b.error
  
  }
  
  # Save Values into List  
  # ---------------------
  return$avg.zero = sum( n.zero )/N   # Save Avg Number of 0 Coefficients 
  return$med.mse = median( mse )      # Save Median MSE 
  return$avg.s = sum( lambda )/N      # Save Avg Regularization Parameter 
  return( return )
}

sim.step = function( n, p, beta, rho, direction, sigma = 1, N = 1000, K = 5 ){
  # ----------------------------------------------------------------------------
  # DESCRIPTION: Obtains various diagnostics from a simulation using the lasso
  # 
  # ARGUMENTS:
  #   n:    number of observations
  #   p:    number of predictors
  #   beta: theoretical coefficients of predictors 
  #   rho:  correlation matrix of predictors 
  #   sigma: noise associated with the error term of lasso 
  #   N:    number of simulations 
  #   K:    number of CV folds 
  #  
  # VALUES:
  #   avg.zero: average number of 0 coefficients 
  #   med.mse:  median MSE of all simulations  
  #   avg.s:    average lambda (fraction) chosen  
  # ----------------------------------------------------------------------------
  
  # Initialize Variables
  # --------------------
  return = list()   # List of Return Values  
  n.zero = c()      # Number of 0 Coefficients in Each Simulation 
  mse = c()         # Vector of MSE 
  
  # Run Simulations 
  # ---------------
  for( i in 1:N ){
    
    # Simulate Observations 
    X = mvrnorm( n, rep( 0, p ), rho )  # Generate RVs
    X = scale( X, TRUE, TRUE)           # Standardize RVs
    X = matrix( X, nrow = n, ncol = p ) # Convert to Matrix 
    
    # Calculate Y = B^T x + sigma * e
    Y = apply( X, 1, function( x ) t( beta ) %*% as.matrix( x ) )  
    Y = Y + sigma * rnorm( n ) 
    
    # Get Null and Full Estimates
    null = lm( Y ~ 1, as.data.frame( X ) )
    full = lm( Y ~ ., as.data.frame( X ) )
    
    # Get Stepwise Regression
    if( direction == "forward" ){
      fit = step( null, scope = list( null, full ), direction = "forward" )
      
    } else if( direction == "backward" ){
      fit = step( full, data = X, direction = "backward" )  
    }
    
    # Extract Coefficients 
    b.hat = coefficients( fit )[-1] 
      # NOTE: MUST BE CHANGED SO THAT FUNCTIONS ARE MATCHED 
    
    # Get Number of 0 Coefficients 
    n.zero[i] = sum( b.hat == 0 ) 
    
    # Get Mean Squared Error  
    b.error = matrix( ( b.hat - beta ), ncol = 1 )
    mse[i] = t( b.error ) %*% cov(X) %*% b.error
#    mse[i] = t( b.error ) %*% rho %*% b.error 
  
  }
  
  # Save Values into List  
  # ---------------------
  return$avg.zero = sum( n.zero )/N   # Save Avg Number of 0 Coefficients 
  return$med.mse = median( mse )      # Save Median MSE 
  return$avg.s = NA  
  
  return( return )
}


# Example 1
# =========
n = 20      # Number of Observations
p = 8       # Number of Variables 
beta = matrix( c( 3, 1.5, 0, 0, 2, 0, 0, 0 ), ncol = 1 )
cor = 0.5 
rho = get.cor( p, cor )

sim.cv.lasso( n, p, beta, rho, sigma = 3, N = 1000, K = 5 ) 
sim.ols( n, p, beta, rho, sigma = 3, N = 1000, K = 5 ) 
sim.ridge( n, p, beta, rho, sigma = 3, N = 1000, K = 5 ) 

# Example 2 
# =========
n = 20      # Number of Observations
p = 8       # Number of Variables 
beta = matrix( rep( 0.85, p ), ncol = 1 )
cor = 0.5 
rho = get.cor( n, cor )

sim.cv.lasso( n, p, beta, rho, sigma = 3, N = 1000, K = 5 ) 
sim.ols( n, p, beta, rho, sigma = 3, N = 1000, K = 5 ) 
sim.ridge( n, p, beta, rho, sigma = 3, N = 1000, K = 5 ) 

# Example 3
# =========
n = 20      # Number of Observations
p = 8       # Number of Variables 
beta = matrix( c( 5, 0, 0, 0, 0, 0, 0, 0 ), ncol = 1 )
cor = 0.5 
rho = get.cor( p, cor )

sim.cv.lasso( n, p, beta, rho, sigma = 2, N = 1000, K = 5 ) 
sim.ols( n, p, beta, rho, sigma = 2, N = 1000, K = 5 ) 
sim.ridge( n, p, beta, rho, sigma = 2, N = 1000, K = 5 ) 

