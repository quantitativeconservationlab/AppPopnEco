
####### load relevant packages ###
library( tidyverse ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI ) #to run RJAGS

#######################    import relevant data   ##############
#clean workspace to improve efficiency:
rm(list = ls() ) 
#load relevant workspace 
load( "OccBayesResults.RData" )

################################################################################
#################### viewing model results ####################################
##################################################################################
#define model object
mr <- m_great

#Total number of iterations ran:
N <- mr$mcmc.info$n.samples

dim(mr$sims.list$lik_yobs)
dim(mr$sims.list$lik_yhat)
########################################################################
## Calculating deviance from likelihood of data conditional on model ###
########################################################################
#pick columns of interest from datadf
evaldf <- detdf %>% 
    dplyr::select( Site.ID, yearid, siteyear_id )
head( evaldf )
#add likelihood values

# Define function to calculate Model deviances 
CalcDevs <- function( lik_yobs = mr$sims.list$lik_yobs,
                      lik_yhat = mr$sims.list$lik_yhat ){
  
  #assign temporary objects
  ll_yobs_i <-  ll_yhat_i <-  matrix( nrow = N, ncol = I )
  Dev_obs_i <<- Dev_hat_i <<- matrix( nrow = N, ncol = I )
  ll_yobs <- ll_yhat <- Dev_obs <- Dev_hat <- rep( NA, N )

#  print( dim( ll_yobs_i ) )
  
  for( n in 1:N ){ #loop over mcmc iterations
    for( i in 1:I ){ #loop over sites
      # print( sum( log( lik_yobs[ n,i, ] ) ) )
       # sum log likelihoods for that site
        ll_yobs_i[ n,i ] <-  sum( log( lik_yobs[ n,i, ] ) )
        ll_yhat_i[ n,i ] <- sum( log( lik_yhat[ n, i, ] ) ) 
        # Calculate species-level deviances
        Dev_obs_i[ n, i ] <<- -2 * ( ll_yobs_i[ n, i ] )
        Dev_hat_i[ n, i ] <<- -2 * ( ll_yhat_i[ n, i ] )  
      } #close I loop
      
    #sum log likelihoods across each mcmc iteration
    #for observed y
    ll_yobs[ n ] <- sum( log( lik_yobs[ n,, ] ) ) 
    #for estimated y
    ll_yhat[ n ] <- sum( log( lik_yhat[ n,, ] ) ) 
    # Calculate overall deviances
    Dev_obs[ n ] <- -2 * ll_yobs[ n ]
    Dev_hat[ n ] <- -2 * ll_yhat[ n ] 
    
  } # close n loop
  
  return( data.frame( Dev_obs, Dev_hat, ll_yobs, ll_yhat ) )
  
} # close function

#######################################################################
# Calculate Bayesian p-value as mean number of times that Deviance of #
# observed data was greater than the deviance of predicted model #####
#########################################################################

Baypvalue <- mean( ModDevs$Dev_obs > ModDevs$Dev_hat )
print( Baypvalue )

#plot Bayesian p value for the model
par( mfrow = c(1,1), cex.axis = 1.5, mar=c(3,4,2,2) )
plot( ModDevs$Dev_obs, ModDevs$Dev_hat, 
      main = paste0( "Bayesian P value = ", 
                     round( Baypvalue, 3 ) ), 
      xlab = "Deviance of observed data", 
      ylab = "Deviance of predicted data", tcl = 0.2, 
      cex.lab = 1.5, bty = "l"  ) 
abline( 0, 1, col = 'black', lwd = 3 )

#############       END OF SCRIPT         ###########################