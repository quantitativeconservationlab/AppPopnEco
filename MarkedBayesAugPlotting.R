######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz to plot results ########
### from Bayesian models of capture recapture including data #
#### agumentation, which allows detection to vary among individuals #
#####################################################################
##### Set up your workspace and load relevant packages -----------
# load packages:
library( tidyverse ) 
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI )
source("https://raw.githubusercontent.com/datavizpyr/data/master/half_flat_violinplot.R")

###################################################################
#### Load or create data -----------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
#load relevant workspace
load( "MarkedAugResultsJAGs.RData" )
################################################################################
#################### viewing model results ####################################
##################################################################################

#define model results to plot
mr <- ma1 #mht
summary(mr)

############## trace plots ############
#plot( mr ) #plots traces and posterior densities for all parameters
par( mfrow = c( 2, 2 ), ask = F, mar = c(3,4,2,2) )
#detection parameters
#intercept 
traceplot( mr, parameters = c( 'int.det') )
#random intercepts for observer
traceplot( mr, parameters = c( 'int.lam') )
#fixed effect
traceplot( mr, parameters = c( 'alpha') )
#for abundance
traceplot( mr, parameters = c( 'beta') )
#random site effects error
traceplot( mr, parameters = c( 'sigma.i', 'sigma.n') )
#random intercept for site
traceplot( mr, parameters = c("eps" ))
# random intercept for individual 
traceplot( mr, parameters = c( 'eta') )

############## whisker plots #############
par( mfrow = c( 2,1 ), ask = F , mar = c(3,4,2,2) )
#for detection
#fixed effects for detection and abundance
whiskerplot( mr, parameters = c( 'alpha', "beta" ) , 
             zeroline = TRUE )
#error for random effects
whiskerplot( mr, parameters = c( 'sigma.i', "sigma.n" ) )
#random effects for site
whiskerplot( m1, parameters = c( "eps" ) )
#random effects for individual
whiskerplot( m2, parameters = c( "eta" ) )
# probability of detection
whiskerplot( m1, parameters = c( "p" ) )

##################################################################
#### visualizing abundance for sites trapped ###################
###
dim( mr$sims.list$N)
#levels: runs, sites, individuals
head( ij_wide )
#create dataframe to store abundance estimates
N_df <- ij_wide %>% 
  dplyr::select( idno, id, SiteID, year )
#add columns for abundance
N_df$N.mean <- N_df$N.lowCI <- N_df$N.highCI <- NA

#number of runs
R <- mr$mcmc.info$n.samples
#create matrix
a <- matrix(data = NA, nrow = R, ncol = I)
#loop over runs and sites
for( i in 1:I){
  for( r in 1:R){
    #sum across individuals
    a[r,i] <- sum( mr$sims.list$N.site[r,i,1:M] )
  }}
for( i in 1:I){
  #average across runs
    N_df[i,"N.mean"] <- mean( a[,i],na.rm = TRUE )
    N_df[i,c("N.lowCI", "N.highCI") ] <- quantile( a[,i],
                                    probs = c(0.025, 0.975),
                                    na.rm = TRUE )
}
#view abundance estimates
head(N_df)
###### plot abundance
#plot abundance between sites and years
N_df %>% 
  ggplot(., aes( y = SiteID, x = N.mean, 
                 color = as.factor(year) ) ) +
  theme_classic( base_size = 12 ) +
  theme( legend.position = "top" ) +
  labs( x = "Annual abundance", 
        y = "Site" ) +
  geom_point( size = 3 ) +
  geom_errorbar( aes( xmin = N.lowCI, xmax = N.highCI ), 
                 size = 2 ) 

######### partial prediction plots #############################
#Define plotting function: ##
estppreds <- function( type = "detection", sl = 100, int, coefs, 
                       rawpreds, labs, nicelabs ){
  # type: detection or abundance
  # sl: length of predictions
  # int: vector of 1s of sl length
  # coefs: fixed betas
  # rawpreds: dataframe containing predictors we want to predict over 
  #on the real scale
  # labs  = column labels in same order as betas
  # nicelabs = if you want to rename include nice labels here 
  
  #intercept for covariate matrix 
  ones <- rep( 1, sl )
  
  #loop over predictors of interest
  for( mp in 1:length( labs ) ){
    # create predicted range for predictor on real scale 
    x <- seq( min( rawpreds[ ,labs[ mp ]  ], na.rm = TRUE ),
              max( rawpreds[ , labs[ mp ]], na.rm = TRUE ),
              length.out = sl )
    #scale predictor
    sclx <- scale( x )
    
    ifelse( length(labs) == 1, 
            fixed <- cbind( int, coefs),
            fixed <-  cbind( int, coefs[,mp] ) )
    ifelse( type == "detection",
            #estimate predicted response
            estpred <- plogis( fixed %*% t( cbind(ones, sclx ) ) ),
            estpred <- exp( fixed %*% t(cbind( ones, sclx) ) ) )
    #calculate mean response:
    m <- apply( estpred, MARGIN = 2, FUN = mean )
    #calculate 95% CIs of response:
    CI <- apply( estpred, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975) )
    #create dataframe with results:
    df <- data.frame( sclx, x, m, t( CI ), 
                      rep( nicelabs[mp], length( m ) ))
    colnames( df ) <- c( "Scaled", "Raw",
                         "Mean", "lowCI", "highCI", "Predictor" )
    #check that it is working
    print( head(df) )
    #create dataframe for 1st species, append data for remaining species:
    ifelse( mp == 1, outdf <- df, outdf <- bind_rows( outdf, df ) )
  } 
  #output of function:
  return( outdf )
}

lam.preds <- estppreds( type = 'abundance',
                        sl = 100, 
                        int = mr$sims.list$int.lam, 
                        coefs = mr$sims.list$beta,
                        rawpreds = ik_df, 
                        labs = "herbaceous",
                        nicelabs = "Herbaceous (%)"
)
#plot predictors for abundance
ggplot( data = lam.preds, aes( x = Raw, y = Mean ) ) +
  theme_classic() +
  ylab( "Ground squirrel abundance" ) +
  xlab( "" ) +
  theme( legend.position = "top", 
         legend.title = element_blank(),
         text = element_text( size = 19 ), 
         axis.line = element_line( size = 1.5 ),
         panel.spacing = unit(1, "lines"),
         strip.background = element_blank() ) +
  geom_line( size = 1.5 ) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI , ymax = highCI ) ) + 
  facet_wrap( ~ Predictor,
              scales = "free_x", ncol = 3,
              strip.position = "bottom" )

p.preds <- estppreds( type = 'detection',
                      sl = 100, 
                      int = mr$sims.list$int.det, 
                      coefs = mr$sims.list$alpha,
                      rawpreds = ikj_df, 
                      labs = c("wind_kmph_st", "tempC_st", "effort_mins" ), 
                      nicelabs = c( "Wind (km/hr)", "Mean temperature (C)",
                                    "Effort (mins)" )
)

ggplot( data = p.preds, aes( x = Raw, y = Mean ) ) +
  theme_classic() +
  ylab( "Detection probability" ) +
  xlab( "" ) +
  theme( legend.position = "top", 
         legend.title = element_blank(),
         text = element_text( size = 19 ), 
         axis.line = element_line( size = 1.5 ),
         panel.spacing = unit(1, "lines"),
         strip.background = element_blank() ) +# ylim( 0.0, 0.2 ) +
  geom_line( linewidth = 1.5 ) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI , ymax = highCI ) ) + 
  facet_wrap( ~ Predictor,
              scales = "free_x", ncol = 3,
              strip.position = "bottom" )

###How do these plots differ from the ones without data augmentation?
# Answer:
# 
#

####################### end of script ################################# 
