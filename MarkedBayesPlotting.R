######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz to plot results ########
### from Bayesian models of capture recapture resembling M[t]    ###
### models that include N as a derived parameter (m1) or one that #
### estimates N conditional on n (observed number of individuals) (m2) #
# both models contain the same predictors for detection and abundance #
# and include a site random intercept in abundance to account for #
# possible overdispersion. #
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
load( "MarkedResultsJAGs2.RData" )
################################################################################
#################### viewing model results ####################################
##################################################################################

#start by quickly comparing model results
summary(m1); summary(m2)

############## whisker plots #############
par( mfrow = c( 2,1 ), ask = F , mar = c(3,4,2,2) )
#for detection
#fixed effects for detection and abundance
whiskerplot( m1, parameters = c( 'alpha', "beta" ) , 
             zeroline = TRUE )

whiskerplot( m2, parameters = c( 'alpha', "beta" ) , 
             zeroline = TRUE )

whiskerplot( m1, parameters = c( "p" ) )
whiskerplot( m2, parameters = c( "p" ) )
whiskerplot( m1, parameters = c( "N" ) )
whiskerplot( m2, parameters = c( "N" ) )
whiskerplot( m1, parameters = c( "eps.i" ) )
whiskerplot( m2, parameters = c( "eps.i" ) )

# calculate Bayesian p values 
#m1 model:
mean(m1$sims.list$fit1_hat > m1$sims.list$fit1 )
#m2 model
mean(m2$sims.list$fit1_hat > m2$sims.list$fit1 )
mean(m2$sims.list$fit2_hat > m2$sims.list$fit2 )

#choose model going forward
mr <- m1
#derived parameters

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

# ### if model does not rely on data augmentation 
for(i in 1:I){
  N_df[i,'N.mean'] <- mean( mr$sims.list$N[,i] )
  N_df[i,c('N.lowCI', 'N.highCI')] <- as.vector( quantile(
    mr$sims.list$N[,i ],probs = c(0.025, 0.975) ) )
}

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
    
    # extract fixed effects 
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
                        nicelabs = c("Herbaceous cover (%)")
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

#### Now replot these for the other model that we run. How do they differ?
# Answer:
# 
#

####################### end of script ################################# 
