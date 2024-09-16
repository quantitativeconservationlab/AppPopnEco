###########################################################################
### This script was developed by Jen Cruz and David Bontrager         ###
####  for use in the Applied Population Ecology class                 ###
###                                                                   ###
### Script import results from OccBayesAnalysis.R estimating          ###
### single season occupancy of owls in Coastal Texas in relation      ###
### to habitat. Here we plot model results                              #
###########################################################################

########### clean workspace and load required packages ####################
###########################################################################
####### load relevant packages ###
library( tidyverse ) #dataframe manipulations.
# set option to see all columns and more than 10 rows
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI ) #to view and run RJAGS code
library( gridExtra ) #plot multiple ggplots
library( ggrepel ) #non overlapping graph labels
#######################    import relevant data   ##############
#####clean workspace to improve efficiency: ###
rm(list = ls() ) 
# load relevant workspace
load( "OccBayesResults.RData" )

################################################################################
#################### viewing model results ####################################
##################################################################################
# Define model results to plot
# Choose either m_barn or m_great
mr <- m_great

#view summary of results
mr

############## trace plots ############
#plot( mr ) #plots traces and posterior densities for all parameters
par( mfrow = c( 3, 3 ), ask = F, mar = c(3,4,2,2) )
#intercept in occupancy submodel
traceplot( mr, parameters = c( 'int.psi') )
#coefficients in occupancy submodel
traceplot( mr, parameters = c( 'beta.psi') )
#intercept in detection submodel
traceplot( mr, parameters = c( 'int.p') )
#coefficients in detection submodel
traceplot( mr, parameters = c( 'alpha.p') )

############## whisker plots #############
par( mfrow = c( 1,1 ), ask = F , mar = c(3,4,2,2) )
#coefficients for occupancy
whiskerplot( mr, parameters = c( 'beta.psi' ) , zeroline = TRUE )
#coefficients for detection
whiskerplot( mr, parameters = c( 'alpha.p' ), zeroline = TRUE )
#intercepts
whiskerplot( mr, parameters = c( 'int.psi', 'int.p' ), zeroline = TRUE )
#derived parameters
whiskerplot( mr, parameters = c( "psi" ) )
whiskerplot( mr, parameters = c( "p" ) )

#############################################################
######### partial prediction plots #############################
# Estimate partial prediction plots (marginal effect plots) for # 
# predictors with 95% CIs not overlapping zero. #

#From the quick look at the whiskerplots...our detection #
#predictors were not that important. We choose to plot #
# the relationship with day of survey as a way of demonstrating #
# the process
head( detdf)
#how many predictions do we want?
n <- 100
# what are the min max wind days:
detdf %>% select( wind1 , wind2 , wind3  ) %>% 
  summarise_all( list(min, max ), na.rm = TRUE )
#use these to define your bounds
wind.pred <- seq( 0, 35, length.out = n )
wind.std <- scale( wind.pred )

#extract detection intercept and relevant coefficient 
fixeddet <- cbind( mr$sims.list$int.p, mr$sims.list$alpha )
#estimate predicted detection
preddet <- plogis( fixeddet %*% t( cbind( rep(1,n), wind.std ) ) )
#calculate mean abundance
mdet <- apply( preddet, MARGIN = 2, FUN = mean )
#calculate 95% Credible intervals for abundance
CIdet <- apply( preddet, MARGIN = 2, FUN = quantile, 
                probs = c(0.025, 0.975) )

#create dataframe combining all predicted values for plotting
p.pred <- data.frame( mdet, t(CIdet),
                     wind.pred, wind.std )
#view
head( p.pred); dim( p.pred )
#rename predicted abundance columns
colnames(p.pred )[1:4] <- c(  "Mean", "lowCI", "highCI",
                             "Wind" )


#plot
ggplot( data = p.pred, aes( x = Wind, y = Mean ) ) +
  theme_classic() +
  ggtitle( "Great Horned Owl: Effects on Detection Probability" ) + # change this for different owl species
  ylab( "Probability of Detection" ) +
  xlab( "Wind speed (miles/hr)" ) +
  theme( legend.position = "top",
         legend.title = element_blank(),
         text = element_text( size = 16 ),
         axis.line = element_line( size = 1.3 ),
         panel.spacing = unit(1, "lines"),
         strip.background = element_blank() ) +
  geom_line( linewidth = 1.2 ) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI , ymax = highCI ) ) 


######## For occupancy predictors we create a function #####
## that allows us to plot all predictors at once as well ####
### as accounting for the different intercepts for year   ###

#Start by creating a vector of predictor labels. Ensure they appear in the 
#same order as the coefficients in submodels
#column labels of predictors used in the model
psilabs <- c("shrub", "ah", "aquatic" )
#Nice labels you want in your plots
nicepsilabs <- c( 'Shrub Cover (%)',
                  "Aplomado Habitat Cover (%)", "Woody wetland (%)" )

#Define function
estppreds <- function( sl = 100, int, coefs, 
                       datadf, labs, nicelabs ){
  # sl: length of predictions
  # int: vector of 1s of sl length
  # coefs: fixed coefficients 
  # datadf: dataframe containing predictors we want to predict over 
  # on the real scale
  # labs  = column labels in same order as betas
  # nicelabs = if you want to rename include nice labels here 
  
  #intercept for covariate matrix 
  ones <- rep( 1, sl )
 
  #loop over predictors of interest
  for( mp in 1:length( labs ) ){
    # create predicted range for predictor on real scale 
    x <- seq( min( datadf[ ,labs[ mp ] ], na.rm = TRUE ),
              max( datadf[ , labs[ mp ] ], na.rm = TRUE ),
              length.out = sl )
    #scale predictor
    sclx <- scale( x )
    # extract fixed effects for year 1
    fixed <-  cbind( int[,1], coefs[,mp] )
    #estimate predicted response for year 1
    estpred <- plogis( fixed %*% t( cbind(ones, sclx ) ) )
    #calculate mean response:
    m <- apply( estpred, MARGIN = 2, FUN = mean )
    #calculate 95% CIs of response:
    CI <- apply( estpred, MARGIN = 2, FUN = quantile, 
                 probs = c(0.025, 0.975) )
    #create dataframe with results:
    df1 <- data.frame( sclx, x, m, t( CI ), 
                      rep( nicelabs[mp], length( m ) ), 
                      rep( "year 1", length(m)) )
    #for year two
    fixed2 <-  cbind( int[,2], coefs[,mp] )
    #estimate predicted response for year 2
    estpred2 <- plogis( fixed2 %*% t( cbind(ones, sclx ) ) )
    #calculate mean response:
    m2 <- apply( estpred2, MARGIN = 2, FUN = mean )
    #calculate 95% CIs of response:
    CI2 <- apply( estpred2, MARGIN = 2, FUN = quantile, 
                 probs = c(0.025, 0.975) )
    df2 <- data.frame( sclx, x, m2, t( CI2 ), 
                      rep( nicelabs[mp], length( m ) ), 
                      rep( "year 2", length(m)) )
    #label columns
    colnames( df1 ) <- colnames( df2 ) <- c( "Scaled", "Raw",
               "Mean", "lowCI", "highCI", "Predictor", "Year" )
    #join dataframes
    df <- bind_rows( df1, df2 )
    #label columns
    colnames( df ) <- c( "Scaled", "Raw",
                "Mean", "lowCI", "highCI", "Predictor", "Year" )
    
    #check that it is working
    print( head(df) )
    #extend dataframe if more than one predictor
    ifelse( mp == 1, outdf <- df, 
            outdf <- bind_rows( outdf, df ) )
  } 
  #output of function:
  return( outdf )
}


#run it for occupancy
occ.preds <- estppreds( 
  sl = 100, 
  int = mr$sims.list$int.psi, 
  coefs = mr$sims.list$beta.psi, 
  datadf = habdf, 
  labs = psilabs, 
  nicelabs = nicepsilabs
)

#plot
ggplot( data = occ.preds, aes( x = Raw, y = Mean, 
         fill = as.factor(Year) ) ) +
  theme_classic() +
  ggtitle( "Great Horned owls: Habitat Effects on Occupancy" ) + 
  ylab( "Probability of Occupancy" ) +
  xlab( "" ) +
  theme( legend.position = "top", 
         legend.title = element_blank(),
         text = element_text( size = 16 ), 
         axis.line = element_line( size = 1.3 ),
         panel.spacing = unit(1, "lines"),
         plot.title = element_text( hjust = 0.5 ),
         strip.background = element_blank() ) + ylim( 0.0, 1.0 ) +
  geom_line( size = 1.2 ) +
  geom_ribbon( alpha = 0.3, aes( ymin = lowCI , ymax = highCI ) ) + 
  facet_wrap( ~ Predictor, 
              scales = "free_x", ncol = length(psilabs),
              strip.position = "bottom" )#

######### save relevant output   ###########################

# How would you save one of your figures?
# Answer:
# 

################ end of script  ##############################