#######################################################################
#######################################################################
##     This script was created by Dr. Jen Cruz as part of            ##
##            the Applied Population Ecology Class                  ###
##                                                                   ##  
## Here we import our cleaned data for 2009 for the point count      ##
#  observations for Piute ground squirrels at the NCA and run a      ##
## closed population N-mixture analysis. The model is hierarchical    #
#  with : (1) an ecological submodel linking abundance to             #
## environmental predictors at each site; (2) an observation submodel #
## linking detection probability to relevant predictors.             ##
##                                                                   ##
# Abundance is expected to be higher in sites with more sagebrush     #
# and lower in those with more cheatgrass.                            #                                        #
# Detection may be related to observer effects and to time of day     #
#######################################################################

##### Set up your workspace and load relevant packages -----------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# Check that you are in the right project folder
getwd()

#install relevant packages
install.packages("jagsUI")

# load packages:
library( tidyverse ) 
options( dplyr.width = Inf, dplyr.print_min = 100 )
library( jagsUI )
###################################################################
#### Load or create data -----------------------------------------
# set directory where your data are:
datadir <- paste( getwd(), "/Data/", sep = "" )
# load cleaned data
closeddf <- read.csv( file = paste( datadir, "closed_counts.csv", sep = ""),
                      header = TRUE )
#view
head( closeddf ); dim( closeddf ) 
#### End of data load -------------
####################################################################
##### Ready data for analysis --------------
# extract broad parameters of interest
#number of sites:
I <- length( closeddf$o.sites )
#number of visits
J = 3

#Create a new dataframe to scale predictors starting with site level predictors #
XI <- closeddf %>% 
  select( cheatgrass, sagebrush ) %>% 
  mutate( cheatgrass = scale( cheatgrass),
          sagebrush = scale( sagebrush) )

#for sitexvisit predictor we create separate objects for each 
#create index of columns of interest
tidx <- grep("time", colnames( closeddf), value = FALSE)
oidx <- grep( "obs", colnames(closeddf), value = FALSE)
#also for our observations
yidx <- grep( "count.j", colnames( closeddf ), value = FALSE)
#now use those indices to automatically select correct columns and scale
time_sc <- scale( closeddf[ ,tidx] )

# for observers, we will include them as a random intercept instead
# of a fixed effect so we convert them to a number to represent index
obvs <- as.matrix( closeddf[,oidx] )
obvs <- as.numeric( as.factor(obvs) )
#turn back into matrix
obvs <- structure( obvs, dim = dim( closeddf[,oidx] ), class = "matrix" )
#number of observers 
S <- max(obvs, na.rm = TRUE)

##########################################################################
####### m1 single season N-mixture abundance model    #
# ecological predictors: cheatgrass and sagebrush #
# detection predictors: time and observer as random intercept #
############################################################################
############## Specify model in bugs language:  #####################
sink( "m1.txt" )
cat( "
     model{
     
      #priors
      #for detection model: 
      #define intercept as mean probs:
      int.det <- log( mean.det / ( 1 - mean.det ) )
      mean.det ~ dbeta( 4, 4 )
      
      #random intercept for observer
      for( s in 1:S ){
        eps.det[s] ~ dnorm( 0, pres.det ) T(-7, 7)
      }
      #associated variance of random intercepts:     
      pres.det <- 1/ ( sigma.det * sigma.det )
      #sigma prior specified as a student t half-normal:
      sigma.det ~ dt( 0, 2.5, 7 ) T( 0, )
      
      #priors for detection coefficients:
      #define as a slightly informative prior
      alpha ~ dnorm( 0, 0.1 ) T(-7, 7 )

      #priors for abundance coefficients:
      for( b in 1:B ){
        #define as a slightly informative prior
        beta[ b ] ~ dnorm( 0, 0.1 ) T(-7, 7 )
      }
      
      #prior for abundance model intercept 
      int.lam ~ dgamma( 0.01, 0.01 ) T(-10, 10 )
      
    
    # ecological model of abundance
    for( i in 1:I ){
      N[ i ] ~ dpois( lambda[ i ] )
      log( lambda[ i ] ) <- int.lam + 
                        inprod( beta, XI[ i, ]  )
    
      #observation model
      for( j in 1:J ){
        logit( p[i,j] ) <- int.det + 
                  #random intercept for observer effect
                  eps.det[ obvs[i,j] ] +
                  #fixed effects
                  alpha * time[i,j]
                  
        #for model evaluation:          
        #observed counts
        y_obs[ i, j ] ~ dbin( p[i,j], N[i] )  
        #expected counts
        y_hat[ i, j ] ~ dbin( p[i,j], N[i] )
                
    } #close J
    } #close I
        
    } #model close
     
     ", fill = TRUE )

sink()
################ end of model specification  #####################################     
modelname <- "m1.txt"
#parameters monitored
params <- c(  'int.det' #intercept for detection
              , 'int.lam' #intercept for lamda
              , 'alpha' #detection coefficients
              , 'eps.det' #random intercepts in detection
              , 'sigma.det' #error for random intercept
              , 'beta' #abundance coefficients
              , 'p' #estimate of detection probability
              , 'y_hat' #predicted observations
              , 'N' #estimates of abundance
)

#initial values defined as max counts
Nst <- apply( closeddf[ ,yidx], 1, max, na.rm = TRUE )
#replace 0s with 1
Nst[which(Nst== 0 )] <- 1

#how many ecological predictors that are fixed effects
B <- dim(XI)[2]
#how many detection predictors that are fixed effects
A <- 1
#define initial parameter values
inits <- function(){ list( beta = rnorm( B ),
                           alpha = rnorm( A ),
                           N = Nst ) }

#define data that will go in the model
str( win.data <- list( y_obs = as.matrix( closeddf[ ,yidx] ),
                       #number of sites, surveys, det predictors, and abund preds
                       I = I, J = J, A = A, B = B, S = S,
                       #site level habitat predictors
                       XI = XI,
                       #observation predictors:
                       time = time_sc,
                       obvs = obvs
) )

#call JAGS and summarize posteriors:
m1 <- autojags( win.data, inits = inits, params, modelname, #
                n.chains = 5, n.thin = 10, n.burnin = 0,
                iter.increment = 10000, max.iter = 500000, 
                Rhat.limit = 1.0,
                save.all.iter = FALSE, parallel = TRUE ) 
###### end m1 ########
mr <- m1
whiskerplot( mr, parameters = "eps.det" , zeroline = TRUE)
sort(unique(countinfo$obsv ) )
whiskerplot( mr, parameters = "alpha", zeroline = TRUE )
whiskerplot( mr, parameters = "beta", zeroline = TRUE )
colnames(XIin)
whiskerplot( m1, parameters = "N", zeroline = TRUE )
i_df$id


#################### end of script ############################################      