#######################################################################
#######################################################################
####   This script was created by Dr. Jen Cruz as part of      ########
####            Applied Population Ecology Class               ########
####
#### In this script we will explore some basic theoretical     ########
#### models of population growth. We will then explore indices  #######
###  for estimating abundance. Lastly, we will assess potential  ######
###  impacts of not accounting for imperfect detection.          #####
###
###  This is the first script. As you read it note its format,    #####
### note the heavy editing, the break up into sections for each  #####
### stage of your work flow. Some of these sections will be      #####
### common to all scripts, like the 'initial set up' section     ####
### Note the script ends with an end of script line to prevent   ####
### readers from accidentally missing some lines of code.        ###
### Note the spacing through out. Don't assume your line of code ###
### worked. Check that each line worked. This helps debugping.   ###
### Ensure that you follow these good coding practices as you    ####
### adapt code to make it your own. Do you have other good        ###
### coding practices?                                             ###
### NEVER RUN BIG CHUNKS OF CODE YOU DO NOT UNDERSTAND. NEVER RUN ###
### SECTIONS OF MULTIPLE SCRIPTS RANDOMLY. THAT REMOVES THE ABILITY ###
### FOR YOUR WORKFLOW TO BE REPLICABLE.                           ###
###                                                               ###
###                   Let's get started!                          ###
#####################################################################
##### Setting up your workspace ------------------------------------
# Clean your workspace to reset your R environment. #
rm( list = ls() )
# This removes existing objects you may have created in previous sessions. #
# Those are often saved automatically even after you close R. #

# You also need to set your "working directory".#
# This is the "Path" in your computer where all code, output, #
# and other components of your current Rstudio session will be stored.#
# 
# If you create a new RStudio project, this project's folder will # 
# become your working directory #
# I recommend you do this. If you don't know where you are check:
getwd()

###################################################################
#### Load or create data -----------------------------------------

# In this section you can load or create the data you will need. # 

#### THEORETICAL MODELS OF POPULATION GROWTH #####################
# We start by simulating population growth based on discrete #
# time increments and ignoring population cohorts, such as age cohorts. #
# We model the size of population in the next time period, N[t+1], to #
# depends on the births, B, immigration, I, deaths, D, and emigration, E, #
# as well as the initial population size, N[t] #

# starting with the initial population size at time t
N <- 500
# we can then define B, I, D, E:
B <- 150; D <- 90; I <- 70; E <- 30
#note semicolons allow you to run multiple commands in the same line

# We create our discrete model to calculate N[t+1]
N[2] <- N[1] + (B - D) + (I - E) #recognize this from your readings?
#check what is your resulting population at N[t+1]?
N[2] #since R indexing starts at zero time = 2 is indexed by 2 in this vector

# If the population is closed, we assume that immigration #
#  and emmigration are negligible. we can simplify the #
# population growth model for closed populations as follows: #

# Start with the same initial population size and B, D values
Nclosed <- N #create a new vector so that we don't overwrite the original
Nclosed[2] <- Nclosed[1] + (B-D) #Do you notice what was changed?
# Check the new estimate of N[t+1]:
Nclosed[2] #Note that you check your output all the time!

# We can then calculate the finite rate of increase, R, as the difference 
# between births and deaths per capita:
R <- ( B / Nclosed[1] ) - ( D / Nclosed[1] ) 
#Notice how we turn total to per capita values

# The population at N[t+1] can be calculated as:
Nclosed[1] +  ( Nclosed[1] * R ) 
# Did we get the same result?
# Answer:
#
#
# Now we move to the classic model of Exponential growth, designed for #
# populations capable of continuous growth with overlapping generations #
# We start by calculating the intrinsic rate of increase, sometimes known #
# also as instantaneous growth rate:
r <- log( Nclosed[2] ) - log( Nclosed[1] )
# Note that this is based on integrating and rearranging a differential #
# equation for instantaneous populations growth: dN / dt = rN
# We can then calculate Nclosed[2] as:
Nclosed[1] * exp( r )
# Did we get the same result?
# Answer:
# 
# Note that the population grows if r > 0, is stationary if r = 0, #
# declines if r < 0. When r > 0, and that the change is exponential.

# If we wanted to calculate growth over multiple time periods we can #
# build a loop. Let's assume the intrinsic rate of growth, r, stays #
# constant through time. #

# Start defining how many periods we want to project to:
T <- 50
#new r 
r <- -0.95
# build a loop that moves iteratively through each time period, t
for( t in 1:(T - 1) ){
   Nclosed[ t + 1 ] <- round( Nclosed[ t ] * exp( r ), digits = 0 )
  # we round our population to whole numbers 
   }
# check
Nclosed
# Plot using base plot
plot(x = 1:T, y = Nclosed ) 
points( x = 1:T, y = Nclosed )

# on the x axis we plot the time periods, N[t] on the y axis

# Now modify your r and/or T to see how that changes the shape of your curve #
# Do this here and add the curves to your plot using lines()
# Answer:

#### End of classic models section ######################
###### INDICES OF ABUNDANCE - PERILS OF IGNORING IMPERFECT DETECTION #####
# Now that we explored some classic theoretical models, how do we #
# get estimates of population size, or abundance to plug into our models? #
# Field work! What techniques can we use to sample abundance in the field? #
# Name some here:
# Answer:
# scat sampling, point counts, plot counts, distance sampling,
# camera traps, trapping, nets, pitfall traps,  howling, hooting calls,
# 1m quadrats, radiotracking, sticky traps, road kills, nest searches,
# 
# Add whether you were able to mark individuals with each technique you listed
# Answer:
# 
# Were those marks unique, permanent? What are the implications for your #
# ability to estimate abundance through time? #
# Answer:
#

# If you can mark individuals uniquely then you have a traditional #
# mark-capture-recapture study in your hands! 
# But you need to trap more than once. Repeat surveys during a short enough #
# time window allow you to assume closure so that you can focus on  #
# how many individuals you are not detecting, as a result of your sampling #
# protocol. This information allows you to estimate both abundance, N, and #
# detection, p. 
# 
# what happens if you cannot tell individuals apart but you can still count them? #
# Historically this gave you an index of abundance, which is the product of #
# detection, p, and abundance, N= Np
# BUT recent statistical advances now allow you to use repeated counts #
# meaning over multiple repeat surveys, to estimate abundance and detection #
# separately using N-mixture models (Royle 2004)

# What if you can only detect signs that the species is present but not #
# of how many individuals you have? If you have repeat surveys you can #
# estimate both occupancy, phi, and detection, p, using occupancy models #
# McKenzie et al. 2002.

# Notice the common parameter across this modeling techniques is p, but #
# to estimate it we need repeat surveys. These 3 model types are the main focus #
# of this class. #

# But let's go back in time before these models were readily available #
# to a time when indices were common. #
# 
# Catch per unit effort, is the number of individuals / effort
# 
# Let's simulate a scenario where we put out 10 nets to catch some fish
# define our number of nets
nets <- 10
# now randomly drawn some fish into each net between the 50 to 300 range:
catch <- sample( 50:300, size = nets, replace = TRUE )
# what was our catch?
catch
# what is our catch per unit effort? #
# Answer:
#
# We calculate it as total catch / number of nets put out:#
sum( catch ) / nets

# We place 10 more nets out but this time our nets end up on the shallows #
# We don't realize it but we are unable to catch as many fish in these nets # 
added.catch <- sample( 10:100, size = nets, replace = TRUE )
# catch per unit effort for this second lot of nets
sum( added.catch ) / ( nets ) 

# If we didn't know these two sets were placed in different locations #
# influencing how many fish we were able to catch we would assume all nets #
# were the same, we would have expected uniform catch throughout the lake:
( sum( catch ) + sum( added.catch ) ) / ( nets + nets )

# So being aware of things that influence our catchability? or in other words #
# our sampling, or ability to detect individuals is important in coming up #
# with reasonable estimates of abundance. 

# Lincoln-Petersen estimate
# Now let's go out and use pitfall traps to catch some sagebrush lizards #
# We marked them using nail polish.
# the true number of lizards is 100
N <- 500
#Let's set seed so that we get the same values
set.seed( 2021 )
# We can trap anywhere from none to all:
n1 <- sample( 0:N, 1 ) 
# How many did we get?
n1
# We go out again on day 2 to check our pitfall traps again #
# how many do we catch this time? 
n2 <- sample( 0:N, 1 )
# How many did we get?
n2
# how many of these were marked at time 1? 
# we can't get more than were marked during time 1 #
# or that were trapped at time 2, so we draw a sample smaller than #
# the smallest of both
m2 <- sample( 1:min( n1, n2 ), 1 )
# check 
m2

# If we assume that capture probability is the same for all animals, #
# that marking doesn't affect their detection, no marks are lost, 
# and each sample is a random sample of the population, we can use 
# the Lincoln-Petersen index:
N_est <- ( n1 * n2 )  / m2
# What is our estimate of abundance of lizards?
N_est

# How can we repeat this process a few times to see how much our 
# N_est will vary when our initial captures, n1, recaptures, m2, 
# and total captures at time 2, n2, vary?
# remember our true abundance:
N
# We can repeat this exercise multiple times by putting all #
# the steps above together in a loop, converting our object to vectors #

# Start by defining how many times we repeat the study. #
# Say that you assumed each repeat was a different site that also #
# had abundance = N #
I <- 30
# Then add values to each vector you created above:
for( i in 2:I){
  n1[i] <- sample( 0:N, 1 )
  n2[i] <- sample( 0:N, 1 )
  m2[i] <- sample( 1:min( n1[i], n2[i] ), 1 )
  N_est[i] <- ( n1[i] * n2[i] )  / m2[i]
}
# What were our estimates:
round( N_est, 0 )
# Note some of our estimates are way off our true N
# let's explore why? 
round( N_est, 0)
round( m2/N_est, 2)
round( n1/N_est, 2)
round( n2/N_est, 2)
# What reasons can you see that result in largely imprecise #
# estimates?
# Answer:
#
###### Schnabel index:
# This index builds up on some of the concepts of keeping track of #
# how many individuals are captured each survey, C[t], marking them, #
# and then keeping track of how many marked and available for recapture, #
# M[t+1], and how many are actually recaptured in the following survey, #
# R[t+1].
# Abundance is estimated as:
# N_est = summation( C[t] * M[t] ) / summation( R[t]), over t surveys

#Let's simulate some data first:
# Let's define what our true population abundance is:
N <- 5000
# our total number of surveys
T <- 30
# let's assume that capture probability is fairly constant among surveys #
# because we are using the same protocol including standardizing weather #
#  conditions, observer, and number of traps #

# We will create a loop to estimate abundance 10 times, we start with #
# drawing capture probabilities in the range between 0.4 and 0.6
P <- runif( n = 10, min = 0.4, max = 0.6 )
P
# Now we simulate data for each capture probability, printing out #
# our estimate of abundance, N_est, for each:
# Start by simulating the number of animals captured in trapping surveys #
# assuming capture probability is p
C <- N * sample( P, T, replace = TRUE ) 
print( C )
# create vector for newly marked individuals, assign newM[1] = C[1]:
newM <- c( C[1], rep(NA, T-1) )
print( newM )
#Set vector of marked individuals, available for capture, with M[1] = 0:
M <- c( 0, rep(NA, T-1) )
# Loop iteratively over following surveys:
for( t in 2:T ){
  # newly marked individuals in later surveys are drawn randomly 
  # to be less than half of those captured, unless
  # M > N at which point there are no new individuals to mark
  newM[t] <- ifelse( M[t-1] >= N, 0, sample(1:(C[t]/2),1 ) ) 
  # when newM makes M > N correct it to make it N:
  newM[t] <- ifelse( (newM[t] + M[t-1]) >= N, (N - M[t-1]), newM[t] )
  # update the number of marked individuals available for capture:
  M[t] <- sum( newM[1:t] )
}
# recaptures are the individuals that were not newly marked each survey:
R <- C - newM 
# We have simulated data for capture probability p.
# Now use simulated values to estimate abundance using the Schanabel index:
N_est <- sum( C * M ) / sum( R )
# print our estimate as we run through the loop
print( N_est )


# Repeat the exercise above in a situation where you have low detection #
# so P between 0.1 and 0.2
# What happens to the estimates of abudance?
# Answer:
# What if you increase the number of sampling sessions to T = 40?
# Answer:
#
# What does this tell you about the influence of low detection on this index?
# Answer:

#############end of indices section #########################
#
########################################################################


###################   END OF SCRIPT  ################################