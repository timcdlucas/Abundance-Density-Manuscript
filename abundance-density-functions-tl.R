
# Function like findDisDistr in MetapopEpi
#   However finds the distr a certain time since invasion
findDisDistrAtTime <- function(pop, final = 1, invadeT, time = 30){
    assertthat::is.count(final)

    invadeTime <- cumsum(pop$sampleWaiting)[invadeT / pop$parameters["sample"]]

    # check that the sim was long enough
    if(sum(pop$sampleWaiting) < invadeTime + time) stop('Did not run sim for long enough')
    
    # Find which sample index is "time" after invadeTime
    indextoSample <- which(cumsum(pop$sampleWaiting) > (invadeTime + time))[1]
    
    finalI <- apply(pop$sample[, , (indextoSample - final):indextoSample], 1, mean)
    finalDis <- sapply(1:pop$parameters["nPathogens"], function(x) finalI[pop$whichClasses[, 
        x]] %>% sum)

    finalDis
}





  #################################
  # Density sim definitions       #
  #################################


# Define our simulation function.
fullSim1 <- function(x){

  # Set seed (this is set within each parallel simulation to prevent reusing random numbers).
  simSeed <- paste0(seed, x)
  set.seed(simSeed)

  # Make the population.
  p <- makePop(space = side[x], 
               transmission = tran[x], 
               meanColonySize = colonySize[x], 
               nColonies = 20, 
               model = 'SIR', 
               events = nEvent, 
               nPathogens = 2, 
               recovery = 0.5,  
               sample = sample, 
               dispersal = 0.01, 
               birth = 0.05, 
               death = 0.05,
               crossImmunity = 0.1, 
               infectDeath = 0,
               maxDistance = 100)

  p$I[1, , 1] <- colonySize[x] - 20

  # Seed endemic pathogen.
  p$I[3, , 1] <- 20
  
  # Recalculate rates of each event after seeding.
  p <- transRates(p, 1)

  # Burn in simulation
  p <- runSim(p, end = invadeT)

  # Seed invading pathogen.
  p$I[2, 1, (invadeT + 1) %% sample] <- 5
  
  # Recalculate rates of each event after seeding.
  p <- transRates(p, (invadeT + 1) %% sample)

  # Continue simulation
  p <- runSim(p, start = invadeT + 1, end = 'end')

  # Was the invasion succesful?
  invasion <- findDisDistrAtTime(p, final = 10, invadeT, time = TotalTime)[1] > 0

  # Save summary stats
  d <- data.frame(transmission = NA)

  d$transmission <- p$parameters['transmission']
  d$dispersal <- p$parameters['dispersal']
  d$nExtantDis <- sum(findDisDistrAtTime(p, final = 2, invadeT, time = TotalTime) > 0)
  d$nPathogens <- p$parameters['nPathogens']
  d$meanK <- sum(p$adjacency != 0 )/p$parameters['nColonies']
  d$maxDistance <- p$parameters['maxDistance']
  d$nEvents <- p$parameters['events']
  d$colonySize <- p$parameters['meanColonySize']
  d$colonyNumber <- p$parameters['nColonies']
  d$pop <- p$parameters['meanColonySize'] * p$parameters['nColonies']
  d$area <- p$parameters['space']^2
  d$dens <- d$pop / d$area
  d$invadeTime <- cumsum(p$sampleWaiting)[invadeT / p$parameters["sample"]]
  d$longenough <- sum(p$sampleWaiting) > d$invadeTime + TotalTime
  #d$path2 <- sum(p$sample[c(2, 4), , dim(p$sample)[3]])



  # Time until extinction
  invadePath <- colSums(p$sample[2,  , (2 + invadeT / sample):(which(cumsum(p$sampleWaiting) > TotalTime)[1])]) + 
                  colSums(p$sample[4,  , (2 + invadeT / sample):(which(cumsum(p$sampleWaiting) > TotalTime)[1])])

  d$invadeTime <- cumsum(p$sampleWaiting)[(2 + invadeT / sample)]
  d$extinctionTime <- cumsum(p$sampleWaiting)[min(which(invadePath == 0)) + (2 + invadeT / sample)]
  d$totalTime <- sum(p$sampleWaiting)
  d$survivalTime <- d$extinctionTime - cumsum(p$sampleWaiting)[(2 + invadeT / sample)]


  write(paste0("finished ", x, ". Invasion: ", invasion ), 'runlog.R', append = TRUE)
  message(paste0("finished ", x, ". Invasion: ", invasion ))
  message(paste('trans:', d$transmission, 'nColonies:', d$colonyNumber, 
    'colonySize:', d$colonySize, 'dens:', d$dens, 'pop:', d$pop))

  message(paste(names(d), collapse = ', '))
  message(paste(d, collapse = ', '))

  
  if(saveData){ 
    file <- paste0('Data/DensSim_', formatC(x, width = 4, flag = '0'), '.RData')
    save(p, file = file)
  }

  #rm(p)

  return(d)

}






  #################################
  # Density2 sim definitions      #
  #################################

# Define our simulation function.
fullSim2 <- function(x){

  # Set seed (this is set within each parallel simulation to prevent reusing random numbers).
  simSeed <- paste0(seed, x)
  set.seed(simSeed)

  # Make the population.
  p <- makePop(space = side[x], 
               transmission = tran[x], 
               meanColonySize = colonySize, 
               nColonies = colonyNumber[x], 
               model = 'SIR', 
               events = nEvent, 
               nPathogens = 2, 
               recovery = 0.5,  
               sample = sample, 
               dispersal = 0.01, 
               birth = 0.05, 
               death = 0.05,
               crossImmunity = 0.1, 
               infectDeath = 0,
               maxDistance = 100)

  p$I[1, , 1] <- colonySize - 20

  # Seed endemic pathogen.
  p$I[3, , 1] <- 20
  
  # Recalculate rates of each event after seeding.
  p <- transRates(p, 1)

  # Burn in simulation
  p <- runSim(p, end = invadeT)

  # Seed invading pathogen.
  p$I[2, 1, (invadeT + 1) %% sample] <- 5
  
  # Recalculate rates of each event after seeding.
  p <- transRates(p, (invadeT + 1) %% sample)

  # Continue simulation
  p <- runSim(p, start = invadeT + 1, end = 'end')

  # Was the invasion succesful?
  invasion <- findDisDistrAtTime(p, final = 10, invadeT, time = TotalTime)[1] > 0

  # Save summary stats
  d <- data.frame(transmission = NA)

  d$transmission <- p$parameters['transmission']
  d$dispersal <- p$parameters['dispersal']
  d$nExtantDis <- sum(findDisDistrAtTime(p, final = 2, invadeT, time = TotalTime) > 0)
  d$nPathogens <- p$parameters['nPathogens']
  d$meanK <- sum(p$adjacency != 0 )/p$parameters['nColonies']
  d$maxDistance <- p$parameters['maxDistance']
  d$nEvents <- p$parameters['events']
  d$colonySize <- p$parameters['meanColonySize']
  d$colonyNumber <- p$parameters['nColonies']
  d$pop <- p$parameters['meanColonySize'] * p$parameters['nColonies']
  d$area <- p$parameters['space']^2
  d$dens <- d$pop / d$area
  d$invadeTime <- cumsum(p$sampleWaiting)[invadeT / p$parameters["sample"]]
  d$longenough <- sum(p$sampleWaiting) > d$invadeTime + TotalTime
  #d$path2 <- sum(p$sample[c(2, 4), , dim(p$sample)[3]])



  # Time until extinction
  invadePath <- colSums(p$sample[2,  , (2 + invadeT / sample):(which(cumsum(p$sampleWaiting) > TotalTime)[1])]) + 
                  colSums(p$sample[4,  , (2 + invadeT / sample):(which(cumsum(p$sampleWaiting) > TotalTime)[1])])

  d$invadeTime <- cumsum(p$sampleWaiting)[(2 + invadeT / sample)]
  d$extinctionTime <- cumsum(p$sampleWaiting)[min(which(invadePath == 0)) + (2 + invadeT / sample)]
  d$totalTime <- sum(p$sampleWaiting)
  d$survivalTime <- d$extinctionTime - cumsum(p$sampleWaiting)[(2 + invadeT / sample)]


  write(paste0("finished ", x, ". Invasion: ", invasion ), 'runlog.R', append = TRUE)
  message(paste0("finished ", x, ". Invasion: ", invasion ))
  message(paste('trans:', d$transmission, 'nColonies:', d$colonyNumber, 
    'colonySize:', d$colonySize, 'dens:', d$dens, 'pop:', d$pop))

  message(paste(names(d), collapse = ', '))
  message(paste(d, collapse = ', '))

  
  if(saveData){ 
    file <- paste0('Data/DensSim_', formatC(x, width = 4, flag = '0'), '.RData')
    save(p, file = file)
  }

  #rm(p)

  return(d)

}







  #################################
  # Population sim definitions    #
  #################################


# Define our simulation function.
fullSim3 <- function(x){

  # Set seed (this is set within each parallel simulation to prevent reusing random numbers).
  simSeed <- paste0(seed, x)
  set.seed(simSeed)

  # Make the population.
  p <- makePop(space = side[x], 
               transmission = tran[x], 
               meanColonySize = colonySize, 
               nColonies = 20, 
               model = 'SIR', 
               events = nEvent, 
               nPathogens = 2, 
               recovery = 0.5,  
               sample = sample, 
               dispersal = 0.01, 
               birth = 0.05, 
               death = 0.05,
               crossImmunity = 0.1, 
               infectDeath = 0,
               maxDistance = 100)

  p$I[1, , 1] <- colonySize - 20

  # Seed endemic pathogen.
  p$I[3, , 1] <- 20
  
  # Recalculate rates of each event after seeding.
  p <- transRates(p, 1)

  # Burn in simulation
  p <- runSim(p, end = invadeT)

  # Seed invading pathogen.
  p$I[2, 1, (invadeT + 1) %% sample] <- 5
  
  # Recalculate rates of each event after seeding.
  p <- transRates(p, (invadeT + 1) %% sample)

  # Continue simulation
  p <- runSim(p, start = invadeT + 1, end = 'end')

  # Was the invasion succesful?
  invasion <- findDisDistrAtTime(p, final = 10, invadeT, time = TotalTime)[1] > 0

  # Save summary stats
  d <- data.frame(transmission = NA)

  d$transmission <- p$parameters['transmission']
  d$dispersal <- p$parameters['dispersal']
  d$nExtantDis <- sum(findDisDistrAtTime(p, final = 2, invadeT, time = TotalTime) > 0)
  d$nPathogens <- p$parameters['nPathogens']
  d$meanK <- sum(p$adjacency != 0 )/p$parameters['nColonies']
  d$maxDistance <- p$parameters['maxDistance']
  d$nEvents <- p$parameters['events']
  d$colonySize <- p$parameters['meanColonySize']
  d$colonyNumber <- p$parameters['nColonies']
  d$pop <- p$parameters['meanColonySize'] * p$parameters['nColonies']
  d$area <- p$parameters['space']^2
  d$dens <- d$pop / d$area
  d$invadeTime <- cumsum(p$sampleWaiting)[invadeT / p$parameters["sample"]]
  d$longenough <- sum(p$sampleWaiting) > d$invadeTime + TotalTime
  #d$path2 <- sum(p$sample[c(2, 4), , dim(p$sample)[3]])



  # Time until extinction
  invadePath <- colSums(p$sample[2,  , (2 + invadeT / sample):(which(cumsum(p$sampleWaiting) > TotalTime)[1])]) + 
                  colSums(p$sample[4,  , (2 + invadeT / sample):(which(cumsum(p$sampleWaiting) > TotalTime)[1])])

  d$invadeTime <- cumsum(p$sampleWaiting)[(2 + invadeT / sample)]
  d$extinctionTime <- cumsum(p$sampleWaiting)[min(which(invadePath == 0)) + (2 + invadeT / sample)]
  d$totalTime <- sum(p$sampleWaiting)
  d$survivalTime <- d$extinctionTime - cumsum(p$sampleWaiting)[(2 + invadeT / sample)]

  write(paste0("finished ", x, ". Invasion: ", invasion ), 'runlog.R', append = TRUE)
  message(paste0("finished ", x, ". Invasion: ", invasion ))
  message(paste('trans:', d$transmission, 'nColonies:', d$colonyNumber, 
    'colonySize:', d$colonySize, 'dens:', d$dens, 'pop:', d$pop))

  message(paste(names(d), collapse = ', '))
  message(paste(d, collapse = ', '))

  
  if(saveData){ 
    file <- paste0('Data/PopSim_', formatC(x, width = 4, flag = '0'), '.RData')
    save(p, file = file)
  }

  #rm(p)

  return(d)

}
