#### Function to generate time intervals from infection times


getTimeIntervals <- function(Inf_times = Inf_times, tcuts = tcuts){
  ## Required inputs:
  ## Inf_times = vector of infection times
  ## tcuts = vector of time points where disease survey observations occurred
  ## Generate matrix of zeros 
  time_intervals <- matrix(0, ncol = 2, nrow = length(Inf_times))
  for(t in 1:length(tcuts)){
    if(tcuts[t] < max(tcuts)){
      for(i in 1:length(Inf_times)){
        if(Inf_times[i] > tcuts[t] && Inf_times[i] <= tcuts[t+1]){
          ## Set lower observation: last time point plant was observed asymptomatic
          time_intervals[i,1] <- tcuts[t]
          ## Set upper observation: first time point plant was observed symptomatic
          time_intervals[i,2] <- tcuts[t+1]
        }
      }
    }
  } 
  return(time_intervals)
}
