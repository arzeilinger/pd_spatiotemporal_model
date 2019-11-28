#### Function to generate time intervals from infection times


getTimeIntervals <- function(Inf_times = Inf_times, tcuts = tcuts, Tmax = Tmax){
  ## Required inputs:
  ## Inf_times = vector of infection times
  ## tcuts = vector of time points where disease survey observations occurred
  ## Generate matrix of zeros 
  time_intervals <- matrix(Tmax, ncol = 2, nrow = length(Inf_times))
  ## t is the upper limit of the interval (i.e., [t-1, t]), so start at index 2
  for(i in 1:length(Inf_times)){
    if(Inf_times[i] < Tmax){
      for(t in 2:length(tcuts)){
        if(Inf_times[i] > tcuts[t-1] & Inf_times[i] <= tcuts[t]){
          ## Set lower observation: last time point plant was observed asymptomatic
          time_intervals[i,1] <- tcuts[t-1]
          ## Set upper observation: first time point plant was observed symptomatic
          time_intervals[i,2] <- tcuts[t]
        }
      }
    }
  } 
  return(time_intervals)
}
