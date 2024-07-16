#' Read and analyse spikes in a single ACQ4 experiment
#'
#' @param base_folder
#' @param query_folder
#' @param file_extension
#' @param threshold
#'
#' @return
#' @export
#' @importFrom mmand gaussianSmooth
#' @importFrom gsignal sgolayfilt gsignal findpeaks
#'
#' @examples
#' \dontrun{
#' }
spikes <- function(base_folder, query_folder, file_extension, threshold) {
  # Get a list of all files in the specified folder
  files <- list.dirs(file.path(base_folder, query_folder))
  files <- sort(files)
  file_list <- files[-1]

  # Initialize a list to store the spike times for each file
  spike_times <- list()

  # Loop through each file and extract the spike times
  for (file in file_list) {
    # Read in the data
    h5_file <- h5read(paste(file, "/Clamp2.ma", sep=""), name="data")

    # Use the second derivative of your data to detect spikes
    sg <- sgolayfilt(gaussianSmooth(h5_file[-(1:300),2]+0.55, 50), m = 2)

    # Find all the local maxima in your data above the threshold
    peaks <- findpeaks(data=gaussianSmooth(sg, 20), MinPeakHeight = threshold, MinPeakDistance = 150, DoubleSided = TRUE)$loc

    # Store the spike times for this file in the list
    spike_times[[sub("^.*?(\\d{4}\\.\\d{2}\\.\\d{2}_\\d+)(.+?)", "\\1\\2\\", file)]] <- peaks
  }

  return(spike_times)
}

#' Read and analyse spikes in all experiments of a single type for a single cell
#'
#' @param base_folder
#' @param experiment_type
#' @param file_extension
#' @param threshold
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' }
spike_times_folder <- function(base_folder,experiment_type,file_extension, threshold) {
  spike_times <- list()
  # Get a list of all files in the specified folder
  files <- list.dirs(file.path(base_folder))
  myfiles=files[grep(paste0(experiment_type, "_\\d{3}/"),list.dirs(file.path(base_folder)))]
  for (file in myfiles) {
    # Read in the data
    h5_file <- h5read(paste(file, "/Clamp2.ma", sep=""), name="data")

    # Use the second derivative of your data to detect spikes
    sg<- sgolayfilt(gaussianSmooth(h5_file[-(1:300),2]+0.55, 50), m = 2)

    # Find all the local maxima in your data above the threshold
    peaks <- findpeaks(data=gaussianSmooth(sg, 50), MinPeakHeight = threshold, MinPeakDistance = 150, DoubleSided = TRUE)$loc

    # Store the spike times for this file in the list
    spike_times[[file]][[sub("^.*?(\\d{4}\\.\\d{2}\\.\\d{2}_\\d+)(.+?)", "\\1\\2\\", file)]] <- peaks
  }
  return(spike_times)
}


#' Read data from an ACQ4 HDF5 file and turn it into an R time series
#'
#' @param files Input acq4 HDF5 files (often end in `.ma` or `.h5`)
#' @param name the name of the HDF5 folder containing the data to read
#' @param bit64conversion What to do with 64 bit integers (this normally comes
#'   up with unsigned 32 bit integers for which R has no native support but can
#'   be encoded without any loss of precision in R's `numeric` (8 byte double)
#'   type). See [h5read()] for details.
#' @param channel Which data channel to import
#' @param ... Additional arguments passed to [h5read()]
#'
#' @return A `ts` or `mts` object
#' @importFrom stats deltat ts
#' @importFrom rhdf5 h5read
#' @export
#'
#' @examples
acq4h52ts <- function(files, name='data', bit64conversion='double', channel=2, ...) {
  tracelist=lapply(files, function(x) h5read(x, name=name, bit64conversion=bit64conversion, ...)[,channel])
  times=h5read(files[1], name='info/1/values', bit64conversion=bit64conversion)
  deltat=diff(times)[1]
  # if(length(deltat)!=1)
  #   stop("Unable to extract unique deltat for time series")
  if(length(tracelist)>1)
    ts(do.call(cbind, tracelist), deltat = deltat, start=0)
  else
    ts(tracelist[[1]], deltat = deltat, start=0)
}

#' Downsample a signal by an integer factor
#'
#' @inheritParams gsignal::downsample
#' @return
#' @export
#' @importFrom stats is.ts is.mts deltat start
#' @importFrom gsignal downsample
#'
#' @examples
#' x=1:10
#' x
#' downsample(x, 2)
#' x3=matrix(1:30, ncol=3, dimnames=list(NULL, c("x", "y", "z")))
#' x3
#' downsample(x3, 2)
downsample <- function(x, n, phase = 0) UseMethod('downsample')

#' @rdname downsample
#' @export
downsample.default <- function(x, n, phase = 0) {
  gsignal::downsample(x, n=n, phase = phase)
}

#' @rdname downsample
#' @export
#' @examples
#' xts=stats::ts(1:10, start=0, deltat=0.1)
#' xts
#' downsample(xts, 2)
#'
#' x3ts=stats::ts(matrix(1:30, ncol=3, dimnames=list(NULL, c("x", "y", "z"))), start=0, deltat=0.1)
#' x3ts
#' downsample(x3ts, 2)
downsample.ts <- function(x, n, phase = 0) {
  stopifnot(is.ts(x))
  dx=if(is.mts(x)) as.matrix(x) else as.vector(x)
  ts(gsignal::downsample(dx, n=n, phase = phase), deltat = deltat(x)*n, start = start(x))
}


#' Generate plots for all the results from a single cell
#'
#' @param path
#' @param stf_list
#' @return
#'
#' @export
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics rug
#'
#' @examples
plot_cell <- function(path, stf_list) {
  # Check if stf_list is a character string and convert it to a list if necessary
  if (is.character(stf_list)) {
    stf_list <- get(stf_list, envir = .GlobalEnv)
  }
  stf_list=stf_list
  # Get all subdirectories
  subdirs <- list.dirs(path, recursive = TRUE, full.names = TRUE)

  for (subdir in subdirs) {
    # List all .ma files in the current subdirectory
    ma_files <- list.files(subdir, pattern = "\\.ma$", full.names = TRUE)

    for (ma_file in ma_files) {
      cat("Processing file:", ma_file, "\n")

      # Extract components of the path
      path_components <- unlist(strsplit(ma_file, "/"))
      cat("Path components:\n", path_components, "\n")

      # Identify the correct level for subdirectory
      base_path <- paste(path_components[1],path_components[2],path_components[3],path_components[4],path_components[5],path_components[6],path_components[7], path_components[8], sep = "/")
      experiment <- path_components[9]
      trial_dir <- path_components[10]
      trial <- path_components[11]  # Trial directory (e.g., 000, 001, etc.)

      # Construct the key correctly without the last component (Clamp2)
      # message("trial:", trial)
      if (length(trial)<1 || !isTRUE(trial=="Clamp2.ma")) {
        message("failed for base_path:", base_path)
        next
      }
      base_stf_key=paste(base_path, experiment, trial_dir, sep = "/")
      stf_key <- paste(base_path, experiment, trial_dir, trial, sep = "/")
      cat("Constructed key:", stf_key, "\n")

      # Check if the constructed key exists in the stf_list
      if (base_stf_key %in% names(stf_list)) {
        cat("Key found in stf_list.\n")
        # Load the data
        op=options(warn = 2)

        # ma_data <- try(h5read(stf_key, name = "data"))
        ma_data <- try(acq4h52ts(stf_key, name = "data"))
        if(inherits(ma_data, 'try-error')){
          message("There was a problem reading file: ", base_stf_key)
          next
        }
        options(op)

        # Get the indices for each stf signal
        stf_indices <- stf_list[[base_stf_key]]
        stf_times=unlist(stf_indices, use.names = F)*deltat(ma_data)

        ds=downsample_ts(ma_data*1e2, 50)
        # ds=downsample_ts(ma_data[,2]*1e2, 50)

        # Plot the data
        tryCatch({
          pdf(paste(base_stf_key,".pdf",sep=""), width = 7, height = 7)
          plot(ds, ylab='Vm /mV', type = "l")#check fix was phys_data
          rug(stf_times)
          },
          finally = dev.off())

        #plot(downsample(ma_data[,2],50),type="l")
       # rug((unname(unlist(stf_list[[base_stf_key]])/50)))
       # legend("topright", legend = paste("Signal", stf_indices), col = 1:length(stf_indices), lwd = 2)
      } else {
        cat("Key not found in stf_list.\n")
      }
    }
  }
}



deriv <- function(x, y) diff(y) / diff(x)
middle_pts <- function(x) x[-1] - diff(x) / 2
# second_d <- deriv(middle_pts(mid), deriv(mid, strike))

