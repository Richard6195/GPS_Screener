
#' GPS_Screener
#'
#' @description This function screens GPS data for anomalies using a combination of GPS-derived metrics and, if provided, acceleration data (e.g., VeDBA or equivalent) within isolation forest models for unsupervised anomaly detection. It supplements isolation forest results with biologically realistic thresholds for additional validation. If plot = TRUE, a dynamic leaflet map visualizing the anomaly detection results is rendered.
#'
#' Key Features:
#' - Innovative Use of GPS Metrics: Identifies anomalies using movement metrics derived from three consecutive pre-processed fixes, including circular distance, outbound and inbound speeds, and vertex angle.
#' - Adaptive Handling of Irregular Fixes: Processes GPS bursts and 1 Hz periods with options for standardizing fix intervals and handling irregular time gaps.
#' - Acceleration Data Integration: Supports optional acceleration data (e.g., VeDBA) to enhance anomaly detection by leveraging activity metrics alongside GPS-derived measures.
#' - Sophisticated Isolation Forest Approach: Employs a multidimensional isolation forest model that detects anomalies through associations between variables, enhanced by post-model threshold-based logic for improved accuracy.
#' 
#' @param GPS_TS POSIXct vector of timestamps for GPS data.
#' @param GPS_longitude Numeric vector of longitude values (in decimal degrees).
#' @param GPS_latitude Numeric vector of latitude values (in decimal degrees).
#' @param ACC_TS OPTIONAL: POSIXct vector of timestamps for VeDBA (or any activity metric) to enahnce anamaly detection. 
#' @param VeDBA OPTIONAL: Numeric vector of VeDBA (or any activity metric), the same length as 'ACC_TS' if supplied. This does not have to match the length of GPS inputs and function deals with merging with GPS data (timestamps do not have to match with GPS data). Mean values are computed per GPS fix (backward in time from fix index) at times when sufficient data is present. This is used to inform isolation forest anomaly detection and additional parameter checks.  
#' @param drop_out Numeric value defining the maximum drop-out time (in seconds) for GPS fixes. Defines when a new "group" begins (default = 600 s). The functionality is obsolete within the function, but output data frame contains this and its useful for identifying large temporal gaps in data sets. 
#' @param burst_method Character specifying how to handle GPS bursts, with options: "median", "mean", "last", or "none" (default = "last"). This filters fixes to a single value per burst before anomaly detection. Default = NULL.
#' @param burst_len Numeric value defining the maximum length of a typical GPS burst (in seconds). Default = 1 (assumes no GPS bursts).
#' @param standardise_time_interval Optional Numeric value defining the stepping range of retaining fixes (in seconds), useful for standardizing fix intervals of 1 Hz GPS periods to typical recording schedule (default = NULL; not applied).
#' @param standardise_universal Logical value (TRUE/FALSE) indicating whether 'standardise_time_interval' is to be applied across all GPS data or just the 1 Hz periods.
#' @param IF_conf Numeric value between 0 and 1 specifying the confidence quantile for isolation forest anomaly detection (default = 0.99).
#' @param iso_sample_size Numeric value specifying the sample size for isolation forest anomaly detection (default = 256). Though I recommmed 500-1000.
#' @param GPS_height Numeric vector specifying the height of GPS fixes (e.g., 'height.above.ellipsoid' in Eobs' GPS files), or NULL if not available. If provided, it must be the same length as GPS_longitude, GPS_latitude, and GPS_TS. This is as an additional variable within the isolation forest anomaly detection.
#' @param start_timestamp (Default: NULL): The earliest timestamp to include in the dataset (in YYYY-MM-DD HH:MM:SS) used for subsetting both GPS and (if supplied) VeDBA data sets.
#' @param end_timestamp (Default: NULL): The latest timestamp to include in the dataset (in YYYY-MM-DD HH:MM:SS) used for subsetting both GPS and (if supplied) VeDBA data sets.
#' @param max_speed (Default: 5): This is a single numeric value specifying the maximum biologically realistic speed the animal can travel between fixes to act as an additional safeguard after unsupervised anomaly detection. 
#' @param GPS_accuracy (Default: 25): This parameter is used in further anomaly detection logic after isolation forests. It specifies approx error radius of fixes under stationary conditions. If unknown, best to leave as default which represents a reasonable estimate (based on empirical data).
#' @param plot Logical value (TRUE/FALSE) indicating whether to generate diagnostic plots of threshold and isolation forest results (default = TRUE).
#'
#' @return Depending on the value of `plot`, either:
#' - A list containing:
#'   - **Thresholds**: histograms depicting various GPS/ACC metrics' distributions.
#'   - **Plot**: Interactive GPS plot using leaflet showing results of anomaly detection.
#'   - **df**: Data frame of results including all calculated metrics and anomaly assessments.
#' 
#' @details
#' **Components of the Function**:
#' 1. **User Defined Movement Thresholds**:
#'   -  Identifies erroneous "spikes" using thresholds for both outgoing/incoming speed, circular distance and turning angles. Optinally includes acceleration and GPS height metrics.
#'
#' 2. **Isolation Forest Models**:
#'   - Implements isolation forest to detect potential anomalies based on various GPS and optiionally ACC metrics. This unsupervised approach can isolate data points using random partitioning.
#'   - Multi-dimensional isolation forests are run, and additional threshold checks are used after the run to make anomaly detection more accurate.

###################################################################################################################################################################################################################################################################################################

#################### START OF FUNCTION ####################

GPS_Screener = function(GPS_TS,  GPS_longitude,  GPS_latitude,  ACC_TS = NULL, VeDBA = NULL, drop_out = 3600, burst_method = "last", standardise_time_interval = NULL, standardise_universal = FALSE, burst_len = 1, IF_conf = 0.99, iso_sample_size = 256, GPS_height = NULL, start_timestamp = NULL,  end_timestamp = NULL, max_speed = 5, GPS_accuracy = 25, plot = TRUE) {
  
  ################################################################################################################################################  
  
  ###Required packages###
  
  required_packages <- c('zoo', 'dplyr', 'tidyr', 'data.table', 'Rcpp', 'ggplot2', 'plotly', 'assertthat', 'leaflet', 'isotree')
  lapply(required_packages, function(pkg) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE, type = "source")
      suppressMessages(library(pkg, character.only = TRUE))
    }
  })
  
  #Check that required packages are installed on the system
  areinstaled = data.frame(installed.packages())
  
  if(all(c('zoo', "dplyr","tidyr", "data.table", "ggplot2", 'leaflet', 'isotree', "assertthat") %in% areinstaled$Package) == FALSE){
    required_packages=c("zoo", "dplyr","tidyr", "data.table", "ggplot2", 'leaflet', 'isotree', "assertthat")
    missing_packages=c("zoo", "dplyr","tidyr", "data.table", "ggplot2", 'leaflet', 'isotree', "assertthat") %in% areinstaled$Package
    stop(paste("The following packages are not installed:", required_packages[which(missing_packages == FALSE)], sep = " "))
  }
  
  ###Input argument checking ###
  
  options(digits.secs = 3) #Specify the number of decimal places of the fractional seconds to show if relevant   
  is.POSIXct = function(x) inherits(x, "POSIXct") #Function to check variable is of POSIXct type
  
  
  validate_inputs <- function(GPS_TS, GPS_longitude, GPS_latitude, ACC_TS = NULL, VeDBA = NULL, 
                              drop_out, burst_method, burst_len, standardise_time_interval = NULL, 
                              standardise_universal, IF_conf, iso_sample_size, GPS_height = NULL, 
                              start_timestamp = NULL, end_timestamp = NULL, max_speed, GPS_accuracy, plot) {
    
    # Ensure GPS_TS is of type POSIXct
    assert_that(is.POSIXct(GPS_TS), msg = "GPS_TS must be of type POSIXct")
    assert_that(all(!is.na(GPS_TS)), msg = "GPS_TS cannot contain missing data")
    
    # Ensure GPS_TS does not contain duplicate timestamps
    assert_that(length(unique(GPS_TS)) == length(GPS_TS), msg = "GPS_TS must not contain duplicate timestamps")
    
    # Ensure GPS_longitude is between -180 and 180
    assert_that(min(GPS_longitude, na.rm = TRUE) >= -180 && max(GPS_longitude, na.rm = TRUE) <= 180, 
                msg = "GPS_longitude must be between -180 and 180 degrees")
    
    # Ensure GPS_latitude is between -90 and 90
    assert_that(min(GPS_latitude, na.rm = TRUE) >= -90 && max(GPS_latitude, na.rm = TRUE) <= 90, 
                msg = "GPS_latitude must be between -90 and 90 degrees")
    
    # Ensure GPS_longitude, GPS_latitude, and GPS_TS have the same length
    assert_that(length(GPS_longitude) == length(GPS_latitude) && length(GPS_longitude) == length(GPS_TS),
                msg = "GPS_longitude, GPS_latitude, and GPS_TS vectors must have the same length")
    
    # Validate the format of GPS_TS
    assert_that(all(format(GPS_TS, "%Y-%m-%d %H:%M:%S") == GPS_TS), 
                msg = "GPS_TS must have a consistent timestamp format (YYYY-MM-DD HH:MM:SS)")
    
    # If ACC_TS is supplied, ensure it is of type POSIXct
    if (!is.null(ACC_TS)) {
      assert_that(is.POSIXct(ACC_TS), msg = "ACC_TS must be of type POSIXct")
    }
    assert_that(all(!is.na(ACC_TS)), msg = "ACC_TS cannot contain missing data")
    
    # If VeDBA is supplied, ensure it is numeric and matches the length of ACC_TS
    if (!is.null(VeDBA)) {
      assert_that(is.numeric(VeDBA), msg = "VeDBA must be a numeric vector")
      assert_that(length(VeDBA) == length(ACC_TS), msg = "VeDBA must have the same length as ACC_TS if both are supplied")
    }
    
    # Ensure drop_out is a single numeric whole number
    assert_that(is.numeric(drop_out), length(drop_out) == 1, drop_out == round(drop_out), 
                msg = "drop_out threshold must be a single whole numeric value")
    
    # Ensure a valid option is selected for burst_method
    valid_methods <- c("median", "mean", "last", "none")
    assert_that(burst_method %in% valid_methods, msg = "burst_method must be one of 'median', 'mean', 'last', or 'none'")
    
    # Ensure burst_len is a single numeric whole number
    assert_that(is.numeric(burst_len), length(burst_len) == 1, burst_len == round(burst_len), 
                msg = "burst_len threshold must be a single whole numeric value")
    
    # If standardise_time_interval is supplied, ensure it is a single numeric value > 1
    if (!is.null(standardise_time_interval)) {
      assert_that(is.numeric(standardise_time_interval), length(standardise_time_interval) == 1, 
                  standardise_time_interval > 1, msg = "standardise_time_interval must be a single numeric value > 1")
    }
    
    # Ensure standardise_universal is a logical value
    assert_that(is.logical(standardise_universal), length(standardise_universal) == 1, msg = "standardise_universal must be a single logical value (TRUE or FALSE)")
    
    # Ensure IF_conf is between 0 and 1 (a confidence level)
    assert_that(is.numeric(IF_conf), length(IF_conf) == 1, IF_conf > 0, IF_conf < 1, 
                msg = "IF_conf must be a numeric value between 0 and 1")
    
    # Ensure iso_sample_size is a single numeric whole number > 1
    assert_that(is.numeric(iso_sample_size), length(iso_sample_size) == 1, iso_sample_size == round(iso_sample_size), 
                iso_sample_size > 1, msg = "iso_sample_size must be a single whole numeric value greater than 1")
    
    # If GPS_height is supplied, ensure it matches the length of GPS data
    if (!is.null(GPS_height)) {
      assert_that(is.numeric(GPS_height), length(GPS_height) == length(GPS_longitude),
                  msg = "GPS_height must be a numeric vector of the same length as GPS_longitude, GPS_latitude, and GPS_TS")
    }
    
    # If start_timestamp is supplied, ensure it is of type POSIXct
    if (!is.null(start_timestamp)) {
      start_timestamp <- as.POSIXct(start_timestamp, format = "%Y-%m-%d %H:%M:%OS")
      assert_that(is.POSIXct(start_timestamp), msg = "start_timestamp must be of type POSIXct")
    }

    # If end_timestamp is supplied, ensure it is of type POSIXct
    if (!is.null(end_timestamp)) {
      end_timestamp <- as.POSIXct(end_timestamp, format = "%Y-%m-%d %H:%M:%OS")
      assert_that(is.POSIXct(end_timestamp), msg = "end_timestamp must be of type POSIXct")
    }
    
    # Ensure max_speed is a single numeric value
    assert_that(is.numeric(max_speed), length(max_speed) == 1, msg = "max_speed must be a single numeric value")
    
    # Ensure GPS_accuracy is a single numeric value > 1
    assert_that(is.numeric(GPS_accuracy), length(GPS_accuracy) == 1, GPS_accuracy > 1, 
                msg = "GPS_accuracy must be a single numeric value > 1")
    
    # Ensure plot is a logical value
    assert_that(is.logical(plot), length(plot) == 1, msg = "plot must be a single logical value (TRUE or FALSE)")
  }
  
  # Input validation
   validate_inputs(GPS_TS, GPS_longitude, GPS_latitude, ACC_TS, VeDBA, 
                   drop_out, burst_method, burst_len, standardise_time_interval,
                   standardise_universal, IF_conf, iso_sample_size, GPS_height, 
                   start_timestamp, end_timestamp, max_speed, GPS_accuracy, plot)
   
  ################################################################################################################################################
  
  ####Required functions###
  
  #Haversine distance formula 
  disty = function(long1, lat1, long2, lat2) { #GPS_longitude and GPS_latitude supplied in degrees
    long1 = long1 * pi/180 ; long2 = long2 * pi/180 ; lat1 = lat1 * pi/180 ; lat2 = lat2 * pi/180 #Function converts to radians
    a = sin((lat2 - lat1) / 2) * sin((lat2 - lat1) / 2) + cos(lat1) * cos(lat2) * sin((long2 - long1) / 2) * sin((long2 - long1) / 2)
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    d1 = 6378137 * c
    return(d1)
  }
  
  #Bearing function --> returns degrees - Great circular bearing between 2D positions
  beary = function(long1, lat1, long2, lat2) { #GPS_longitude and GPS_latitude supplied in degrees
    long1 = long1 * pi/180 ; long2 = long2 * pi/180 ; lat1 = lat1 * pi/180 ; lat2 = lat2 * pi/180 #Function converts to radians
    a = sin(long2 - long1) * cos(lat2)
    b = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(long2 - long1)
    c = ((atan2(a, b) / pi) * 180)  #Units returned in degrees (-180 to +180 degree scale)
    return(c)
  }
  
  ################################################################################################################################################
  message("Processing GPS data...")
  ######################################### (1) Create the initial data frame, and process data  ########################################
  
  Observation <- seq_along(GPS_TS) # Row number used for indexing and merging data frames
  
  df = data.frame(Observation, GPS_TS, GPS_longitude, GPS_latitude) ; colnames(df) = c("Observation", "Timestamp", "GPS_longitude", "GPS_latitude")
  
  # Check if GPS_height is provided and add it to the data frame if it's not NULL
  if (!is.null(GPS_height)) {
    df$GPS_height <- GPS_height
  }
  # Check if VeDBA is present and add as pseudo column
  if (!is.null(VeDBA)) {
    df$VeDBA <- NA
  }
  
  #Remove NA fixes
  df$GPS_longitude = ifelse(df$GPS_longitude == 0, NA, df$GPS_longitude) # In case missing coordinates are filled as zeros, replace with NA's
  df$GPS_latitude = ifelse(df$GPS_latitude == 0, NA, df$GPS_latitude) # In case missing coordinates are filled as zeros, replace with NA's
  df = df[!with(df, is.na(GPS_longitude) | is.na(GPS_latitude)) ,] # Remove NA fixes
  
  # Initialize processed GPS columns with original values
  df$GPS_longitude_filtered <- df$GPS_longitude
  df$GPS_latitude_filtered <- df$GPS_latitude
  
  df <- df %>%
    arrange(Timestamp) %>% # Arrange timestamp 
    mutate(
      Time_diff = c(0, as.numeric(diff(Timestamp))), # Create a time difference (s) between values
      Fix_number = cumsum(c(TRUE, Time_diff[-1] > burst_len)), # Make a group column which increments by one each GPS burst.
      Window_group = cumsum(c(TRUE, Time_diff[-1] >= drop_out)) # Make a group column ('Window_group') which increments by one each time the GPS drops-out (missing fixes) >= the user-defined 'drop_out' threshold (s)
    ) %>%
    group_by(Fix_number) %>%
    mutate(
      orig_burst_length = n(),
      process_burst = orig_burst_length <= burst_len
    ) %>%
    ungroup()
  
  # Handle burst processing as per selected method
  if (burst_method == "last") {
    df <- df %>%
      group_by(Fix_number) %>%
      mutate(
        GPS_latitude_filtered = if_else(
          process_burst & row_number() == n(),
          last(GPS_latitude),
          if_else(process_burst, NA_real_, GPS_latitude_filtered)
        ),
        GPS_longitude_filtered = if_else(
          process_burst & row_number() == n(),
          last(GPS_longitude),
          if_else(process_burst, NA_real_, GPS_longitude_filtered)
        )
      ) %>%
      ungroup()
  } else if (burst_method %in% c("mean", "median")) {
    df <- df %>%
      group_by(Fix_number) %>%
      mutate(
        middle_index = ceiling(n() / 2),
        GPS_latitude_filtered = if_else(
          process_burst & row_number() == middle_index,
          if (burst_method == "mean") mean(GPS_latitude) else median(GPS_latitude),
          if_else(process_burst, NA_real_, GPS_latitude_filtered)
        ),
        GPS_longitude_filtered = if_else(
          process_burst & row_number() == middle_index,
          if (burst_method == "mean") mean(GPS_longitude) else median(GPS_longitude),
          if_else(process_burst, NA_real_, GPS_longitude_filtered)
        )
      ) %>%
      ungroup() %>%
      select(-middle_index)
  }
  
  # Standardise GPS intervals?
  if (!is.null(standardise_time_interval)) {
    # Standardise just 1 Hz periods 
    if(standardise_universal == FALSE){
      # Subset to rows with non-NA filtered GPS data
      df_filtered <- df %>%
        filter(!is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered)) %>%
        group_by(Fix_number) %>%
        mutate(
          # Calculate cumulative time from the first fix within each Fix_number group
          cumulative_time = cumsum(c(0, as.numeric(diff(Timestamp)))), # Time from start of burst
          # Assign each fix to the nearest time interval
          interval_marker = floor(cumulative_time / standardise_time_interval),  # Nearest interval
          time_diff_to_interval = abs(cumulative_time - interval_marker * standardise_time_interval)  # Difference
        ) %>%
        arrange(interval_marker, time_diff_to_interval) %>%  # Order within intervals by proximity
        group_by(Fix_number, interval_marker) %>%
        slice_head(n = 1) %>%  # Keep only the closest fix to each interval
        ungroup()
    } else {
      # Standardize whole GPS dataset to a lower resolution
      df_filtered <- df %>%
        filter(!is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered)) %>%
        mutate(
          cumulative_time = cumsum(c(0, as.numeric(diff(Timestamp)))), # Time from start of dataset
          interval_marker = floor(as.numeric(difftime(Timestamp, min(Timestamp), units = "secs")) / standardise_time_interval),
          time_diff_to_interval = abs(as.numeric(difftime(Timestamp, min(Timestamp), units = "secs")) - 
                                        interval_marker * standardise_time_interval)
        ) %>%
        arrange(interval_marker, time_diff_to_interval) %>% # Order by proximity
        group_by(interval_marker) %>%
        slice_head(n = 1) %>%  # Keep only the closest fix
        ungroup()
    }
      # Update original data with filtered results
      df <- df %>%
        left_join(
          df_filtered %>% select(Timestamp, GPS_longitude_filtered, GPS_latitude_filtered),
          by = "Timestamp",
          suffix = c("", "_new")
        ) %>%
        mutate(
          GPS_longitude_filtered = ifelse(!is.na(GPS_longitude_filtered_new), GPS_longitude_filtered_new, NA),
          GPS_latitude_filtered = ifelse(!is.na(GPS_latitude_filtered_new), GPS_latitude_filtered_new, NA)
        ) %>%
        dplyr::select(-GPS_longitude_filtered_new, -GPS_latitude_filtered_new)
      
      # Recompute Fix_number for non-NA filtered fixes
      df <- df %>%
        mutate(Fix_number = cumsum(!is.na(GPS_longitude_filtered)))
    }
  
  ############### Activity processing --> #############  VeDBA supplied? #####################
  
  # Create an ACC dataframe with the same column structure as df
  if (!is.null(VeDBA) & !is.null(ACC_TS) & length(VeDBA) == length(ACC_TS)) {
    message("Processing Acceleration data...")
    acc_df <- data.frame(matrix(NA, nrow = length(ACC_TS), ncol = ncol(df)))
    colnames(acc_df) <- colnames(df)  # Set column names to match GPS dataframe
    # Fill specific columns with ACC data
    acc_df$Timestamp <- ACC_TS
    acc_df$VeDBA <- VeDBA
    
    # Combine the two dataframes
    df <- bind_rows(df, acc_df) %>%
      arrange(Timestamp)  # Sort by time
    
    # Create a column for fix groups
    df <- df %>%
      mutate(
        fix_group = cumsum(!is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered))
      )
    
    # Step 1: Reverse the data
    df_reversed <- df %>%
      arrange(desc(Timestamp)) %>% # Reverse the order
      mutate(
        fix_group = cumsum(!is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered)) # Group by GPS fixes
      )
    
    # Step 2: Compute the mean VeDBA for each group
    df_reversed <- df_reversed %>%
      group_by(fix_group) %>%
      mutate(
        mean_VeDBA_interval = ifelse(
          !is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered),
          mean(VeDBA, na.rm = TRUE), # Calculate mean VeDBA for the backward interval
          NA
        )
      ) %>%
      ungroup()
    
    # Step 3: Restore the original order
    df <- df_reversed %>%
      arrange(Timestamp) %>% # Restore original order
      select(-fix_group) # Drop the temporary grouping column
    
    df = subset(df, !is.na(df$Observation)) # Revert df back to normal after ACC processing
    df$VeDBA <- NULL
    
  }
  
  # Subset data between start and end timestamps, if provided
  message("Subsetting data based on provided start and/or end timestamps...")
  if (!is.null(start_timestamp) || !is.null(end_timestamp)) {
    df <- df %>%
      filter(
        if (!is.null(start_timestamp)) Timestamp >= start_timestamp else TRUE,
        if (!is.null(end_timestamp)) Timestamp <= end_timestamp else TRUE
      )
    
    if (nrow(df) == 0) {
      stop("No data remains after subsetting with the given timestamps.")
    }
  }
  
  ################################################################################################################################################
  message("Computing movement metrics to use in anomaly detection...")
  ######################################### (2) Speed and distance parameters  ########################################
  
  # Create a subset with retained filtered GPS values
  df_subset <- df %>%
    filter(!is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered)) %>%
    mutate(Time_diff = c(0, as.numeric(diff(Timestamp)))) %>%
    dplyr::rename(
      Longitude = GPS_longitude_filtered,
      Latitude = GPS_latitude_filtered
    )
  
  ### Identify erroneous spikes in your GPS data based on speed thresholds (incoming and outgoing (m/s)) and the turning angle (0-180°) between three consecutive points.
  
  df$Time_diff_filtered <- NA
  df$Ang_vertex <- NA
  df$Outgoing_speed <- NA
  df$Incoming_speed <- NA
  df$Dist_circular <- NA
  
  df_subset <- df_subset %>%
    mutate(
      # Shift longitude and latitude rows forwards (lead) and backwards (lag) by one row per group
      Longitude.lag = dplyr::lag(Longitude, n = 1, default = NA),
      Latitude.lag = dplyr::lag(Latitude, n = 1, default = NA),
      Longitude.lead = dplyr::lead(Longitude, n = 1, default = NA),
      Latitude.lead = dplyr::lead(Latitude, n = 1, default = NA),
      # Calculate bearing (angle) between 'Longitude'/'Latitude' and lag/lead coordinates using the beary function
      Ang.lag = beary(Longitude, Latitude, Longitude.lag, Latitude.lag),
      Ang.lead = beary(Longitude, Latitude, Longitude.lead, Latitude.lead),
      # Ensure angles are between 0 and 360 degrees
      Ang.lag = ifelse(Ang.lag < 0, Ang.lag + 360, Ang.lag),
      Ang.lead = ifelse(Ang.lead < 0, Ang.lead + 360, Ang.lead),
      # Calculate turning angle (0-180°) between every combination of three consecutive GPS fixes
      Ang.vertex = abs((Ang.lead - Ang.lag + 360) %% 360 - 180),
      # Calculate the incoming and outgoing speeds either side of angle vertex
      Dist.lag = disty(Longitude, Latitude, Longitude.lag, Latitude.lag),
      Dist.lead = disty(Longitude, Latitude, Longitude.lead, Latitude.lead),
      # Lead of time difference
      Time.diff.lag = as.numeric(difftime(Timestamp, dplyr::lag(Timestamp), units = "secs")),
      Time.diff.lead = dplyr::lead(Time.diff.lag, n = 1, default = NA),
      # Outgoing and incoming speed (m/s)
      Outgoing.speed = Dist.lag / Time.diff.lag,
      Incoming.speed = Dist.lead / Time.diff.lead
    ) %>%
    ungroup()
  
  # Experimental metrics
  df_subset <- df_subset %>%
    mutate(
      # Distance between inbound and outbound locations
      Dist_circular = disty(Longitude.lag, Latitude.lag, Longitude.lead, Latitude.lead),
    )
  
  # Merge computed values from df_subset to df where GPS coordinates are not NA
  df$Time_diff_filtered[!is.na(df$GPS_latitude_filtered)] <- df_subset$Time_diff
  # Forward-fill time differences within groups defined by non-NA filtered fixe
  df$Time_diff_filtered = zoo::na.locf(df$Time_diff_filtered, na.rm = FALSE) # Forward-fill
  df$Time_diff_filtered = zoo::na.locf(df$Time_diff_filtered, fromLast = TRUE) # Backward-fill for unfiltered rows
  df$Ang_vertex[!is.na(df$GPS_latitude_filtered)] <- df_subset$Ang.vertex
  df$Outgoing_speed[!is.na(df$GPS_latitude_filtered)] <- df_subset$Outgoing.speed
  df$Incoming_speed[!is.na(df$GPS_latitude_filtered)] <- df_subset$Incoming.speed
  df$Dist_circular[!is.na(df$GPS_latitude_filtered)] <- df_subset$Dist_circular
  
  
  # Keep only the desired columns 
  # Define the columns that will always be included
  base_columns <- c("Observation", "Timestamp", "Time_diff", "GPS_longitude", "GPS_latitude", 
                    "Fix_number", "Window_group", "orig_burst_length", 
                    "GPS_longitude_filtered", "GPS_latitude_filtered",
                    "Time_diff_filtered", "Ang_vertex", "Outgoing_speed", 
                    "Incoming_speed", "Dist_circular")
  
  # Check if GPS_height and VeDBA are present and include it if available
  if (!is.null(GPS_height)) {
    base_columns <- c(base_columns, "GPS_height")
  }
  if (!is.null(VeDBA)) {
    base_columns <- c(base_columns, "mean_VeDBA_interval")
  }
  
  # Replace NaN with NA in the entire data frame
  df <- df %>%
    mutate(across(where(is.numeric), ~ ifelse(is.nan(.), NA, .)))
  
  # Keep only the desired columns
  df <- df %>%
    select(all_of(base_columns))
  
  ################################################################################################################################################
  message("Beginning Isolation forest anomaly detection...")
  ######################################### (3) Isolation forest component of function  ########################################
  
  # Subset data to include just the fixes
  df_subset <- df %>%
    filter(!is.na(GPS_longitude_filtered) & !is.na(GPS_latitude_filtered))
  
  # Check if default isolation forest sample size settings need to be adjusted
  if (nrow(df_subset) < iso_sample_size) {
    sample_size <- nrow(df_subset)
  } else {
    sample_size <- iso_sample_size
  }
  
  # Begin anomaly detection
  message(paste("Using sample size:", sample_size, "for Isolation Forest anomaly detection"))
  if (!all(diff(df_subset$Timestamp) > 0)) {
    stop("Timestamps may be out of order - Ensure that timestamps are in ascending order")
  }
  
  # List of potential variables
  potential_vars <- c("Ang_vertex", "Outgoing_speed", "Incoming_speed", 
                      "GPS_height", "mean_VeDBA_interval")
  
  # Identify available variables in df
  available_vars <- intersect(potential_vars, names(df_subset))
  
  df_sub_data <- df_subset[, available_vars, drop = FALSE]
  
  # Remove columns with all NA values
  df_sub_data <- df_sub_data[, colSums(!is.na(df_sub_data)) > 0, drop = FALSE]
  
  # Check the number of available variables
  available_vars_count <- ncol(df_sub_data)
  
  # Determine ndim value
  ndim_value <- max(1, round(sqrt(available_vars_count)+1)) 
  
  # Build the Isolation Forest model
  iforest_model <- isolation.forest(
    df_sub_data,
    ntry=20,
    ndim = ndim_value,        #Multidimensioanl splines
    ntrees = 1000,            # Number of trees
    sample_size = iso_sample_size,       # Sub-sampling
    nthreads = 0,     
    standardize_data = TRUE, # Standardize data (default)
    missing_action = "impute", # Handle missing values
    min_gain = .25, # Minimum gain threshold
    prob_pick_pooled_gain = 0.75, # Quite aggressive
    penalize_range = TRUE
  )
  
  # Predict anomaly scores
  anomaly_scores <- predict(iforest_model, df_sub_data, type = "score")
  
  # Add anomaly scores to your original dataframe
  df_subset$IF_anomaly_score <- anomaly_scores
  
  # Label results
  df_subset$Verdict_IF <- ifelse(df_subset$IF_anomaly_score >= 
                                   as.numeric(quantile(df_subset$IF_anomaly_score, IF_conf, na.rm = TRUE)), 
                                 "Anomalous", "Not Anomalous")
  
  # Merge results back to main df
  df$Verdict_IF <- "Not Anomalous" 
  df$Verdict_IF <- ifelse(
    !is.na(df$GPS_latitude_filtered), 
    df_subset$Verdict_IF[match(df$Observation, df_subset$Observation)], 
    "Not Anomalous"
  )
  
  message("Isolation forest anomaly detection finished...")
  ############################
  # Biological relevant corrections
  
  if (!is.null(VeDBA)) {
    # Correct for genuine high, directed travelling movement 
    df$Verdict_IF <- ifelse(
      df$Verdict_IF == "Anomalous" & 
        !is.na(df$mean_VeDBA_interval) & 
        df$mean_VeDBA_interval >= quantile(df$mean_VeDBA_interval, 0.95, na.rm = TRUE) &
        (!is.na(df$Outgoing_speed) & df$Outgoing_speed >= quantile(df$Outgoing_speed, 0.95, na.rm = TRUE) | 
           !is.na(df$Incoming_speed) & df$Incoming_speed >= quantile(df$Incoming_speed, 0.95, na.rm = TRUE)) &
        !is.na(df$Ang_vertex) & df$Ang_vertex < 150,  # Additional condition for Ang_vertex
      "Not Anomalous",
      df$Verdict_IF
    )
    # Correct for genuine spikes in GPS during rest
    df$Verdict_IF <- ifelse(
      df$Verdict_IF == "Not Anomalous" & 
        !is.na(df$Ang_vertex) & df$Ang_vertex >= 170 & 
        !is.na(df$Dist_circular) & df$Dist_circular <= GPS_accuracy & 
        !is.na(df$mean_VeDBA_interval) & df$mean_VeDBA_interval < quantile(df$mean_VeDBA_interval, 0.25, na.rm = TRUE),
      "Anomalous",
      df$Verdict_IF
    )
  } else {
    # Correct for genuine high, directed travelling movement (without VeDBA)
    df$Verdict_IF <- ifelse(
      df$Verdict_IF == "Anomalous" & 
        (!is.na(df$Outgoing_speed) & df$Outgoing_speed >= quantile(df$Outgoing_speed, 0.95, na.rm = TRUE) | 
           !is.na(df$Incoming_speed) & df$Incoming_speed >= quantile(df$Incoming_speed, 0.95, na.rm = TRUE)) &
        !is.na(df$Ang_vertex) & df$Ang_vertex < 150,  # Additional condition for Ang_vertex
      "Not Anomalous",
      df$Verdict_IF
    )
    # Correct for genuine spikes during rest (without VeDBA)
    df$Verdict_IF <- ifelse(
      df$Verdict_IF == "Not Anomalous" & 
        !is.na(df$Ang_vertex) & df$Ang_vertex >= 170 & 
        !is.na(df$Dist_circular) & df$Dist_circular <= GPS_accuracy & 
        (!is.na(df$Outgoing_speed) & df$Outgoing_speed <= quantile(df$Outgoing_speed, 0.05, na.rm = TRUE) |
           !is.na(df$Incoming_speed) & df$Incoming_speed <= quantile(df$Incoming_speed, 0.05, na.rm = TRUE)),
      "Anomalous",
      df$Verdict_IF
    )
  }
  ###########################################################################
  
  # Ensure extreme spikes are anomalous
  df$Verdict_IF <- ifelse(
    df$Verdict_IF == "Not Anomalous" & 
      !is.na(df$Ang_vertex) & df$Ang_vertex >= 175 & 
      !is.na(df$Dist_circular) & df$Dist_circular <= quantile(df$Dist_circular, 0.05, na.rm = TRUE),
    "Anomalous",
    df$Verdict_IF
  )
  # Lastly, ensure no completely unrealistic fix speeds remained after above alterations
  df$Verdict_IF <- ifelse(
    df$Verdict_IF == "Not Anomalous" & 
      !is.na(df$Outgoing_speed) & df$Outgoing_speed >= max_speed,
    "Anomalous",
    df$Verdict_IF
  )
  
  ################################################################################################################################################
  message("Preparing summary plots...")
  ######################################### (4) Plot results and return data  ########################################
  
  old.par <- par(no.readonly = TRUE) # Save old graphical parameters
  plot_results <- function(df_subset, plot = TRUE, VeDBA = NULL, GPS_height = NULL) {
    
    if (plot == TRUE) {
      
      # Determine which variables are to be plotted
      plots_to_include <- list(
        Ang_vertex = list(data = df_subset$Ang_vertex, title = "Angle Between 3 Fixes (°)", xlab = "Angle Between 3 Fixes (°)"),
        Speeds = list(data = c(df_subset$Outgoing_speed, df_subset$Incoming_speed), title = "Outgoing/Incoming Speed (m/s)", xlab = "Speed (m/s)"),
        VeDBA = if (!is.null(VeDBA)) list(data = df_subset$mean_VeDBA_interval, title = "VeDBA", xlab = "VeDBA (g)") else NULL,
        GPS_height = if (!is.null(GPS_height)) list(data = df_subset$GPS_height, title = "GPS Height (m)", xlab = "Height (m)") else NULL
      )
      
      # Filter non-NULL plots
      plots_to_include <- Filter(Negate(is.null), plots_to_include)
      
      # Determine layout dynamically
      num_plots <- length(plots_to_include)
      num_rows <- ceiling(sqrt(num_plots))
      num_cols <- ceiling(num_plots / num_rows)
      
      # Set up plot layout
      par(mfrow = c(num_rows, num_cols), mar = c(4, 4, 2, 1)) # Adjust margins for better readability
      
      # Helper function for plotting
      plot_histogram <- function(data, title, xlab) {
        hist(data, breaks = "Scott", main = title, xlab = xlab, ylab = "Density",
             cex.lab = 1.2, freq = FALSE, col = "lightblue")
        abline(v = quantile(data, c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE), 
               col = c("red", "green", "green", "red"), lwd = 2, lty = 2)
        legend("topright", legend = c("1%", "5%", "95%", "99%"),
               col = c("red", "green", "green", "red"), lty = 2, lwd = 2, bty = "n", cex = 0.8)
      }
      
      # Loop through plots and generate
      for (plot_info in plots_to_include) {
        plot_histogram(plot_info$data, plot_info$title, plot_info$xlab)
      }
      
      # Restore old graphical parameters
      par(old.par)
      
      # Save the summary plot
      Thresholds <- recordPlot()
      return(Thresholds)
    }
  }
  
  Thresholds <- plot_results(df_subset, plot = TRUE, VeDBA = df_subset$mean_VeDBA_interval, GPS_height = df_subset$GPS_height)
  
  #####################################################################################################################
  # Dynamic leaflet map
  
  if (plot == TRUE) {
    ###########
    # Create a color palette for the Verdict_IF column
    verdict_palette <- colorFactor(
      palette = c("red", "green"),
      domain = c("Not Anomalous", "Anomalous")
    )
    
    df.sub.proc <- subset(df, !is.na(GPS_latitude_filtered) & !is.na(GPS_longitude_filtered))
    
    # Base map
    leaflet_map <- leaflet(data = df, options = leafletOptions(maxZoom = 22)) %>%
      addProviderTiles(providers$OpenStreetMap) %>%
      # Add a scale bar
      addScaleBar(
        position = "bottomleft",
        options = scaleBarOptions(
          maxWidth = 100,   # Maximum width of the scale bar
          metric = TRUE,    # Display in metric units
          imperial = FALSE  # Do not display in imperial units
        )
      ) %>%
      # Add the original, unfiltered GPS track (yellow points and lines)
      addPolylines(
        lng = ~GPS_longitude,
        lat = ~GPS_latitude,
        color = "black",
        weight = 0.7,
        opacity = 0.7,
        group = "Unfiltered GPS Track"
      ) %>%
      addCircleMarkers(
        lng = ~GPS_longitude,
        lat = ~GPS_latitude,
        color = "black",
        stroke = FALSE,
        fillOpacity = 0.7,
        radius = 0.7, # Smaller circles
        popup = ~paste(
          "Timestamp: ", as.POSIXct(Timestamp, origin = "1970-01-01"), "<br>",
          "Observation: ", Observation, "<br>",
          "Longitude: ", GPS_longitude, "<br>",
          "Latitude: ", GPS_latitude, "<br>"
        ),
        group = "Unfiltered GPS Track"
      ) %>%
      # Add the filtered track for "Not Anomalous" (purple points and lines)
      addPolylines(
        data = df.sub.proc %>% filter(Verdict_IF == "Not Anomalous"),
        lng = ~GPS_longitude_filtered,
        lat = ~GPS_latitude_filtered,
        color = "orange",
        weight = 1,
        opacity = 0.7,
        group = "Filtered Track (Not Anomalous)"
      ) %>%
      addCircleMarkers(
        data = df.sub.proc %>% filter(Verdict_IF == "Not Anomalous"),
        lng = ~GPS_longitude_filtered,
        lat = ~GPS_latitude_filtered,
        color = "orange",
        stroke = FALSE,
        fillOpacity = 0.7,
        radius = 1.7, # Smaller circles
        popup = ~paste(
          "Timestamp: ", as.POSIXct(Timestamp, origin = "1970-01-01"), "<br>",
          "Observation: ", Observation, "<br>",
          "Longitude: ", round(GPS_longitude_filtered, 5), "<br>",
          "Latitude: ", round(GPS_latitude_filtered, 5), "<br>",
          "Outgoing speed (m/s): ", round(Outgoing_speed, 3), "<br>",
          "Ang vertex (°): ", round(Ang_vertex, 1), "<br>",
          "Incoming speed (m/s): ", round(Incoming_speed, 3), "<br>",
          "Dist_circular (m): ", round(Dist_circular, 3), "<br>",
          if (!is.null(GPS_height)) {
            paste("GPS height (m): ", round(GPS_height, 1), "<br>")
          } else {
            ""
          },
          if (!is.null(VeDBA)) {
            paste("VeDBA (g): ", round(mean_VeDBA_interval, 3), "<br>")
          } else {
            ""
          },
          "Verdict: ", Verdict_IF, "<br>"
        ),
        group = "Filtered Track (Not Anomalous)"
      ) %>%
      
      # Add a legend for the Verdict_IF column
      addLegend(
        "bottomright",
        pal = verdict_palette,
        values = ~Verdict_IF,
        title = "Anomaly Verdict",
        opacity = 1
      ) %>%
      
      # Add the filtered track (cyan lines and points, colored by Verdict_IF)
      addPolylines(
        data = df.sub.proc,
        lng = ~GPS_longitude_filtered,
        lat = ~GPS_latitude_filtered,
        color = "blue",
        weight = 1,
        opacity = 0.7,
        group = "Filtered GPS Track"
      ) %>%
      addCircleMarkers(
        data = df.sub.proc,
        lng = ~GPS_longitude_filtered,
        lat = ~GPS_latitude_filtered,
        color = ~verdict_palette(Verdict_IF),
        stroke = FALSE,
        fillOpacity = 0.7,
        radius = 1.7, # Smaller circles
        popup = ~paste(
          "Timestamp: ", as.POSIXct(Timestamp, origin = "1970-01-01"), "<br>",
          "Observation: ", Observation, "<br>",
          "Longitude: ", round(GPS_longitude_filtered, 5), "<br>",
          "Latitude: ", round(GPS_latitude_filtered, 5), "<br>",
          "Outgoing speed (m/s): ", round(Outgoing_speed, 3), "<br>",
          "Ang vertex (°): ", round(Ang_vertex, 1), "<br>",
          "Incoming speed (m/s): ", round(Incoming_speed, 3), "<br>",
          "Dist_circular (m): ", round(Dist_circular, 3), "<br>",
          if (!is.null(GPS_height)) {
            paste("GPS height (m): ", round(GPS_height, 1), "<br>")
          } else {
            ""
          },
          if (!is.null(VeDBA)) {
            paste("VeDBA (g): ", round(mean_VeDBA_interval, 3), "<br>")
          } else {
            ""
          },
          "Verdict: ", Verdict_IF, "<br>"
        ),
        group = "Filtered Track (Not Anomalous)"
      ) %>%
      # Add layer controls for toggling
      addLayersControl(
        overlayGroups = c(
          "Unfiltered GPS Track",
          "Filtered GPS Track",
          "Filtered Track (Not Anomalous)"
        ),
        options = layersControlOptions(collapsed = FALSE)
      )
    
    # Print the map
    print(leaflet_map)
    
    # Return results
    
    par(old.par) # Return plotting parameters back to normal
  }
  
  df$Verdict_IF <- ifelse(
    is.na(df$GPS_latitude_filtered), 
    "No Filtered GPS data", # Make rows with no GPS data not reflecting 'Not Anomalous' fixes
    df$Verdict_IF
  )
  rm(df_subset, df_sub_data, df.sub.proc)
  
  message("Done!")
  return(df) # Return main data frame
  
  gc() #Clear memory

}

# End of function
###################################################################################################################################################################################################################################
###################################################################################################################################################################################################################################

#Example:
#
# Load GPS data
#df <- read.csv("xxx.csv")
#df = subset(df, df$sensor.type == "gps")
#df = subset(df, is.na(df$timestamp) != TRUE)
#df$timestamp = as.POSIXct(df$timestamp, format = "%Y-%m-%d %H:%M:%S")
#head(df$timestamp, 1) ; tail(df$timestamp, 1)
#Remove duplicated timestamps
#df<-df[!duplicated(df[c("timestamp")]),]
# Ensure GPS data are sorted by time
#df <- df %>% arrange(timestamp) 

# Load acc data -->  Use the Eobs_Data_Reader function to processs and extract acceleration data from Eobs' devices
#setwd("xxx")
#df.acc <- read_fst("xxx.fst") #library(fst)
#df.acc$interpolated_timestamp = as.POSIXct(df.acc$interpolated_timestamp, format = "%Y-%m-%d %H:%M:%OS")
#head(df.acc$interpolated_timestamp, 1) ; tail(df.acc$interpolated_timestamp, 1)
# Filter out duplicate times and short bursts
#df.acc = subset(df.acc, df.acc$duplicate_times == FALSE & df.acc$standardized_burst_duration >= 5)

# Group by standardized_burst_id and compute mean VeDBA and median interpolated_timestamp
#df.acc_summary <- df.acc %>%
#  group_by(standardized_burst_id) %>% # Mean VeDBA value per burst ID (here, about 10 s bursts every 1 min)
 # summarise(
 #   mean_VeDBA = mean(VeDBA, na.rm = TRUE), # Mean VeDBA
 #   timestamp = median(interpolated_timestamp, na.rm = TRUE) # Median interpolated timestamp of the burst
 # ) %>% arrange(timestamp) %>% # Ensure timestamps are in correct chronological order
 # ungroup() 

#rm(df.acc) # Remove main acc data



#test = GPS_Screener(GPS_TS = df$timestamp,
              #   GPS_longitude = df$location.long, 
              #   GPS_latitude = df$location.lat, 
              #   ACC_TS = df.acc_summary$timestamp,
              #   VeDBA = df.acc_summary$mean_VeDBA,
              #   drop_out = 1800, 
              #   burst_method = "last",
              #   standardise_time_interval = 240, #4 min intervals
              #   standardise_universal = FALSE,
              #   burst_len = 10, #bursts were 10 s
              #   IF_conf = 0.99, #99 percentile for anomaly detection
              #   iso_sample_size = 500,
              #   GPS_height = df$height.above.ellipsoid,
              #   start_timestamp = "2024-03-14 21:40:37",
              #   end_timestamp = "2024-06-06 23:27:11",
              #   max_speed = 5,
              #   GPS_accuracy = 25,
              #   plot = TRUE)

      