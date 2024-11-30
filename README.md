
# GPS_Screener

The `GPS_Screener` function is a fast and versatile R-based tool for automatically detecting anomalies in GPS data. using a combination of movement metrics. It combines movement metrics, optional acceleration data (e.g., VeDBA), and isolation forest models to handle a wide range of GPS data complexities, including burst modes and 1 Hz recordings. Additionally, it offers visualization capabilities to inspect anomaly detection results interactively.

![Dynamic Plot Example](GPSScreenerLeaflet.png)

**[This function is in active development. Feedback, testing, and suggestions for improvement are highly encouraged!]**

---

## Key Features

### Core Functionality
- **Data Pre-Processing**:
  - Processes GPS burst modes and standardizes GPS intervals to a consistent sampling frequency, ensuring robust anomaly detection for diverse datasets (e.g., burst modes or 1 Hz recordings).
- **Anomaly Detection**:
  - Derives GPS movement metrics by analyzing three consecutive fixes, including outbound and inbound speed, vertex angles, and distance changes between the first and third fixes.
  - Combines within-function defined thresholds with unsupervised multi-dimensional isolation forest models.

- **Acceleration Data Integration** (see *'Eobs_Data_Reader'* function to process acceleration data from Eobs devices):
  - Leverages VeDBA (or equivalent) as an optional additional metric to improve anomaly detection. 
  - Matches acceleration data to GPS fixes, even with irregular or non-aligned timestamps, and calculates mean VeDBA for preceding intervals. (See the *Eobs_Data_Reader* function for pre-processing acceleration data from e-obs devices).
  - Requires both GPS and acceleration timestamps to be in the same POSIXct format (Y-m-d H:M:S or Y-m-d H:M:OS).

- **Interactive Visualization**:
  - Provides dynamic Leaflet maps to explore anomaly detection results, including filtering layers for intuitive data inspection.

### Additional Features
- Adaptive handling of missing or irregular data.
- Post-modelling refinement using biologically meaningful thresholds to improve accuracy.
- Optional integration of GPS height as an additional parameter for anomaly detection and timestamp subsetting.
- Detailed outputs, including processed data and diagnostic visuals, for comprehensive analysis.

---

## Input Parameters

### Required Inputs
1. **`GPS_TS`**: POSIXct vector of GPS timestamps (must be in Y-m-d H:M:S format). Missing values (NA) are not allowed.
2. **`GPS_longitude`**: Numeric vector of GPS longitude values (decimal degrees). May include NAs or zeros.
3. **`GPS_latitude`**: Numeric vector of GPS latitude values (decimal degrees). May include NAs or zeros.

### Optional Inputs
1. **`ACC_TS`**: POSIXct vector of acceleration timestamps. If supplied, it must match the format of GPS_TS (Y-m-d H:M:S or Y-m-d H:M:OS) and cannot contain NAs (default = NULL).
2. **`VeDBA`**: Numeric vector of VeDBA (or equivalent activity metric) values matching the length of `ACC_TS` (default = NULL).
3. **`drop_out`**: Numeric value (in seconds) defining the maximum temporal gap for identifying large drop-outs in GPS data (default = 3600). Although its functionality is not actively used in the function, the output includes a corresponding column for reference.
4. **`burst_method`**: Character string specifying the method for processing GPS bursts. Options include "median", "mean", "last", or "none" (default = "last").
5. **`burst_len`**: Numeric value (seconds) for defining the typical maximum burst durations (default = 1).
6. **`standardise_time_interval`**: Numeric value (seconds) for standardizing GPS fix intervals to a consistent time step (default = NULL).
7. **`standardise_universal`**: Logical value (TRUE/FALSE) indicating whether **`'standardise_time_interval'`** is applied across all GPS data or restricted to continuous 1 Hz periods (default = FALSE).
7. **`IF_conf`**: Numeric confidence level for isolation forest anomaly detection (default = 0.99).
8. **`iso_sample_size`**: Numeric value for the sample size used in the isolation forest model (default = 256).
9. **`GPS_height`**: Numeric vector of GPS height values matching the length of GPS data (default = NULL). This is used as an additional parameter within the anomaly detection.
10. **`start_timestamp`**: POSIXct value for subsetting data, defining the start of the time range (optional; but must be used with `end_timestamp`).
11. **`end_timestamp`**: POSIXct value for subsetting data, defining the end of the time range (optional; but must be used with `start_timestamp`).
12. **`max_speed`**: Numeric value defining the maximum biologically plausible speed (in m/s) for the species of interest (default = 5).
13. **`GPS_accuracy`**: Numeric value estimating the GPS error radius (in meters) under stationary conditions (default = 25). A value of 25 m is recommended based on empirical data for forested habitats. Recommended to not make this value smaller than 10 m.
14. **`plot`**: Logical value (TRUE/FALSE) indicating whether to generate diagnostic plots (default = TRUE).

---

## Outputs

### Core Output
A **data frame** containing:
1. **`Observation`**: A sequential integer vector that uniquely identifies each row in the dataset. Used as a reference index (relative to initial supplied data) for processing and merging data.
2. **`Timestamp`**: Timestamp of the fix
3. **`Time_diff`**: Time difference (in seconds) between consecutive GPS fixes.
4. **`GPS_longitude`** & **`GPS_latitude`**: Raw GPS coordinates (in decimal degrees).
5. **`Fix_number`**:  Relative order of each kept processed GPS fix per burst or continuous 1 Hz data collection session (unless **`standardise_time_interval`** was used).
6. **`Window_group`**: Grouping variable which incremented each time the cumualtive time exceeded the **`drop_out`** threshold.
7. **`orig_burst_length`**: The original length (in seconds) of each burst or continuous GPS session prior to processing.
8. **`GPS_longitude_filtered`** & **`GPS_latitude_filtered`**: GPS coordinates (in decimal degrees), corresponding to the processed values.
9. **`Time_diff_filtered`**: Represents the time interval (s) between filtered GPS fixes but is replicated across all rows associated with each fix group.
10. **`Ang_vertex`**: The turning angle (in degrees) at each fix, calculated using three consecutive GPS points. Values range from 0° to 180°, where higher angles indicate sharper turns.
11. **`Outgoing_speed`**: Movement speeds (in meters per second) calculated between the current GPS fix and the following one.
12. **`Incoming_speed`**: Movement speeds (in meters per second) calculated between the current GPS fix and the preceding one.
13. **`Dist_circular`**: The straight-line distance (in meters) between the GPS fix preceding the current one and the fix following the current one. Useful for identifying circular or looping movements.
14. **`GPS_height (optional)`**: vertical altitude values (in meters) associated with each GPS fix. This field is included only if the user provides height data as an input.
15. **`mean_VeDBA_interval (optional)`**: The mean VeDBA (or equivalent supplied metric) over the time interval leading up to each GPS fix. This field is included only if acceleration data is provided.
16. **`Verdict_IF`**: A categorical variable indicating whether a GPS fix is classified as **`"Anomalous"`** or **`"Not Anomalous"`** based on isolation forest analysis and custom rules.

### Optional Outputs (if `plot = TRUE`)
- **Interactive Map**:
  - Dynamic leaflet map with color-coded anomaly results.
  - Toggleable layers for unfiltered, filtered, and "Not Anomalous" tracks.
  - Scale bar for distance estimation.
- **Summary Plots**:
  - Histograms showing distributions of key metrics, annotated with quantile thresholds.

---

## Example Workflow

### Step 1: Load Your Data
```r
library(dplyr)
# GPS data
df <- read.csv("C:/Users/richard/xxxxxx/Gandalf.csv")
# Covert to POSIXct
df$timestamp = as.POSIXct(df$timestamp, format = "%Y-%m-%d %H:%M:%S")
head(df$timestamp, 1) ; tail(df$timestamp, 1)
# Ensure no duplicated time stamps
df<-df[!duplicated(df[c("timestamp")]),]
# make sure no NAs in timestamp
df = subset(df, !is.na(df$timestamp))
# Ensure GPS data are sorted by time
df <- df %>% arrange(timestamp) 

# Load ACC VeDBA data -->  Use the Eobs_Data_Reader function to processs and extract acceleration data from Eobs' devices
setwd("C:/Users/richard/xxxxxx/Eobs data [All studies]/Invisible networks/Processed data")
library(fst) # Fast, and easy way to serialize data frames when writing and reading data.
df.acc <- read_fst("Gandalf.processed.acc.fst")
# Covert to POSIXct
df.acc$interpolated_timestamp = as.POSIXct(df.acc$interpolated_timestamp, format = "%Y-%m-%d %H:%M:%OS")
head(df.acc$interpolated_timestamp, 1) ; tail(df.acc$interpolated_timestamp, 1)
# Filter out potentially duplicated timestamps and short ACC bursts  
df.acc = subset(df.acc, df.acc$duplicate_times == FALSE & df.acc$standardized_burst_duration >= 5)
# make sure no NAs in timestamp
df.acc = subset(df.acc, !is.na(df.acc$interpolated_timestamp))

# Group by standardized_burst_id and compute mean VeDBA and median interpolated_timestamp to obtain a single mean value per burst
df.acc_summary <- df.acc %>%
  group_by(standardized_burst_id) %>% # Mean VeDBA value per burst ID (here, about 10 s bursts every 1 min)
  summarise(
    mean_VeDBA = mean(VeDBA, na.rm = TRUE), # Mean VeDBA
    timestamp = median(interpolated_timestamp, na.rm = TRUE) # Median interpolated timestamp of the burst
  ) %>% arrange(timestamp) %>% # Ensure timestamps are in correct chronological order
  ungroup() 
```

### Step 2: Run the Function
```r
results <- GPS_Screener(GPS_TS = df$timestamp,
                 GPS_longitude = df$location.long, 
                 GPS_latitude = df$location.lat, 
                 ACC_TS = df.acc_summary$timestamp,
                 VeDBA = df.acc_summary$mean_VeDBA,
                 drop_out = 1800, 
                 burst_method = "last",
                 standardise_time_interval = 240, #4 min intervals
                 standardise_universal = FALSE,
                 burst_len = 10, #bursts were 10 s
                 IF_conf = 0.99, #99 percentile for anomaly detection
                 iso_sample_size = 500,
                 GPS_height = df$height.above.ellipsoid,
                 start_timestamp = "2024-05-04 10:31:00",
                 end_timestamp = "2024-05-31 19:59:11",
                 max_speed = 5,
                 GPS_accuracy = 25,
                 plot = TRUE)
)
```

---

## Required R Packages

The function installs and uses the following R packages:
- **zoo**
- **dplyr**
- **tidyr**
- **data.table**
- **ggplot2**
- **leaflet**
- **assertthat**
- **isotree**

---

## Limitations
- Input timestamps must consistently use the POSIXct format.
- Results are heavily influenced by the quality and consistency of the input data.
- A previously C++ implemented "distance from median" rolling time function—designed to calculate median fixes within a dynamically adjusted temporal window—was removed. This metric proved ineffective for datasets with temporally sparse fix intervals (e.g., >1 minute), especially for highly dynamic animal movement patterns. For people with high-res data sets, this maybe a useful metric. Contact for more info.

---

## License

This project is licensed under the MIT License.

## Contact

For questions, bug reports, suggestions, or contributions, please contact:
- Richard Gunner
- Email: rgunner@ab.mpg.de
- GitHub: [Richard6195](https://github.com/Richard6195)
