#Add some metadata to the netCDF files
#Two variables: instrument_name and instrument_mfr are requested by readMSData from xcms,
#but they are missing from the original netCDF files so we add them to avoid errors.
#We only need to run this script once. The modified netCDF files are in the folder
#C:\\Users\\andre\\Desktop\\SUMMER_COURSE\\data
#The unmodified ones are in C:\\Users\\andre\\Desktop\\SUMMER_COURSE

library(ncdf4)

# Open the existing netCDF file
data_files <- list.files(path = "C:\\Users\\andre\\Desktop\\SUMMER_COURSE\\data", pattern = ".cdf", all.files = TRUE, full.names = TRUE)

for (i in 1:length(data_files))
{

nc_file <- nc_open(data_files[i], write = TRUE)

# Define a new variable for instrument_name
dim1 <- ncdim_def(name = "char_dim", units = "", vals = 1:10, create_dimvar = FALSE)
var1 <- ncvar_def(name = "instrument_name", units = "", dim = dim1, prec = "char")

dim2 <- ncdim_def(name = "char_dim", units = "", vals = 1:10, create_dimvar = FALSE)
var2 <- ncvar_def(name = "instrument_mfr", units = "", dim = dim2, prec = "char")

# Add the new variables to the file
nc_file <- ncvar_add(nc_file, var1)
nc_file <- ncvar_add(nc_file, var2)

# Write place-holder data to the new variable
ncvar_put(nc_file, var1, "Unknown")
ncvar_put(nc_file, var2, "Unknown")

# Close the file
nc_close(nc_file)

# Now try reading the modified file with readMSData
ms_data <- readMSData(data_files[i], mode = "onDisk")
print(data_files[i])
print(ms_data)
}
