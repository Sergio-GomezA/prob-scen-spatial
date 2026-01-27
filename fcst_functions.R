

readwindyr <- function(year,
                       path = "../data/"
                       # skip.lines = 23,
                       # column_names = c("time","forecast","actuals",
                       # "load","hydro","thermal","netinter")
){
  extensions <- c(".xls",".xlsx")
  if (year>=2022) fext <- extensions[2] else fext <- extensions[1]
  file = paste(path,"WindGenTotalLoadYTD_",year,fext,sep = "")
  
  if (year>=2022){
    excerpt = read_xlsx(file, n_max = 30, col_names = F, .name_repair = "unique_quiet")
  } else {
    excerpt = read_xls(file, n_max = 30, col_names = F, .name_repair = "unique_quiet")
  }
  
  skip.lines <- which(excerpt[,1]=="Date/Time")-1
  
  if (year>=2022){
    table <- read_xlsx(file, skip = skip.lines,
                       col_names = T)
  } else{
    table <- read_xls(file, skip = skip.lines, sheet = 1,
                      col_names = T) %>% # read first semester
      bind_rows(
        read_xls(file, skip = skip.lines, sheet = 2,
                 col_names = T) # read second semester
      )
  }
  fix.names <- data.frame(long.name = colnames(table)) %>% 
    left_join(catalog,
              by = c("long.name"="colname"))
  colnames(table) <- fix.names$tag
  
  table <- table %>% 
    mutate(time=as.POSIXct(time, tz = "PST", format = "%m/%d/%y %H:%M")) %>% # time format
    group_by(time) %>% 
    summarise_all(first)
  table
}

list.plants.to.date <- function(end.date){
  plant.list <- capacity.list %>% 
    filter(date<=end.date) %>% 
    # remove extra characters
    mutate(
      plant=sub("Phase.*|PH.*|\\(phase.*|\\(OUT.*|@.*","",plant) %>% 
        trimws()
    ) %>% 
    group_by(plant) %>% 
    summarise(nameplate.capacity=sum(nameplate.capacity)) %>% 
    filter(nameplate.capacity>0) %>% 
    arrange(desc(nameplate.capacity))
  plant.list
}


get.wind.closes.avg <- function(date){
  #
  file <- sprintf("MERRA2_400.tavg1_2d_slv_Nx.%d%02d%02d.nc4.nc4",year(date),month(date),day(date))
  # read file
  nc <- nc_open(paste(merrapath,file,sep = ""))
  
  # get wind speed
  east <- ncvar_get(nc, "U50M") # store the data in a 3-dimensional array
  north <- ncvar_get(nc, "V50M") # store the data in a 3-dimensional array
  w.speed <- sqrt(east^2+north^2)
  nc_close(nc)
  
  # extract all data (inefficient memory usage)
  #Create 2D matrix of long, lat and time
  lonlattime <- as.matrix(expand.grid(lon,lat,t)) # this might take several seconds
  #reshape whole lswt_array
  lswt_vec_long <- as.vector(w.speed)
  length(lswt_vec_long) # 
  #Create data.frame
  lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
  rm(lswt_vec_long,lonlattime)
  colnames(lswt_obs) <- c("lon","lat","t","w.speed")
  
  # get wind speed for closest points
  daily.wind <- expand_grid(merra.grid.close,t) %>% 
    left_join(lswt_obs, by=c("lon","lat","t")) %>% 
    mutate(time = date + minutes(t)) %>% 
    # average wind speed
    group_by(time) %>% 
    summarise(w.speed = mean(w.speed))
  
  daily.wind
}




get.wind.interp <- function(date){
  #
  file <- merra.names$fname[merra.names$date==date]
  # read file
  nc <- nc_open(paste(merrapath,file,sep = ""))
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat", verbose = F)
  t <- ncvar_get(nc, "time")
  # get wind speed
  east <- ncvar_get(nc, "U50M") # store the data in a 3-dimensional array
  north <- ncvar_get(nc, "V50M") # store the data in a 3-dimensional array
  w.speed <- sqrt(east^2+north^2)
  nc_close(nc)
  
  # # extract all data (inefficient memory usage)
  # #Create 2D matrix of long, lat and time
  # lonlattime <- as.matrix(expand.grid(lon,lat,t)) # this might take several seconds
  # #reshape whole lswt_array
  # lswt_vec_long <- as.vector(w.speed)
  # length(lswt_vec_long) # 
  # #Create data.frame
  # lswt_obs <- data.frame(cbind(lonlattime, lswt_vec_long))
  # rm(lswt_vec_long,lonlattime)
  # colnames(lswt_obs) <- c("lon","lat","t","w.speed")
  
  # bilinear interpolation from fields package
  # wind speed for all longitudes, latitudes at time 0:00 
  control_dat <- list(x=lon, y=lat, z=w.speed[,,1])
  
  # interpolates w.speed. Gets 21 farms for 1st hour in nc file
  interp_var <- interp.surface(
    control_dat, 
    cbind(plants.bpa.loc$xlong,plants.bpa.loc$ylat)
    # plant.loc[,c("xlong","ylat")]
  )
  
  # bilinear interpolation
  interp_series <- apply(w.speed, MARGIN = 3, 
                         FUN = \(wind.array) {
                           control_dat <- list(x=lon, y=lat, z=wind.array)
                           interp_var <- interp.surface(control_dat, 
                                                        cbind(plants.bpa.loc$xlong,plants.bpa.loc$ylat))
                         },
                         simplify = T)
  
  # simple average
  ws.locm <- c(apply(interp_series, MARGIN = 2,
                     mean))
  
  # weighted average. capacity is fixed
  ws.locw <- c(apply(interp_series, MARGIN = 2,
                     FUN = \(wloc) {
                       sum(wloc*plants.bpa.loc$nameplate.capacity)/sum(plants.bpa.loc$nameplate.capacity)
                     }))
  
  daily.wind <- data.frame(
    t = t,
    ws.m = ws.locm,
    ws.w = ws.locw) %>% 
    mutate(time = as.Date(date) + minutes(t)) %>% 
    select(time,ws.m,ws.w)
  
  daily.wind
}