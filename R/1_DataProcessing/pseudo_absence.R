#' Author: TJ Sipin
#' Date: March 8, 2024
#' Purpose:
  #' Get pseudo-absence points for each month

set.seed(030824)
setwd("/home/tjsipin/Dengue_CA")

# Packages ----------------------------------------------------------------

library(terra)
library(raster)
library(purrr)
library(dplyr)
library(flexsdm)
tidymodels::tidymodels_prefer()
library(foreach)
library(doParallel)
library(ggplot2)
library(tidyterra)
library(readr)
library(tidyr)
library(sf)

# Data --------------------------------------------------------------------

# read in presence data
aedes_p_spdf = read.csv("Data/1_DataProcessing/Aedes/aedes_presence_mixedfemales.csv") %>% 
  select(longitude, latitude, Month) %>% 
  # set ID and presence values
  mutate(
    id = row_number(),
    pr_ab = 1
  ) %>% 
  # turn to spatial points data frame
  SpatialPointsDataFrame(
    coords = data.frame(.$longitude, .$latitude),
    data = .,
    proj4string = CRS("+proj=longlat +datum=WGS84")
  ) %>% 
  spTransform(CRS("EPSG:3310"))

aedes_p = aedes_p_spdf %>% 
  as_tibble() %>% 
  # make species column for thinning function
  mutate(species = "aedes") %>% 
  select(longitude, latitude, Month, id, pr_ab, species)

# path to environmental rasters
env_raster_path = "Data/1_DataProcessing/rasters/env/"
env_rasters = map(list.files(env_raster_path, full.names = T), rast)
env_stack = rast(env_rasters)[[c(5, 4, 8, 1, 9, 7, 6, 2, 12, 11, 10, 3)]] #%>% 
  # terra::project(crs(aedes_p_spdf))
  # terra::project(crs("+proj=longlat +datum=WGS84"))

# given month, get spatially thinned presences and write pseudo-absence locations as RDS
psuedoabsenceFunc = function(mon){
  tryCatch(
    {
      t0 = Sys.time()
      this_dir = paste0("Data/1_DataProcessing/pseudoabs/geo_env/", mon, "/")
      dir.create(this_dir)
      # get environmental raster for month
      this_env_rast = env_stack[[mon]]
      
      # subset presence points to month
      this_aedes_p = aedes_p %>% 
        filter(Month == mon) %>% 
        tibble()
      
      this_aedes_pf = spThin::thin(
        loc.data = this_aedes_p,
        lat.col = "latitude",
        long.col = "longitude",
        spec.col = "species",
        thin.par = 0.27,
        reps = 10,
        write.files = F,
        locs.thinned.list.return = T,
        out.dir = paste0(this_dir, "thinned"),
        out.base = paste0(mon, "_pf")
      ) %>% 
        bind_rows() %>% 
        distinct() %>% 
        mutate(pr_ab = 1) %>% 
        SpatialPointsDataFrame(
          coords = data.frame(.$Longitude, .$Latitude),
          data = .,
          proj4string = CRS("+proj=longlat") ## try1
          # proj4string = CRS("EPSG:3310") ## try2
        ) %>% 
        spTransform(CRS("EPSG:3310")) %>% 
        as.data.frame() %>% 
        rename(x = coords.x1, y = coords.x2)
      
      # create partition blocks
      ### should these parameters be changed?
      print(paste0(mon, ": occurrence partition"))
      this_occ_part <- this_aedes_pf %>%
        part_sblock(
          data = .,
          env_layer = this_env_rast,
          pr_ab = "pr_ab",
          x = "x", ##try2
          y = "y", 
          prop = 0.5,
          n_part = 3,
          min_res_mult = 200,
          max_res_mult = 4000,
          num_grids = 100
        )
      
      saveRDS(this_occ_part, paste0(this_dir, mon, "_occ_part.rds"))
      
      this_part = this_occ_part$part %>% 
        mutate(species = "aedes")
      
      this_part_f = this_aedes_pf %>% 
        select(x, y) %>% 
        left_join(
          # this_part %>% select(Longitude=x, Latitude=y, .part, pr_ab) ##try1
          this_part %>% select(x, y, .part, pr_ab) ##try2
        ) %>% 
        distinct()
      
      print(paste0(mon, ": calibration"))
      
      cali2 = maps::map('state', plot=F, fill=T, regions="california") %>%
        st_as_sf() %>%
        # st_transform("+proj=longlat +datum=WGS84")
        st_transform(CRS("EPSG:3310"))
      
      cali2_vect = terra::vect(cali2)
      
      this_ca <-
        calib_area(
          data = this_part_f,
          # x = 'Longitude', ##try1
          # y = 'Latitude',
          x = "x", ##try2
          y = "y",
          # defined by buffered minimum convex polygon
          method =  c("mask", cali2_vect, "ID"),
          # method = "mcp",
          crs = crs(this_env_rast)
        ) # create a calibration area with 15 km buffer around occurrence points
      
      
      print(paste0(mon, ": block_layer"))
      this_block_layer = get_block(
        env_layer = this_env_rast, best_grid = this_occ_part$grid
      )
      
      saveRDS(this_block_layer, paste0(this_dir, mon, "_block_layer.rds"))
      
      this_aedes_psa = lapply(sort(unique(this_part_f$.part)), function(x){
        sample_pseudoabs(
          # data = this_part_f %>% rename(x=Longitude, y=Latitude), ##try1
          data = this_part_f, ##try2
          x = "x",
          y = "y",
          n = sum(this_part_f$.part == x),
          method = c("geo_env_const", width = "10000", env = this_env_rast), #10km buffer, is this ok?
          rlayer = this_block_layer,
          # maskval = x,
          calibarea = this_ca
        )
      }) %>% 
        bind_rows() %>% 
        mutate(Month = unique(this_aedes_p$Month))
      
      saveRDS(this_part_f, paste0(this_dir, mon, "_partition_filtered.rds"))
      saveRDS(this_aedes_psa, paste0(this_dir, mon, "_aedes_psa.rds"))
      
      this_aedes_p_psa = this_aedes_pf %>%
        mutate(Month = unique(this_aedes_p$Month)) %>%
        select(x, y, pr_ab, Month) %>%
        # full_join(this_aedes_psa, by = c("longitude"="x", "latitude"="y", "Month", "pr_ab")) %>% ##try1
        # rename(Longitude=longitude, Latitude=latitude)
        rbind(this_aedes_psa) %>%   ##try2
        mutate(species = "aedes")
      
      this_aedes_p_psa_spdf = SpatialPointsDataFrame(
        coords = data.frame(this_aedes_p_psa$x, this_aedes_p_psa$y),
        data = this_aedes_p_psa,
        proj4string = CRS("EPSG:3310")
      ) %>%
        sp::spTransform(crs("+proj=longlat")) %>% 
        as_tibble() %>%
        select(
          # Month, pr_ab, Longitude, Latitude, species ##try1
          Month, x, y, pr_ab, species, Longitude=coords.x1, Latitude=coords.x2
        )
      # 
      this_aedes_p_psa_thinned = spThin::thin(
        loc.data = this_aedes_p_psa_spdf,
        lat.col = "Latitude", ##try1
        long.col = "Longitude",
        spec.col = "species",
        thin.par = 0.27,
        reps = 10,
        write.files = F,
        locs.thinned.list.return = T
      ) %>%
        bind_rows() %>%
        distinct() %>%
        left_join(
          this_aedes_p_psa_spdf %>% select(Longitude, Latitude, Month, pr_ab)
        ) %>%
        SpatialPointsDataFrame(
          data = .,
          coords = cbind(.$Longitude, .$Latitude),
          proj4string = CRS("+proj=longlat +datum=WGS84") ##try1
          # proj4string = CRS("EPSG:3310")
        ) %>%
        spTransform(crs("epsg:3310")) %>% 
        as_tibble() %>% 
        rename(x = coords.x1, y = coords.x2)
      
      saveRDS(this_aedes_p_psa_thinned, paste0(this_dir, mon, "_thinned_p_psa.rds"))
      
      message(paste0(mon, " has finished."))
      message(paste0(mon, ": ", t0 - Sys.time()))
      
      
      #### DELETE ####
      
      # t0 = Sys.time()
      # this_dir = paste0("Data/1_DataProcessing/pseudoabs/geo_env/", mon, "/")
      # this_aedes_psa = readRDS(paste0(this_dir, mon, "_aedes_psa.rds"))
      # # get environmental raster for month
      # this_env_rast = env_stack[[mon]]
      # 
      # # subset presence points to month
      # this_aedes_p = aedes_p %>% 
      #   filter(Month == mon) %>% 
      #   tibble()
      
      
      # this_aedes_p_psa = this_aedes_pf %>%
      #   mutate(Month = unique(this_aedes_p$Month)) %>%
      #   select(x, y, pr_ab, Month) %>%
      #   # full_join(this_aedes_psa, by = c("longitude"="x", "latitude"="y", "Month", "pr_ab")) %>% ##try1
      #   # rename(Longitude=longitude, Latitude=latitude)
      #   rbind(this_aedes_psa) %>%   ##try2
      #   mutate(species = "aedes")
      # 
      # this_aedes_p_psa_spdf = SpatialPointsDataFrame(
      #     # coords = data.frame(this_aedes_p_psa$Longitude, this_aedes_p_psa$Latitude), ##try1
      #     coords = data.frame(this_aedes_p_psa$x, this_aedes_p_psa$y),
      #     data = this_aedes_p_psa,
      #     proj4string = CRS("EPSG:3310")
      #   ) %>%
      #   sp::spTransform(crs("+proj=longlat")) %>% 
      #   as_tibble() %>%
      #   select(
      #     # Month, pr_ab, Longitude, Latitude, species ##try1
      #     Month, x, y, pr_ab, species, Longitude=coords.x1, Latitude=coords.x2
      #   )
      # # 
      # this_aedes_p_psa_thinned = spThin::thin(
      #   loc.data = this_aedes_p_psa_spdf,
      #   lat.col = "Latitude", ##try1
      #   long.col = "Longitude",
      #   # lat.col = "y", ##try2
      #   # long.col = "x",
      #   spec.col = "species",
      #   thin.par = 0.27,
      #   reps = 10,
      #   write.files = F,
      #   locs.thinned.list.return = T
      # ) %>%
      #   bind_rows() %>%
      #   distinct() %>%
      #   # rename(x = Longitude, y = Latitude) %>%
      #   left_join(
      #     this_aedes_p_psa_spdf %>% select(Longitude, Latitude, Month, pr_ab)
      #   ) %>%
      #   SpatialPointsDataFrame(
      #     data = .,
      #     coords = cbind(.$Longitude, .$Latitude),
      #     proj4string = CRS("+proj=longlat +datum=WGS84") ##try1
      #     # proj4string = CRS("EPSG:3310")
      #   ) %>%
      #   spTransform(crs("epsg:3310")) %>% 
      #   as_tibble() %>% 
      #   rename(x = coords.x1, y = coords.x2)
      
    },
    error = function(cond){
      message(paste0(mon, " has failed. See month directory for more details."))

      message(cond)
    }
  )
}



registerDoParallel(detectCores() - 1)

# run in parallel
foreach(m = 1:12)%dopar%{
  psuedoabsenceFunc(m)
}

stopImplicitCluster()


# get all pseudo-absence files and reorder
## As of March 10, only able to get pseudo-absences for April to November
pseudoabs_files = list.files(
  path = "Data/1_DataProcessing/pseudoabs/geo_env/",
  recursive = T, pattern = "psa.rds", full.names = T
)[c(1, 4, 5, 6, 7, 8, 9, 10, 2, 3)]

pseudoabs_records = map(pseudoabs_files, readRDS) %>%
  map(bind_rows) %>%
  map2(.x = ., .y = c(1, 3:11),
       function(x, .y){
         x = x %>%
           mutate(Month = .y)
       }) %>%
  bind_rows()

presence_files = list.files(
  path = "Data/1_DataProcessing/pseudoabs/geo_env/",
  recursive = T, pattern = "filtered.rds", full.names = T
)[c(3:9, 1, 2)]

presence_records = map(presence_files, readRDS) %>% 
  map(bind_rows) %>% 
  map2(.x = ., .y = c(3:11),
       function(x, .y){
         x = x %>%
           mutate(Month = .y)
       }) %>%
  bind_rows()

# combine the presences and pseudo-absences
aedes_p_psa = presence_records %>%
  full_join(pseudoabs_records, by = c("x", "y", "Month", "pr_ab")) %>%
  select(Month, pr_ab, x, y) %>%
  arrange(Month) %>% 
  tibble()

split_p_psa = aedes_p_psa %>% 
  SpatialPointsDataFrame(
    coords = data.frame(aedes_p_psa$x, aedes_p_psa$y), 
    data = .,
    proj4string = CRS("EPSG:3310")
  ) %>% 
  spTransform(CRS("+proj=longlat")) %>% 
  as_tibble() %>% 
  select(
    Month, pr_ab, x, y, Longitude=coords.x1, Latitude=coords.x2
  ) %>% 
  mutate(species = "aedes") %>% 
  group_by(Month) %>% 
  group_split()


t0=Sys.time()
thinned_split_p_psa = for(i in 1:length(split_p_psa)){
  this_month = split_p_psa[[i]] %>% 
    pull(Month) %>% unique()
  print(this_month)
  this = spThin::thin(
    loc.data = split_p_psa[[i]],
    lat.col = "Latitude",
    long.col = "Longitude",
    spec.col = "species",
    thin.par = 0.27,
    reps = 20,
    write.files = F,
    locs.thinned.list.return = T,
    out.dir = paste0(this_dir, "thinned"),
    out.base = paste0(mon, "_pf")
  ) %>% 
    bind_rows() %>% 
    distinct() %>% 
    # mutate(pr_ab = 1) %>% 
    SpatialPointsDataFrame(
      coords = data.frame(.$Longitude, .$Latitude),
      data = .,
      proj4string = CRS("+proj=longlat +datum=WGS84") ## try1
    ) %>% 
    spTransform(CRS("EPSG:3310")) %>% 
    as.data.frame() %>% 
    rename(x = coords.x1, y = coords.x2) %>% 
    select(Longitude, Latitude) %>% 
    left_join(split_p_psa[[i]])
  
  write.csv(this, paste0("Data/1_DataProcessing/pseudoabs/geo_env/", this_month, "/", this_month, "_thinned_p_psa.csv"))
}

Sys.time()-t0

# temp thin split psa --------------------------------------------------------

pseudoabs_records_split = pseudoabs_records %>% 
  group_by(Month) %>% 
  group_split()

thinned_psa = for(i in 1:length(psuedoabs_records_split)){
  this_pseudoabs = psuedoabs_records_split[[i]] %>% 
    mutate(species = "aedes") %>% 
    SpatialPointsDataFrame(
      coords = data.frame(.$x, .$y),
      data = .,
      proj4string = CRS("EPSG:3310")
    ) %>% 
    sp::spTransform(
      CRS("+proj=longlat +datum=WGS84") ##try1
    ) %>% 
    as.data.frame() %>% 
    rename(Longitude = coords.x1, Latitude = coords.x2)
  
  this_month = this_pseudoabs %>% 
    pull(Month) %>% 
    unique()
  
  print(this_month)
  
  this_thin = spThin::thin(
    loc.data = this_pseudoabs,
    lat.col = "Latitude",
    long.col = "Longitude",
    spec.col = "species",
    thin.par = 0.27,
    reps = 10,
    write.files = F,
    locs.thinned.list.return = T
  ) %>% 
    bind_rows() %>% 
    distinct() %>% 
    mutate(pr_ab = 0) %>%
    SpatialPointsDataFrame(
      coords = data.frame(.$Longitude, .$Latitude),
      data = .,
      proj4string = CRS("+proj=longlat +datum=WGS84") ## try1
    ) %>% 
    sp::spTransform(
      CRS("EPSG:3310") ##try1
    ) %>% 
    as.data.frame() %>% 
    rename(x = coords.x1, y = coords.x2) %>% 
    select(x, y, pr_ab) %>% 
    mutate(Month = this_month)
    
  saveRDS(this_thin, paste0("Data/1_DataProcessing/pseudoabs/geo_env/", this_month, "/", this_month, "_thinned_psa.csv"))
}

thinned_psa_files = list.files(
  "Data/1_DataProcessing/pseudoabs/geo_env/", 
  pattern = "thinned_psa.csv", 
  recursive = T, full.names = T
)

thinned_psa_records = map(
  .x = thinned_psa_files,
  .f = readRDS
) %>% bind_rows()

p_psa_files = list.files(
  "Data/1_DataProcessing/pseudoabs/geo_env/", 
  pattern = "thinned_p_psa.rds", 
  recursive = T, full.names = T
)

thinned_p_psa_data = map2(
  .x = p_psa_files,
  .y = c(10, 11, 4:9),
  ~ readRDS(.x) %>% 
    mutate(Month = .y)
)[c(3:8, 1, 2)] %>% 
  bind_rows() 
  

# readr::write_csv(thinned_p_psa_data, "Data/1_DataProcessing/pseudoabs/geo_env/aedes_thinned_presence_psa.csv")
# thinned_split_p_psa = read_csv("Data/1_DataProcessing/pseudoabs/geo_env/aedes_thinned_presence_psa.csv")

# convert to SpatVector for visualization
aedes_p_psa_vec = terra::vect(
  thinned_p_psa_data,
  geom=c("Longitude", "Latitude"),
  crs = crs("epsg:4326")
  # crs="+proj=longlat +datum=WGS84"
  # geom=c("x", "y"),
  # crs = crs("EPSG:3310")
) #%>% 
  #terra::project(CRS("epsg:4326"))

library(ggmap)
register_google("AIzaSyCPfXOob9nGUTPq8tBL7TrZiH1dAHsl4WM")

# base map of CA
ggmap(
  get_map(location = "California, CA", zoom = 6)
) +
  # add SpatVector
  geom_spatvector(
    data = aedes_p_psa_vec,
    aes(color = as.factor(pr_ab)),
    # set inherit = F to not run into error
    inherit = F,
    size = 0.2
  ) +
  # faceted on month
  facet_wrap(~Month)

cali = map_data("state") %>% 
  filter(region=="california") %>%
  select(long, lat) %>%
  as.matrix(ncol = 2) %>%
  list() %>%
  sf::st_polygon() %>%
  as(., "Spatial")
crs(cali) = "epsg:4326"
cali = spTransform(cali, "epsg:3310")

aedes_p_psa_sf = SpatialPointsDataFrame(
  data = thinned_p_psa_data,
  coords = cbind(thinned_p_psa_data$x, thinned_p_psa_data$y),
  proj4string = CRS("epsg:3310")
) %>% 
  st_as_sf()

cali_sf =  st_as_sf(cali, coords = c("long", "lat"), crs = CRS("epsg:4326")) %>% 
  # st_convex_hull() %>% 
  st_transform(CRS("epsg:3310"))

# Transform projection to EPSG:32610
# cali_sf <- st_transform(cali_sf, crs = "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs +type=crs")

  # terra::vect(
  #   .,
  #   geom=c("long", "lat"),
  #   # crs = crs(env_stack, proj=T)
  #   crs = "+proj=longlat"
  # ) %>%
  # as.polygons()

plot(cali_sf)
plot(aedes_p_psa_sf)
  
ggplot() +
  geom_polygon(
    data = cali,
    aes(x = long, y = lat, group = group)
  ) +
  # geom_sf(
  #   data = cali_sf
  # ) +
  # add SpatVector
  geom_point(
    data = thinned_p_psa_data,
    aes(x = x, y = y, color = as.factor(pr_ab)),
    size = 0.1
  ) +
  # faceted on month
  facet_wrap(~Month)




# TRY ---------------------------------------------------------------------

t0 = Sys.time()
this_dir = paste0("Data/1_DataProcessing/pseudoabs/geo_env/", mon, "/")

thinned_p_psa_data

aedes_p_psa_spdf = SpatialPointsDataFrame(
  # coords = data.frame(this_aedes_p_psa$Longitude, this_aedes_p_psa$Latitude), ##try1
  coords = data.frame(thinned_p_psa_data$x, thinned_p_psa_data$y),
  data = thinned_p_psa_data,
  proj4string = CRS("EPSG:3310")
) %>% st_as_sf()


cali2 = maps::map('state', plot=F, fill=T, regions="california") %>%
  st_as_sf() %>%
  # st_transform("+proj=longlat +datum=WGS84")
  st_transform(CRS("EPSG:3310"))


ggplot() +
  geom_sf(
    data = cali2
  ) +
  # add SpatVector
  geom_sf(
    data = aedes_p_psa_spdf %>% st_transform(CRS("EPSG:3310")),# %>% st_transform("+proj=longlat +datum=WGS84"),
    aes(color = as.factor(pr_ab)),
    size = 0.1
  ) +
  # faceted on month
  facet_wrap(~Month)

env_stack_df = env_stack %>% 
  as.data.frame(xy=T) 

env_stack_df_longer = env_stack_df %>% 
  pivot_longer(-c(x, y))

ggplot() + 
  geom_sf(
    data = cali_sf
  ) +
  stars::geom_stars(
    data = stars::st_as_stars(env_stack$tmn2010sep),
    alpha = 0.3
  ) #+ 
  # facet_wrap(name ~ .) +
  # coord_equal()


# Appendix ----------------------------------------------------------------

#' Author: TJ Sipin
#' Date: March 8, 2024
#' Purpose:
#' Get pseudo-absence points for each month

# set.seed(030824)
# setwd("/home/tjsipin/Dengue_CA")
# 
# # Packages ----------------------------------------------------------------
# 
# library(terra)
# library(raster)
# library(purrr)
# library(dplyr)
# library(flexsdm)
# tidymodels::tidymodels_prefer()
# library(foreach)
# library(doParallel)
# library(ggplot2)
# library(tidyterra)
# library(readr)
# 
# # Data --------------------------------------------------------------------
# 
# # path to environmental rasters
# env_raster_path = "Data/1_DataProcessing/rasters/env/"
# env_rasters = map(list.files(env_raster_path, full.names = T), rast)
# env_stack = rast(env_rasters)[[c(5, 4, 8, 1, 9, 7, 6, 2, 12, 11, 10, 3)]]
# 
# # read in presence data
# aedes_p = read.csv("Data/1_DataProcessing/Aedes/aedes_presence_mixedfemales.csv") %>% 
#   select(longitude, latitude, Month) %>% 
#   # set ID and presence values
#   mutate(
#     id = row_number(),
#     pr_ab = 1
#   ) %>% 
#   # turn to spatial points data frame
#   SpatialPointsDataFrame(
#     coords = data.frame(.$longitude, .$latitude),
#     data = .,
#     proj4string = CRS("+proj=longlat")
#   ) %>% 
#   # transfom crs to NAD83
#   spTransform(crs(env_stack)) %>%
#   # revert back to data frame
#   as.data.frame() %>% 
#   # rename to x, y coords
#   rename(
#     x = coords.x1,
#     y = coords.x2
#   ) %>% 
#   # make species column for thinning function
#   mutate(species = "aedes")
# 
# # given month, get spatially thinned presences and write pseudo-absence locations as RDS
# psuedoabsenceFunc = function(mon){
#   tryCatch(
#     {
#       t0 = Sys.time()
#       this_dir = paste0("Data/1_DataProcessing/pseudoabs/geo_env/", mon, "/")
#       dir.create(this_dir)
#       # get environmental raster for month
#       this_env_rast = env_stack[[mon]]
#       
#       # subset presence points to month
#       this_aedes_p = aedes_p %>% 
#         filter(Month == mon) %>% 
#         tibble() %>% 
#         mutate(species = "aedes")
#       
#       # this_aedes_pf = this_aedes_p %>% 
#       #   occfilt_geo(
#       #     data = .,
#       #     x = "x",
#       #     y = "y",
#       #     env_layer = this_env_rast,
#       #     method = c("defined", d = 0.27),
#       #     prj = crs(this_env_rast, proj=T)
#       #   )
#       
#       this_aedes_pf = spThin::thin(
#         loc.data = this_aedes_p,
#         lat.col = "latitude",
#         long.col = "longitude",
#         spec.col = "species",
#         thin.par = 0.27,
#         # reps = this_aedes_p %>% 
#         #   select(x, y) %>% 
#         #   distinct() %>% 
#         #   nrow(),
#         reps = 20,
#         write.files = T,
#         locs.thinned.list.return = T,
#         out.dir = paste0(this_dir, "thinned"),
#         out.base = paste0(mon, "_pf")
#       ) %>% 
#         bind_rows() %>% 
#         distinct() %>% 
#         mutate(pr_ab = 1) %>% 
#         SpatialPointsDataFrame(
#           coords = data.frame(.$Longitude, .$Latitude),
#           data = .,
#           proj4string = CRS("+proj=longlat")
#         ) %>% 
#         spTransform(crs(this_env_rast)) %>% 
#         as.data.frame() %>% 
#         rename(
#           x = coords.x1,
#           y = coords.x2
#         )
#       
#       # print(paste0(mon, ": calibration"))
#       # this_ca <-
#       #   calib_area(
#       #     data = this_aedes_pf,
#       #     x = 'Longitude',
#       #     y = 'Latitude',
#       #     # defined by buffered minimum convex polygon
#       #     method =  c("bmcp", width=15000),
#       #     crs = crs(this_env_rast)
#       #   ) # create a calibration area with 15 km buffer around occurrence points
#       
#       # create partition blocks
#       ### should these parameters be changed?
#       print(paste0(mon, ": occurrence partition"))
#       this_occ_part <- this_aedes_p %>%
#         part_sblock(
#           data = .,
#           env_layer = this_env_rast,
#           pr_ab = "pr_ab",
#           x = "x",
#           y = "y",
#           prop = 0.5,
#           n_part = 3,
#           min_res_mult = 200,
#           max_res_mult = 4000,
#           num_grids = 100
#         )
#       
#       saveRDS(this_occ_part, paste0(this_dir, mon, "_occ_part.rds"))
#       
#       this_part = this_occ_part$part %>% 
#         mutate(species = "aedes")
#       # mutate(species = "aedes") %>% 
#       # SpatialPointsDataFrame(
#       #   data = .,
#       #   coords = data.frame(.$x, .$y)
#       # )
#       
#       this_part_f = this_aedes_pf %>% 
#         select(x, y) %>% 
#         left_join(
#           this_part %>% select(x, y, .part, pr_ab)
#         ) %>% 
#         distinct()
#       
#       print(paste0(mon, ": calibration"))
#       this_ca <-
#         calib_area(
#           data = this_part_f,
#           x = 'x',
#           y = 'y',
#           # defined by buffered minimum convex polygon
#           method =  c("bmcp", width=15000),
#           crs = crs(this_env_rast)
#         ) # create a calibration area with 15 km buffer around occurrence points
#       
#       # this_part_f =
#       # spThin::thin(
#       #   loc.data = this_part,
#       #   lat.col = "y",
#       #   long.col = "x",
#       #   spec.col = "species",
#       #   thin.par = 0.27,
#       #   # reps = this_aedes_p %>%
#       #   #   select(x, y) %>%
#       #   #   distinct() %>%
#       #   #   nrow(),
#       #   reps = 20,
#       #   write.files = F,
#       #   locs.thinned.list.return = T,
#       #   out.dir = paste0(this_dir, "thinned"),
#       #   out.base = paste0(mon, "_pf")
#       # ) %>%
#       # bind_rows() %>%
#       # distinct() %>%
#       # mutate(pr_ab = 1) %>%
#       # rename(
#       #   x = Longitude,
#       #   y = Latitude
#       # )
#       # SpatialPointsDataFrame(
#       #   coords = data.frame(.$x, .$y),
#       #   data = .,
#       #   # proj4string = CRS(this_env_rast)
#       # ) %>% 
#       # # spTransform(crs(this_env_rast)) %>%
#       # as.data.frame() %>% 
#       # select(x, y)
#       
#       print(paste0(mon, ": block_layer"))
#       this_block_layer = get_block(
#         env_layer = this_env_rast, best_grid = this_occ_part$grid
#       )
#       
#       saveRDS(this_block_layer, paste0(this_dir, mon, "_block_layer.rds"))
#       
#       
#       this_aedes_psa = lapply(sort(unique(this_part_f$.part)), function(x){
#         sample_pseudoabs(
#           data = this_part_f,
#           x = "x",
#           y = "y",
#           n = sum(this_part_f$.part == x),
#           method = c("geo_env_const", width = "10000", env = this_env_rast), #10km buffer, is this ok?
#           rlayer = this_block_layer,
#           # maskval = x,
#           calibarea = this_ca
#         )
#       }) %>% 
#         bind_rows() 
#       
#       saveRDS(this_part_f, paste0(this_dir, mon, "_partition_filtered.rds"))
#       saveRDS(this_aedes_psa, paste0(this_dir, mon, "_aedes_psa.rds"))
#       
#       message(paste0(mon, " has finished."))
#       message(paste0(mon, ": ", t0 - Sys.time()))
#     },
#     error = function(cond){
#       message(paste0(mon, " has failed. See month directory for more details."))
#       
#       message(cond)
#     }
#   )
# }
# 
# 
# 
# registerDoParallel(detectCores() - 1)
# 
# # run in parallel
# foreach(m = 1:12)%dopar%{
#   psuedoabsenceFunc(m)
# }
# 
# stopImplicitCluster()
# 
# 
# # get all pseudo-absence files and reorder
# ## As of March 10, only able to get pseudo-absences for April to November
# pseudoabs_files = list.files(
#   path = "Data/1_DataProcessing/pseudoabs/geo_env/",
#   recursive = T, pattern = "psa.rds", full.names = T
# )[c(1, 4, 5, 6, 7, 8, 9, 10, 2, 3)]
# 
# psuedoabs_records = map(pseudoabs_files, readRDS) %>%
#   map(bind_rows) %>%
#   map2(.x = ., .y = c(1, 3:11),
#        function(x, .y){
#          x = x %>%
#            mutate(Month = .y)
#        }) %>%
#   bind_rows()
# 
# presence_files = list.files(
#   path = "Data/1_DataProcessing/pseudoabs/geo_env/",
#   recursive = T, pattern = "filtered.rds", full.names = T
# )[c(3:9, 1, 2)]
# 
# presence_records = map(presence_files, readRDS) %>% 
#   map(bind_rows) %>% 
#   map2(.x = ., .y = c(3:11),
#        function(x, .y){
#          x = x %>%
#            mutate(Month = .y)
#        }) %>%
#   bind_rows()
# 
# # combine the presences and pseudo-absences
# aedes_p_psa = presence_records %>%
#   full_join(psuedoabs_records) %>%
#   select(Month, pr_ab, x, y) %>%
#   arrange(Month) %>% 
#   tibble()
# 
# readr::write_csv(aedes_p_psa, "Data/1_DataProcessing/pseudoabs/geo_env/aedes_thinned_presence_psa.csv")
# 
# # convert to SpatVector for visualization
# aedes_p_psa_vec = terra::vect(
#   aedes_p_psa,
#   geom=c("x", "y"),
#   crs = crs(env_stack, proj=T)
# ) %>%
#   #reproject to long lat to overlay onto base map
#   terra::project("+proj=longlat")
# 
# library(ggmap)
# register_google("AIzaSyCPfXOob9nGUTPq8tBL7TrZiH1dAHsl4WM")
# 
# # base map of CA
# ggmap(
#   get_map(location = "California, CA", zoom = 6)
# ) +
#   # add SpatVector
#   geom_spatvector(
#     data = aedes_p_psa_vec,
#     aes(color = as.factor(pr_ab)),
#     # set inherit = F to not run into error
#     inherit = F,
#     size = 0.2
#   ) +
#   # faceted on month
#   facet_wrap(~Month)
# 
# 
# 
# # Appendix ----------------------------------------------------------------
# 
# # # produce all available pseudoabsences
# psuedoabsenceFunc_v2 = function(mon){
#   tryCatch(
#     {
#       this_dir = paste0("Data/1_DataProcessing/pseudoabs/", mon, "/")
#       # get environmental raster for month
#       this_env_rast = env_stack[[mon]]
#       
#       # subset presence points to month
#       this_aedes_p = aedes_p %>%
#         filter(Month == mon) %>%
#         tibble()
#       
#       this_ca <-
#         calib_area(
#           data = this_aedes_p,
#           x = 'x',
#           y = 'y',
#           method =  c("bmcp", width=15000),
#           crs = crs(this_env_rast)
#         ) # create a calibration area with 15 km buffer around occurrence points
#       
#       print(paste0(mon, ": this_occ_part"))
#       
#       this_occ_part = readRDS(paste0(this_dir, mon, "_occ_part.rds"))
#       this_aedes_p_part = this_occ_part$part
#       this_block_layer = readRDS(paste0(this_dir, mon, "_block_layer.rds"))
#       
#       print(paste0(mon, ": this_aedes_psa start"))
#       t0 = Sys.time()
#       this_aedes_psa = lapply(sort(unique(this_aedes_p_part$.part)), function(x){
#         t0 = Sys.time()
#         sample_pseudoabs(
#           data = this_aedes_p_part,
#           x = "x",
#           y = "y",
#           n = sum(this_aedes_p_part$.part == x),
#           method = c("geo_env_const", width = "10000", env = this_env_rast), #10km buffer, is this ok?
#           rlayer = this_block_layer,
#           maskval = x,
#           calibarea = this_ca
#         )
#       })
#       t0 - Sys.time()
#       print(paste0(mon, ": this_aedes_psa finish"))
#       
#       saveRDS(this_aedes_psa, paste0(this_dir, mon, "_aedes_psa.rds"))
#       
#       message(paste0(mon, " has finished."))
#       return(NULL)
#     },
#     error = function(cond){
#       message(cond)
#       return(NULL)
#     },
#     warning = function(cond){
#       message(cond)
#       return(NULL)
#     }
#   )
# }
# 
# # for(m in c(7, 11)){
# #   psuedoabsenceFunc(m)
# # }
