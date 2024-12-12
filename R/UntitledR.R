

raw_vect = raw %>% 
    filter(coordinate_precision!="Unknown") %>% 
    filter(calculated_state == "California") %>% 
    mutate(
        collection_date = mdy(collection_date), 
        Year = year(collection_date),
        Month = month(collection_date),
        woy = week(collection_date)
    ) %>% 
    select(collection_id, code, Year, Month, woy, species, longitude, latitude) %>% 
    vect(geom=c("longitude", "latitude"), crs = crs("EPSG:4326")) %>% 
    terra::project(crs("EPSG:3310"))

ggplot() + 
    geom_spatvector(data = CA) +
    geom_spatvector(
        data = raw_vect %>% 
            filter(species == "Aedes albopictus"), 
        aes(color = species), 
        size = 0.5, alpha = 0.5
    )
