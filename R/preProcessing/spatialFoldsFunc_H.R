library(sf) # to read in shapefiles
library(geosphere) # for distance functions
library(terra) # to read in rasters
# library(rgeos) # TODO
library(dplyr) # always
library(tidyr) # always
library(cluster) # TODO
library(lubridate) # for date functions
library(purrr) # for map function
library(sp) # for spatial points data frame
library(readr)
library(rsample)
library(spatialsample)
tidymodels::tidymodels_prefer()

# set seed
set.seed(123123)

# spatial folds function to join with the training data
getSpatialFolds = function(training){
    # Feature space then spatial ----------------------------------------------
    
    avg_sil_pt1 = factoextra::fviz_nbclust(
        training %>% 
            select(-c(Month, species, pr_ab)), 
        FUN = kmeans, method = "silhouette", k.max = min(round(nrow(training)*0.5), 10)
    )
    
    dist_training = proxy::dist(
        training %>% 
            select(-c(Month, species, pr_ab))
    ) %>% 
        hclust()
    
    # cut tree into opt_num_clusts clusters
    location_cutree <- cutree(dist_training, k = which.max(avg_sil_pt1$data$y[1:length(avg_sil_pt1$data$y)]))
    
    training_folds = training %>% 
        cbind(init_fold = location_cutree) %>% 
        tibble()
    
    getObs = function(thisFold){
        t1 = Sys.time()
        print(thisFold)
        ids = training_folds %>% 
            filter(init_fold == thisFold) %>% 
            pull(id)
        
        obs = training_folds %>% 
            filter(init_fold == thisFold)
        
        this_fold_spdf = SpatialPointsDataFrame(
            coords = obs[, c("Longitude", "Latitude")],
            data = obs,
            proj4string = CRS("WGS84")
        )
        
        if(length(this_fold_spdf) <= 2){
            res = data.frame(
                id = this_fold_spdf$id,
                block_fold = as.character(thisFold)
            )
            
            return(res)
        }
        
        mdistAllID <- geosphere::distm(this_fold_spdf)
        
        # cluster all points using a hierarchical clustering approach
        hcAllID <- hclust(as.dist(mdistAllID), method="complete")
        
        avg_sil = factoextra::fviz_nbclust(obs, FUN = factoextra::hcut, method = "silhouette", k.max = 10)
        
        
        # opt_num_clusts = which.max(avg_sil$data$y[3:length(avg_sil$data$y)])+2
        opt_num_clusts = which.max(avg_sil$data$y[1:length(avg_sil$data$y)])
        
        # cut tree into opt_num_clusts clusters
        this_fold_spdf$block <- cutree(hcAllID, k = opt_num_clusts)
        
        res = data.frame(
            id = ids,
            block_fold = paste(thisFold, this_fold_spdf$block, sep = "_")
        )
        
        res
        
    }
    
    t0 = Sys.time()
    id_block_folds = map(training_folds$init_fold %>% unique(), getObs) %>% 
        bind_rows()
    Sys.time() - t0
    
    block_fold_fold_dict = data.frame(
        block_fold = sort(unique(id_block_folds$block_fold)),
        fold = 1:length(sort(unique(id_block_folds$block_fold)))
    )
    
    training_folds_v2 = training_folds %>%
        left_join(id_block_folds) %>% 
        left_join(block_fold_fold_dict) %>% 
        select(-block_fold, -init_fold) %>% 
        # full_join(clust_spdf %>% st_as_sf()) %>% 
        # select(clust, block, fold, geometry, x, y) %>% 
        filter(!is.na(fold)) %>% 
        st_as_sf(coords = c("round_x", "round_y"), crs = "EPSG:3310") %>%
        # st_transform(crs = "EPSG:3310") %>%
        # Get number of points in cluster to appropriately handle center
        dplyr::group_by(fold) %>%
        mutate(
            pr_ab = as.factor(pr_ab),
            numPoints = n()
        ) %>% 
        st_as_sf()
    
    # Cluster with 1 location
    points <- training_folds_v2 %>%
        filter(numPoints == 1) %>%
        dplyr::summarise()
    
    
    # Clusters with 2 locations
    lines <- training_folds_v2 %>%
        filter(numPoints == 2) %>%
        dplyr::summarise() %>%
        st_cast("LINESTRING")
    
    # CLusters with 3 or more
    polys <- training_folds_v2 %>%
        filter(numPoints > 2) %>%
        dplyr::summarise() %>%
        # Need to cast to polygon and create convex hull around points
        st_cast("POLYGON") %>%
        st_convex_hull()
    
    
    # Get lines and polys to make circumscribing circle
    allFeats <- points %>% 
        bind_rows(lines) %>%
        bind_rows(polys)
    
    plot(allFeats)
    
    # ggplot(training_folds_v2) + #, aes(x = longitude, y = latitude)) + 
    #     # geom_polygon(aes(
    #     #   color = fold %>% as.factor)) +
    #     geom_sf(
    #         aes(color = fold %>% as.factor())) 
    
    mapview::mapview(training_folds_v2, z = "fold")
    
    # join together the spatial folds with the training to output
    training %>% 
        left_join(training_folds_v2) 
}

getSpatialFolds_v2 = function(training){
    sf = st_as_sf(
        training, 
        coords = c("x", "y"),
        crs = 3310
    )
    
    
    spatial_clusters = spatial_clustering_cv(sf, v = 4)
    spatial_blocks = spatial_block_cv(sf, v = 4)
    
    return(
        list(
            spatial_clusters = spatial_clusters,
            spatial_blocks = spatial_blocks
        )
    )
}

