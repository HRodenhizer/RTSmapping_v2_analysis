##### Fix the NaturalEarth dataset
### align_coords
# This function takes coordinates within 0.115 degrees of 180/-180 and replaces them with 180/-180 depending on the sign of the original coordinate.
# This allows polygons on either side of the 180th meridian to later be unioned to remove pesky lines in the middle of countries/landmasses.
align_coords <- function(x, dist) {
  case_when(
    -180 - x > dist * -1 & sign(x) == -1 ~ -180,
    180 - x < dist & sign(x) == 1 ~ 180,
    TRUE ~ x
  )
}

### fix_geometry_at_origin
# This function is designed to remove small, purposefully offsets around the 180th meridian in the natural earth datasets.
# The naturalearth datasets must first be converted to a projected coordinate system with units of m.
# It takes coordinates smaller than 2 m and replaces them with 0.
# This allows polygons on either side of the 180th meridian to later be unioned to remove pesky lines in the middle of countries/landmasses.
fix_geometry_at_origin <- function(geometry, dist) {
  # Replace incorrect values
  new_values <- geometry[[1]][[1]][, 1]
  new_values <- align_coords(new_values, dist)
  geometry[[1]][[1]][, 1] <- new_values

  return(geometry)
}

##### Lat and Lon Lines for maps
### long_lines_from_coords
long_lines_from_coords <- function(lon, lat_min, lat_max) {
  line <- expand_grid(lon = lon, lat = c(lat_min, lat_max)) %>%
    st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326") %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("LINESTRING") %>%
    smoothr::densify(max_distance = 1) %>%
    mutate(label = lon)

  return(line)
}

### lat_lines_from_coords
lat_lines_from_coords <- function(lat, lon_min, lon_max) {
  line <- expand_grid(lon = c(lon_min, lon_max), lat = lat) %>%
    st_as_sf(coords = c("lon", "lat"), crs = "EPSG:4326") %>%
    summarise(geometry = st_combine(geometry)) %>%
    st_cast("LINESTRING") %>%
    smoothr::densify(max_distance = 1) %>%
    mutate(label = lat)

  return(line)
}
