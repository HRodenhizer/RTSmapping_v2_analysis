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

##### Fix invalid geometries in mountain glaciers dataset
# clean_mp
# this function removes multipolygon parts with 2 points or less
clean_mp <- function(geom) {
  parts <- lapply(unclass(geom), function(poly) {
    if (nrow(poly[[1]]) < 3) {
      return(NULL)
    } # collapsed exterior: drop the whole part
    Filter(function(ring) nrow(ring) >= 3, poly) # collapsed holes: drop just the ring
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) {
    return(st_multipolygon())
  }
  st_multipolygon(parts)
}

# fix_nested_shell
# this function takes the symmetric difference of multipolygon parts that fall within the shell to turn them into holes
fix_nested_shell <- function(g, crs) {
  p <- st_cast(st_sfc(g, crs = crs), "POLYGON") |> as.list()
  Reduce(
    function(a, b) {
      st_sym_difference(st_sfc(a, crs = crs), st_sfc(b, crs = crs))[[1]]
    },
    p
  )
}

# holes_outside
# this function determines if there are multipolygon parts that are defined in the wrong shell
holes_outside <- function(g) {
  sapply(unclass(g), function(poly) {
    if (length(poly) < 2) {
      return(0L)
    }
    shell <- st_polygon(poly[1])
    sum(
      !vapply(
        poly[-1],
        function(h) {
          st_within(st_polygon(list(h)), shell, sparse = FALSE)[1, 1]
        },
        logical(1)
      )
    )
  })
}

# rebuild_mp
# this function reassigns holes in multipolygons to the correct shell
rebuild_mp <- function(g) {
  rings <- unlist(unclass(g), recursive = FALSE) # flat list of ring matrices
  polys <- st_sfc(lapply(rings, \(r) st_polygon(list(r)))) # no CRS: planar tests
  cont <- st_contains_properly(polys, polys, sparse = FALSE) # cont[i, j]: i contains j
  depth <- colSums(cont)
  area <- as.numeric(st_area(polys))

  ext <- which(depth %% 2 == 0)
  hol <- which(depth %% 2 == 1)

  parent <- vapply(
    hol,
    function(h) {
      cand <- ext[cont[ext, h]]
      cand[which.min(area[cand])]
    },
    integer(1)
  )

  st_multipolygon(lapply(ext, \(e) c(rings[e], rings[hol[parent == e]])))
}

# map_hex_grid
map_hex_grid <- function(
  data,
  domain,
  variable,
  color_scale = scale_color_viridis,
  fill_scale = scale_fill_viridis,
  label,
  transformation
) {
  ggplot() +
    geom_sf(data = world_north) +
    geom_sf(
      data = data,
      aes(color = .data[[variable]], fill = .data[[variable]])
    ) +
    color_scale(
      name = label,
      trans = transformation,
      na.value = "transparent"
    ) +
    fill_scale(
      name = label,
      trans = transformation,
      na.value = "transparent"
    ) +
    geom_sf(
      data = world_north,
      fill = "transparent"
    ) +
    geom_sf(
      data = domain,
      color = "black",
      fill = "transparent",
      linewidth = 0.5
    ) +
    geom_sf(
      data = crop_poly,
      color = "black",
      fill = "transparent"
    ) +
    # coord_sf(expand = FALSE) +
    theme_void()
}
