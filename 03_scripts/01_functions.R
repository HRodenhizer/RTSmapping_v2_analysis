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

##### Mapping
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

##### Functions for Bivariate Mapping
# classify_counts
# This classifies a continous variable for bivariate mapping
classify_counts <- function(count_data, breaks) {
  new_count_data <- case_when(
    count_data <= breaks[1] ~ 1,
    count_data <= breaks[2] ~ 2,
    count_data <= breaks[3] ~ 3,
    count_data <= breaks[4] ~ 4,
    count_data > breaks[4] ~ 5
  )
}

# custom_bi_legend
custom_bi_legend <- function(
  data,
  var1,
  var2,
  bi_dim,
  xbreaks,
  ybreaks,
  color_palette,
  xlab,
  ylab
) {
  bi_legend(
    pal = color_palette,
    dim = bi_dim,
    xlab = xlab,
    ylab = ylab,
    # size = 5,
    breaks = list(
      "bi_x" = c(
        xbreaks,
        round(max(data |> pull(!!var1), na.rm = TRUE))
      ),
      "bi_y" = c(
        ybreaks,
        round(max(data |> pull(!!var2), na.rm = TRUE))
      )
    )
  ) +
    coord_fixed(
      xlim = c(0.5, 5.5),
      ylim = c(0.5, 5.5),
      clip = "off"
    ) +
    theme(
      plot.margin = margin(6, 5, 6, 10),
      plot.background = element_rect(
        fill = fill_alpha("white", 0)
      )
    )
}

# map_hex_grid_bi
map_hex_grid_bi <- function(
  data,
  bi_var,
  var1,
  var2,
  domain,
  bi_dim,
  xbreaks,
  ybreaks,
  color_palette,
  legend = TRUE,
  legend_loc,
  xlab,
  ylab
) {
  bi_map <- ggplot() +
    geom_sf(data = world_north) +
    geom_sf(
      data = data,
      aes(fill = .data[[bi_var]]),
      color = "transparent",
      alpha = 1
    ) +
    bi_scale_fill(
      pal = color_palette,
      dim = bi_dim,
      guide = "none"
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

  if (legend == TRUE) {
    bi_legend <- custom_bi_legend(
      data,
      var1,
      var2,
      bi_dim,
      xbreaks,
      ybreaks,
      color_palette,
      xlab,
      ylab
    )

    bi_map <- bi_map +
      inset_element(
        bi_legend,
        left = legend_loc["left"],
        bottom = legend_loc["bottom"],
        right = legend_loc["right"],
        top = legend_loc["top"]
      )
  }

  bi_map
}

##### RTS Region Smoothing
# These functions are designed to smooth the edges of polygons derived from
# raster datasets such that shared boundaries remain shared (eliminating gaps
# or overlap after smoothing).

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
# Area (numeric) of a single sfg
area_of <- function(g) as.numeric(st_area(st_sfc(g)))

# Interior rings (holes) of the dissolved footprint of an sfc, as an sfc
get_holes <- function(g) {
  crs <- st_crs(g)
  u <- st_cast(st_union(g), "POLYGON")
  holes <- list()
  for (k in seq_along(u)) {
    for (r in unclass(u[[k]])[-1]) {
      holes[[length(holes) + 1]] <- st_polygon(list(r))
    }
  }
  st_sfc(holes, crs = crs)
}

# Smooth one arc's coordinate matrix (Gaussian-weighted moving average).
# Open arcs: endpoints pinned, window tapers near the ends.
# circular = TRUE: closed ring with no junctions, smoothed with wrap-around.
smooth_coords <- function(m, k, circular) {
  kern <- function(kk) {
    j <- -kk:kk
    w <- dnorm(j, sd = max(kk, 1) / 2)
    w / sum(w)
  }

  if (circular) {
    m <- m[-nrow(m), , drop = FALSE]
    n <- nrow(m)
    k <- min(k, (n - 1) %/% 2)
    if (k < 1) {
      return(rbind(m, m[1, ]))
    }
    out <- m
    for (i in seq_len(n)) {
      idx <- ((i - k):(i + k) - 1) %% n + 1
      out[i, ] <- colSums(m[idx, , drop = FALSE] * kern(k))
    }
    return(rbind(out, out[1, ]))
  }

  n <- nrow(m)
  out <- m
  if (n < 3) {
    return(out)
  }
  for (i in 2:(n - 1)) {
    kk <- min(k, i - 1, n - i)
    idx <- (i - kk):(i + kk)
    out[i, ] <- colSums(m[idx, , drop = FALSE] * kern(kk))
  }
  out
}

# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------
# x            : sf polygons in a projected CRS (holes are allowed)
# window       : half-width of smoothing window, in vertices
# dens         : max segment length (CRS units) used to segmentize arcs before
#                smoothing; ~ half to one raster cell size is a good start.
# iter         : number of smoothing passes
# smooth_outer : TRUE = smooth the dataset's outer perimeter too;
#                FALSE = leave outer arcs exactly as they are
# min_cover    : a face is assigned to an original polygon only if at least
#                this fraction of its area lies inside the original footprint.
#                Faces below this are treated as holes (if they sit in an
#                original hole; reported by message) or as unexpected
#                (reported by warning).
smooth_shared <- function(
  x,
  window = 5,
  dens = NULL,
  iter = 1,
  smooth_outer = TRUE,
  min_cover = 0.5
) {
  stopifnot(inherits(x, "sf"))
  crs <- st_crs(x)
  if (is.na(crs)) {
    stop("`x` has no CRS. Assign a projected CRS first.", call. = FALSE)
  }
  if (isTRUE(st_is_longlat(x))) {
    stop(
      "`x` has a geographic (lon/lat) CRS. Reproject to a projected CRS ",
      "with st_transform() first.",
      call. = FALSE
    )
  }

  area_floor <- 1e-6 # CRS units^2; overlaps smaller than this count as none

  x <- x |> st_make_valid() |> mutate(.id = row_number())
  geom_x <- st_geometry(x)
  holes_orig <- get_holes(geom_x)

  # 1. Unique arcs: node all boundaries, merge only at degree-2 nodes
  arcs <- geom_x |>
    st_boundary() |>
    st_union() |>
    st_line_merge() |>
    st_cast("LINESTRING")

  # Flag outer arcs: midpoint touches fewer than 2 original polygons.
  # (Hole boundaries are also "outer" in this sense: they border one polygon.)
  mid <- st_line_sample(arcs, sample = 0.5) |> st_cast("POINT")
  n_touch <- lengths(st_is_within_distance(mid, geom_x, dist = 1e-6))
  is_outer <- n_touch < 2

  to_smooth <- if (smooth_outer) rep(TRUE, length(arcs)) else !is_outer

  # Arcs touching another arc have junctions -> endpoints must be pinned
  touches_other <- lengths(st_intersects(arcs)) > 1

  # 2 + 3. Smooth each arc once, endpoints pinned
  sm <- lapply(seq_along(arcs), function(i) {
    g <- arcs[[i]]
    if (!to_smooth[i]) {
      return(g)
    }
    if (!is.null(dens)) {
      g <- st_segmentize(g, dens)
    }
    m <- st_coordinates(g)[, 1:2, drop = FALSE]
    closed <- all(m[1, ] == m[nrow(m), ])
    circular <- closed && !touches_other[i]
    for (p in seq_len(iter)) {
      m <- smooth_coords(m, window, circular)
    }
    st_linestring(m)
  }) |>
    st_sfc(crs = crs)

  # 4. Rebuild faces from the smoothed arcs
  faces <- sm |>
    st_union() |>
    st_polygonize() |>
    st_collection_extract("POLYGON") |>
    st_cast("POLYGON")

  # Overlap of every face with every original polygon
  faces_sf <- st_sf(.face = seq_along(faces), geometry = faces)
  x_id <- x[".id"]
  st_agr(faces_sf) <- "constant"
  st_agr(x_id) <- "constant"

  inter <- st_intersection(faces_sf, x_id)
  inter$.a <- as.numeric(st_area(inter))
  best <- inter |>
    st_drop_geometry() |>
    filter(.a > area_floor) |>
    group_by(.face) |>
    summarise(covered = sum(.a), .id = .id[which.max(.a)], .groups = "drop")

  # A face is assigned only if enough of it lies inside the original footprint
  face_area <- as.numeric(st_area(faces))
  covered <- best$covered[match(seq_along(faces), best$.face)]
  covered[is.na(covered)] <- 0
  keep <- covered / face_area >= min_cover
  unmatched <- which(!keep)

  # Of the unmatched faces, which sit inside original holes (expected)?
  in_hole <- rep(FALSE, length(faces))
  if (length(unmatched) > 0 && length(holes_orig) > 0) {
    hu <- st_union(holes_orig)[[1]]
    frac <- vapply(
      unmatched,
      function(k) {
        area_of(st_intersection(faces[[k]], hu)) / face_area[k]
      },
      numeric(1)
    )
    in_hole[unmatched] <- frac >= 0.5
  }
  hole_faces <- unmatched[in_hole[unmatched]]
  dropped <- setdiff(unmatched, hole_faces)

  if (length(hole_faces) > 0) {
    message(sprintf(
      "%d face(s) (total area %.4g CRS units^2) lie within original holes and were left as holes.",
      length(hole_faces),
      sum(face_area[hole_faces])
    ))
  }
  if (length(dropped) > 0) {
    warning(
      sprintf(
        "%d face(s) lie outside the original footprint and outside any original hole, and were dropped (total area %.4g CRS units^2). This leaves gaps. Check `attr(result, 'dropped_faces')`.",
        length(dropped),
        sum(face_area[dropped])
      ),
      call. = FALSE
    )
  }

  out <- faces_sf |>
    inner_join(
      best[best$.face %in% which(keep), c(".face", ".id")],
      by = ".face"
    ) |>
    group_by(.id) |>
    summarise(.groups = "drop") # dissolve faces to one feature per row

  lost <- setdiff(x$.id, out$.id)
  if (length(lost) > 0) {
    warning(
      sprintf(
        "%d original feature(s) received no face and are missing from the output (row numbers: %s).",
        length(lost),
        paste(lost, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- out |>
    left_join(st_drop_geometry(x), by = ".id") |>
    arrange(.id) |>
    select(-.id)

  attr(out, "dropped_faces") <- faces[dropped] # unexpected only
  attr(out, "hole_faces") <- faces[hole_faces] # expected holes
  out
}

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
# original, smoothed : sf polygon objects
# id          : column identifying features in both objects. If NULL, rows are
#               matched by position (requires equal row counts).
# area_tol    : areas below this (CRS units^2) are treated as numerical noise
# buffer      : distance (CRS units) around each original polygon used for the
#               misassignment check. NULL = 5% of sqrt(median original area).
# outside_tol : flag features with more than this % of area outside the
#               buffered original
smooth_diagnostics <- function(
  original,
  smoothed,
  id = NULL,
  area_tol = 1e-3,
  buffer = NULL,
  outside_tol = 1
) {
  stopifnot(inherits(original, "sf"), inherits(smoothed, "sf"))
  dropped <- attr(smoothed, "dropped_faces")

  # --- align features -------------------------------------------------------
  if (is.null(id)) {
    if (nrow(original) != nrow(smoothed)) {
      stop(
        "Row counts differ; supply `id` (a column present in both objects).",
        call. = FALSE
      )
    }
    key_o <- key_s <- seq_len(nrow(original))
  } else {
    for (nm in c("original", "smoothed")) {
      if (!id %in% names(get(nm))) {
        stop(sprintf("Column '%s' not found in `%s`.", id, nm), call. = FALSE)
      }
    }
    key_o <- original[[id]]
    key_s <- smoothed[[id]]
    if (anyDuplicated(key_o) || anyDuplicated(key_s)) {
      stop(
        sprintf("Column '%s' must be unique in both objects.", id),
        call. = FALSE
      )
    }
  }
  common <- intersect(key_o, key_s)
  if (length(common) == 0) {
    stop(
      sprintf(
        "No matching feature ids between `original` and `smoothed`. First original ids: %s; first smoothed ids: %s",
        paste(head(key_o, 3), collapse = ", "),
        paste(head(key_s, 3), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  missing <- setdiff(key_o, key_s)
  o <- st_geometry(original)[match(common, key_o)]
  s <- st_geometry(smoothed)[match(common, key_s)]

  # --- validity -------------------------------------------------------------
  n_invalid <- sum(!st_is_valid(s))

  # --- overlaps -------------------------------------------------------------
  su <- st_union(s)[[1]]
  total_area <- sum(as.numeric(st_area(s)))
  union_area <- area_of(su)
  ov <- st_relate(s, s, pattern = "T********", sparse = TRUE) # interiors meet
  ov <- data.frame(
    i = rep(seq_along(ov), lengths(ov)),
    j = as.integer(unlist(ov))
  )
  ov <- ov[ov$i < ov$j, , drop = FALSE]
  ov$area <- vapply(
    seq_len(nrow(ov)),
    function(k) area_of(st_intersection(s[[ov$i[k]]], s[[ov$j[k]]])),
    numeric(1)
  )
  ov <- ov[ov$area > area_tol, , drop = FALSE]
  overlaps <- data.frame(
    id_i = common[ov$i],
    id_j = common[ov$j],
    area = ov$area
  )

  # --- holes: unexpected gaps and filled holes ------------------------------
  holes_o <- get_holes(o)
  holes_s <- get_holes(s)
  hu_o <- if (length(holes_o)) st_union(holes_o)[[1]] else NULL

  hs_area <- vapply(
    seq_along(holes_s),
    function(k) area_of(holes_s[[k]]),
    numeric(1)
  )
  hs_in_orig <- vapply(
    seq_along(holes_s),
    function(k) {
      if (is.null(hu_o)) {
        return(0)
      }
      area_of(st_intersection(holes_s[[k]], hu_o)) / hs_area[k]
    },
    numeric(1)
  )
  gap_area <- hs_area[hs_in_orig < 0.5 & hs_area > area_tol] # not in an original hole

  ho_area <- vapply(
    seq_along(holes_o),
    function(k) area_of(holes_o[[k]]),
    numeric(1)
  )
  ho_uncovered <- vapply(
    seq_along(holes_o),
    function(k) {
      area_of(st_difference(holes_o[[k]], su)) / ho_area[k]
    },
    numeric(1)
  )
  holes_filled <- data.frame(area = ho_area[ho_uncovered < 0.5])

  # --- dropped faces (unexpected only) --------------------------------------
  n_dropped <- length(dropped)
  dropped_area <- if (n_dropped) sum(as.numeric(st_area(dropped))) else 0

  # --- parts per feature ----------------------------------------------------
  n_parts <- function(g) {
    vapply(
      g,
      function(p) if (inherits(p, "MULTIPOLYGON")) length(p) else 1L,
      integer(1)
    )
  }
  p_o <- n_parts(o)
  p_s <- n_parts(s)

  # --- adjacency (shared boundary of positive length) -----------------------
  adj_pairs <- function(g) {
    r <- st_relate(g, g, pattern = "F***1****", sparse = TRUE)
    d <- data.frame(
      i = rep(seq_along(r), lengths(r)),
      j = as.integer(unlist(r))
    )
    d <- d[d$i < d$j, , drop = FALSE]
    paste(d$i, d$j, sep = "-")
  }
  decode <- function(p) {
    if (!length(p)) {
      return(data.frame(id_i = common[0], id_j = common[0]))
    }
    ij <- do.call(rbind, strsplit(p, "-"))
    data.frame(
      id_i = common[as.integer(ij[, 1])],
      id_j = common[as.integer(ij[, 2])]
    )
  }
  adj_o <- adj_pairs(o)
  adj_s <- adj_pairs(s)
  lost_adj <- decode(setdiff(adj_o, adj_s))
  gained_adj <- decode(setdiff(adj_s, adj_o))

  # --- per-feature change and misassignment check ---------------------------
  a_o <- as.numeric(st_area(o))
  a_s <- as.numeric(st_area(s))
  if (is.null(buffer)) {
    buffer <- 0.05 * sqrt(median(a_o))
    message(sprintf("`buffer` not supplied; using %.4g CRS units.", buffer))
  }
  buffer <- as.numeric(buffer)
  outside <- vapply(
    seq_along(s),
    function(k) {
      area_of(st_difference(s[[k]], st_buffer(o[[k]], buffer)))
    },
    numeric(1)
  )

  per_feature <- data.frame(
    id = common,
    area_orig = a_o,
    area_new = a_s,
    area_change_pct = 100 * (a_s - a_o) / a_o,
    parts_orig = p_o,
    parts_new = p_s,
    pct_outside = 100 * outside / a_s
  )
  extra_parts <- per_feature[per_feature$parts_new > per_feature$parts_orig, ]
  far_off <- per_feature[per_feature$pct_outside > outside_tol, ]

  # --- summary --------------------------------------------------------------
  # Helpers guarantee exactly one value per check; a check whose input came
  # back empty shows NA / "n/a" instead of breaking data.frame().
  cnt <- function(v) if (length(v) == 1) as.numeric(v) else NA_real_
  fmt <- function(f, ...) {
    r <- tryCatch(sprintf(f, ...), error = function(e) character(0))
    if (length(r) == 1) r else "n/a"
  }

  summary <- data.frame(
    check = c(
      "Invalid geometries",
      "Overlapping pairs",
      "Unexpected gaps (holes not in original)",
      "Original holes filled",
      "Dropped faces (outside footprint and holes)",
      "Features missing from output",
      "Features with more parts",
      "Lost adjacencies",
      "Gained adjacencies",
      sprintf("Features >%g%% outside buffered original", outside_tol)
    ),
    n_issues = c(
      cnt(n_invalid),
      cnt(nrow(overlaps)),
      cnt(length(gap_area)),
      cnt(nrow(holes_filled)),
      cnt(n_dropped),
      cnt(length(missing)),
      cnt(nrow(extra_parts)),
      cnt(nrow(lost_adj)),
      cnt(nrow(gained_adj)),
      cnt(nrow(far_off))
    ),
    detail = c(
      "",
      fmt(
        "total overlap area %.4g; sum(areas) - union = %.4g",
        sum(overlaps$area),
        total_area - union_area
      ),
      fmt(
        "%d hole(s) in original, %d in output; unexpected area %.4g",
        length(holes_o),
        length(holes_s),
        sum(gap_area)
      ),
      fmt("total area %.4g", sum(holes_filled$area)),
      fmt("total area %.4g", dropped_area),
      "",
      "",
      "",
      "",
      fmt("buffer = %.4g", buffer)
    ),
    stringsAsFactors = FALSE
  )
  summary$status <- ifelse(
    !is.na(summary$n_issues) & summary$n_issues == 0,
    "ok",
    "CHECK"
  )

  cat(sprintf(
    "Total area change (info only): %.3f%%\n\n",
    100 * (sum(a_s) - sum(a_o)) / sum(a_o)
  ))
  print(summary, row.names = FALSE)

  invisible(list(
    summary = summary,
    overlaps = overlaps,
    gaps = gap_area,
    holes_filled = holes_filled,
    dropped_faces = dropped,
    missing_features = missing,
    extra_parts = extra_parts,
    lost_adjacency = lost_adj,
    gained_adjacency = gained_adj,
    far_off = far_off,
    per_feature = per_feature
  ))
}

# Example
# smoothed <- smooth_shared(polys, window = 4, dens = 5e4, iter = 1,
#                           smooth_outer = FALSE)
# diag <- smooth_diagnostics(polys, smoothed, id = "my_id_column", buffer = 1e5)
