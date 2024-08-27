#' Obtain NOAA World Ocean Atlas dataset
#'
#' Retrieves data from the NOAA World Ocean Atlas.
#'
#' Functions to retrieve data from the
#' [NOAA World Ocean Atlas](https://www.ncei.noaa.gov/products/world-ocean-atlas)
#' . Data is a 3D array (longitude, latitude, and depth) and is loaded as a
#' [`stars`][stars::st_as_stars()] object. Check [`NOAA_data`] for available
#' variables, respective units and their citations. The function can
#' automatically cache the extracted files (default: `cache = FALSE`). The cached file will
#' then reside in the package's `extdata` directory.
#'
#' @param var The chemical or physical variable of interest (possible choices:
#'  `"temperature"`, `"phosphate"`, `"nitrate"`, `"silicate"`, `"oxygen"`,
#'  `"salinity"`, `"density"`).
#' @param spat_res Spatial resolution, either 0.25, 1, or 5-degree grid-cells (numeric).
#' @param av_period Temporal resolution, either `"annual"`, specific seasons
#'  (e.g. `"winter"`), or month (e.g. `"August"`).
#' @param cache Caching the extracted NOAA file in the package's `extdata`
#'  directory (default = `FALSE`). Size of individual files is around 12 Mb. Use
#'  [list_NOAA()] to list cached data resources.
#'
#' @return [`stars`][stars::st_as_stars()] object or path.
#' @export
#'
#' @examples
#'
#' # path to NOAA server or local data source
#' url_parser("oxygen", 1, "annual")
#'
#' if (interactive()) {
#'
#' # retrieve NOAA data
#' get_NOAA("oxygen", 1, "annual")
#'
#' }
get_NOAA <- function(var, spat_res, av_period, cache = FALSE) {

  # Abbreviate variable
  if (var == "silicate") {
    v <- strsplit(var, "")[[1]][2]
  } else if (var == "density") {
    v <- "I"
  } else {
    v <- strsplit(var, "")[[1]][1]
  }

  # Stat to extract
  if (spat_res == 5) {
    stat <- paste(v, "mn", sep = "_")
  } else {
    stat <- paste(v, "an", sep = "_")
  }

  # Where is package
  pkg_path <- fs::path_package("oceanexplorer")

  # Path
  NOAA_path <- url_parser(var, spat_res, av_period, cache = cache)

  if (!"external" %in% names(NOAA_path)) {
    # Get data and make stars object
    NOAA <- readRDS(fs::path(pkg_path, NOAA_path$local)) |> stars::st_as_stars()
  } else {
    # Get netcdf
    NOAA <- read_NOAA(NOAA_path$external, stat)

    if (isTRUE(cache)) {
      # Write stars object if extracted from NOAA server
      # Create dir
      fs::dir_create(pkg_path, fs::path_dir(NOAA_path$local))
      # Create file
      saveRDS(NOAA, fs::path(pkg_path, NOAA_path$local))
    }
  }

  # Return object
  NOAA
}

#' @rdname get_NOAA
#' @export
url_parser <- function(var, spat_res, av_period, cache = FALSE) {

  # Temporal resolution
  averaging_periods <- c("annual", month.name, "winter", "spring", "summer", "autumn")
  stopifnot(av_period %in% averaging_periods)

  # Base path to NCEI server
  base_path <- "https://data.nodc.noaa.gov/thredds/dodsC/ncei/woa"

  # Grouped variables
  chem <- c("phosphate", "nitrate", "silicate", "oxygen")
  misc <- c("temperature", "salinity", "density")

  # Recording range
  if (var %in% chem) {
    deca <- "all"
  } else if (var %in% misc) {
    deca <- "decav"
  }

  # Abbreviate variable
  if (var == "silicate") {
    v <- strsplit(var, "")[[1]][2]
  } else if (var == "density") {
    v <- "I"
  } else {
    v <- strsplit(var, "")[[1]][1]
  }

  # Averaging period
  tp <- grep(av_period, averaging_periods, ignore.case = TRUE) - 1
  tp <- sprintf(fmt = "%02.0f", tp)

  # Grid-cell size
  gr <- if (spat_res == 5) {
    "5d"
  } else if (spat_res == 1) {
    "01"
  } else {
    "04" # Quarter-degree resolution
  }

  # Complete file name
  file <- paste0(paste("woa18", deca, paste0(v, tp), gr, sep = "_"), ".nc")

  # Complete file path
  file_path <- paste(var, deca, if (spat_res == 5) "5deg" else if (spat_res == 1) "1.00" else "0.25", file, sep = "/")

  if (isTRUE(cache)) {
    # Where is package
    pkg_path <- fs::path_package("oceanexplorer")
    # Create extdata if not already existing
    if (!fs::dir_exists(fs::path(pkg_path, "extdata"))) {
      fs::dir_create(fs::path(pkg_path, "extdata"))
    }
  }

  # Check whether the file exists locally
  local_path <- fs::path("extdata", fs::path_ext_remove(file_path), ext = "rds")

  # If not, create external path to server
  noaa_path <- try(fs::path_package("oceanexplorer", local_path), silent = TRUE)
  if (!inherits(noaa_path, "fs_path")) {
    external_path <- paste(base_path, file_path, sep = "/")
    pt <- list(external = external_path)
    # Caching also adds local_path
    if (isTRUE(cache)) {
      pt <- append(pt, list(local = local_path))
    }
    pt
  } else {
    list(local = local_path)
  }
}

# Read the NOAA netcdf
read_NOAA <- function(conn, var) {

  # Make connection
  nc <- RNetCDF::open.nc(conn)

  # Variable
  lat <- RNetCDF::var.get.nc(nc, "lat")
  lon <- RNetCDF::var.get.nc(nc, "lon")
  depth <- RNetCDF::var.get.nc(nc, "depth")
  lat_bnds <- RNetCDF::var.get.nc(nc, "lat_bnds")
  lon_bnds <- RNetCDF::var.get.nc(nc, "lon_bnds")
  depth_bnds <- RNetCDF::var.get.nc(nc, "depth_bnds")
  attr <- RNetCDF::var.get.nc(nc, var)

  # Close connection
  RNetCDF::close.nc(nc)

  st <- stars::st_as_stars(attr) |>
    stars::st_set_dimensions(
      which = 1,
      offset = min(lon_bnds),
      delta = unique(diff(lon)),
      refsys = sf::st_crs(4326),
      names = "lon"
    ) |>
    stars::st_set_dimensions(
      which = 2,
      offset = min(lat_bnds),
      delta = unique(diff(lat)),
      refsys = sf::st_crs(4326),
      names = "lat"
    ) |>
    stars::st_set_dimensions(
      which = 3,
      values = depth_bnds[1, ],
      names = "depth"
    )

  # Variable name
  names(st) <- var
  st
}

#-------------------------------------------------------------------------------
# Auxiliary functions
#-------------------------------------------------------------------------------
list_NOAA <- function() {

  # Directory
  cache_pkg <- try(
    fs::path_package(package = "oceanexplorer", "extdata"),
    silent = TRUE
  )

  if (inherits(cache_pkg, "try-error")) {
    character(0)
  } else {
    # List files and delete them
    list.files(cache_pkg, full.names = TRUE)
  }
}

clean_cache <- function(...) {
  list_NOAA() |> fs::file_delete()
}


