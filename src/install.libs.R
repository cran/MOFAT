files <- Sys.glob(paste0("*", SHLIB_EXT))
destination <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
if (length(files)) {
    copied <- file.copy(files, destination, overwrite = TRUE)
    if (!all(copied)) stop("Failed to install one or more MOFAT shared libraries.")
}
if (file.exists("symbols.rds")) {
    file.copy("symbols.rds", destination, overwrite = TRUE)
}
