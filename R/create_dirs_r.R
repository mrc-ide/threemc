#' @title Recursively Create Missing Directories
#' @description Function to recursively create directories if any of the
#' directories in a provided path are missing. Somewhat equivalent to
#' \code{mkdir -p} in bash.
#' @param dir_path Path to a file or directory which you want to generate.
#' @export
#'
create_dirs_r <- function(dir_path) {

  # split by "/"
  dirs <- stringr::str_split(dir_path, "/")[[1]]

  # check if dir_path is for dir (required) or specific file
  # does our string end in "/"? Then likely a dir
  cond <- substr(dir_path, nchar(dir_path), nchar(dir_path))
  cond <- grepl("/", cond)
  # does last word contain "."? Likely a file
  cond <- cond & !grepl(".", last(dirs))
  if (!cond) {
    dirs <- dirs[-length(dirs)] # remove file name
  }
  if (length(dirs) == 0) stop("No missing directories created")
  # loop through each directory level in target directory, create any missing
  for (i in seq_along(dirs)) {
    # "recursively" check directories
    if (i == 1) {
      spec_dir <- dirs[i]
    } else {
      spec_dir <- paste(dirs[1:i], collapse = "/")
    }
    # create if missing
    if (!dir.exists(spec_dir)) {
      dir.create(spec_dir)
    }
  }
}
