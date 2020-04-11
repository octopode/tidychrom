#' Unlist Columns
#'
#' Take a tibble, unlist all list columns as able (e.g. S3 objects will remain wrapped in a list.)
#' @param tbl_in
#' @keywords unlist list columns tibble
#' @export
#' @examples
#' # effectively add an annotation column to peaks_matched:
#' peaks_matched <- unlist_columns(peaks_matched, singletons = T)
