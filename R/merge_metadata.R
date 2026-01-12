#' Synchronize Metadata Across Single-Cell Objects
#'
#' Merge metadata from multiple sub-objects to a main object based on cell name matching.
#' This is particularly useful for propagating cell type annotations from stratified
#' analyses back to a main combined object.
#'
#' @param main_obj The main object (typically a Seurat object or similar) to receive metadata
#' @param sub_objects A list of sub-objects containing the metadata to merge
#' @param metadata_col Character string. The name of the metadata column to merge.
#'   If NULL, will merge all columns except the specified exclusion columns.
#' @param by Character string. How to match cells: "cell_names" (default, uses colnames),
#'   or "rownames" (uses rownames for matching)
#' @param na.fill Value to fill cells not found in any sub_object. 
#'   Default is "Unclassified"
#' @param as.factor Logical. Should the result be converted to factor? Default is TRUE
#' @param verbose Logical. Print matching statistics? Default is TRUE
#' @param ... Additional arguments passed to methods
#'
#' @return The main_obj with updated metadata column
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' merged_obj <- merge_metadata(
#'   main_obj = main_seurat_obj,
#'   sub_objects = list(sub_obj1, sub_obj2, sub_obj3),
#'   metadata_col = "celltype",
#'   na.fill = "Unclassified"
#' )
#' }
#'
#' @export
merge_metadata <- function(
    main_obj,
    sub_objects,
    metadata_col = "celltype",
    by = "cell_names",
    na.fill = "Unclassified",
    as.factor = TRUE,
    verbose = TRUE,
    ...) {
  
  # Input validation
  if (!is.list(sub_objects)) {
    sub_objects <- list(sub_objects)
  }
  
  if (length(sub_objects) == 0) {
    stop("sub_objects must contain at least one object")
  }
  
  # Initialize metadata column as character
  if (!(metadata_col %in% colnames(main_obj@meta.data))) {
    main_obj@meta.data[[metadata_col]] <- NA_character_
  } else {
    main_obj@meta.data[[metadata_col]] <- as.character(main_obj@meta.data[[metadata_col]])
  }
  
  # Get cell identifiers based on 'by' parameter
  if (by == "cell_names") {
    main_cells <- colnames(main_obj)
  } else if (by == "rownames") {
    main_cells <- rownames(main_obj@meta.data)
  } else {
    stop("by must be 'cell_names' or 'rownames'")
  }
  
  # Track matching statistics
  total_matched <- 0
  sub_matched <- integer(length(sub_objects))
  
  # Iterate through sub_objects and merge metadata
  for (i in seq_along(sub_objects)) {
    sub_obj <- sub_objects[[i]]
    
    # Get cell identifiers from sub_object
    if (by == "cell_names") {
      sub_cells <- colnames(sub_obj)
    } else {
      sub_cells <- rownames(sub_obj@meta.data)
    }
    
    # Check if metadata column exists in sub_object
    if (!(metadata_col %in% colnames(sub_obj@meta.data))) {
      warning(sprintf("Object %d does not have '%s' column. Skipping.", i, metadata_col))
      next
    }
    
    # Match cells
    idx <- match(sub_cells, main_cells)
    
    # Find valid matches
    valid_matches <- !is.na(idx)
    n_matched <- sum(valid_matches)
    sub_matched[i] <- n_matched
    total_matched <- total_matched + n_matched
    
    # Assign metadata
    main_obj@meta.data[[metadata_col]][idx[valid_matches]] <- 
      as.character(sub_obj@meta.data[[metadata_col]][valid_matches])
    
    if (verbose) {
      cat(sprintf("Object %d: matched %d cells\n", i, n_matched))
    }
  }
  
  # Fill unmatched cells
  na_count <- sum(is.na(main_obj@meta.data[[metadata_col]]))
  main_obj@meta.data[[metadata_col]][is.na(main_obj@meta.data[[metadata_col]])] <- na.fill
  
  # Convert to factor if requested
  if (as.factor) {
    main_obj@meta.data[[metadata_col]] <- factor(main_obj@meta.data[[metadata_col]])
  }
  
  # Print summary statistics
  if (verbose) {
    cat("\n=== Merge Summary ===\n")
    cat(sprintf("Total cells in main_obj: %d\n", ncol(main_obj)))
    cat(sprintf("Total matched from sub_objects: %d\n", total_matched))
    cat(sprintf("Filled with '%s': %d\n", na.fill, na_count))
    cat(sprintf("Final distribution:\n"))
    print(table(main_obj@meta.data[[metadata_col]]))
  }
  
  return(main_obj)
}
