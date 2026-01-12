# metaSync: Synchronize Metadata Across Single-Cell Objects

An R package to merge and synchronize metadata (such as cell type annotations) from multiple single-cell objects to a main object using cell name matching.

## Problem

When working with single-cell RNA-seq data, you often have:
- A main integrated/combined object (`new_all_obj`)
- Multiple sub-objects with stratified analyses (`C0_obj`, `C1_obj`, etc.)
- Metadata annotations (like `celltype`) in sub-objects that need to be propagated back to the main object

Manually matching and merging these can be tedious and error-prone.

## Solution

`metaSync` provides a simple, reliable way to synchronize metadata across objects:

```r
# Install (when ready)
# devtools::install_github("xiaoqqjun/metaSync")

library(metaSync)

# Merge celltype annotations from multiple sub-objects to main object
merged_obj <- merge_metadata(
  main_obj = new_all_obj,
  sub_objects = list(C0_obj, C1_obj, C2_obj, C3_obj, C4_obj, C5_obj),
  metadata_col = "celltype",
  na.fill = "Unclassified",
  as.factor = TRUE,
  verbose = TRUE
)
```

## Features

- **Simple API**: One function does the heavy lifting
- **Cell name matching**: Automatically matches cells by colnames
- **Error reporting**: Detailed statistics on matching success
- **Flexible NA handling**: Fill unmatched cells with custom values
- **Factor conversion**: Optional automatic conversion to factor
- **Verbose output**: Track exactly what's being merged

## Parameters

- `main_obj`: The main Seurat object to receive metadata
- `sub_objects`: List of sub-objects containing metadata to merge
- `metadata_col`: Name of the metadata column to merge (default: "celltype")
- `by`: How to match cells ("cell_names" or "rownames")
- `na.fill`: Value for unmatched cells (default: "Unclassified")
- `as.factor`: Convert result to factor? (default: TRUE)
- `verbose`: Print statistics? (default: TRUE)

## Workflow Example

```r
# Your typical workflow:
# 1. Create main combined object
new_all_obj <- merge(C0_obj, y = c(C1_obj, C2_obj, ...))

# 2. Stratify and annotate in sub-objects
# (do analysis, annotate celltypes in each sub_obj)

# 3. Sync metadata back to main object with metaSync
new_all_obj <- merge_metadata(
  main_obj = new_all_obj,
  sub_objects = list(C0_obj, C1_obj, C2_obj, C3_obj, C4_obj, C5_obj),
  metadata_col = "celltype"
)

# 4. Proceed with downstream analysis
```

## Development

This package is actively developed. Feedback and suggestions are welcome!

### TODO
- [ ] Support for multiple metadata columns at once
- [ ] Conflict detection (same cell with different annotations)
- [ ] Batch metadata merging
- [ ] Better handling of sparse matching scenarios

## Citation

If you use metaSync in your research, please cite:

```
xiaoqqjun. (2025). metaSync: Synchronize Metadata Across Single-Cell Objects
```

## License

MIT License - See LICENSE file for details
