# metaSync - 快速测试脚本
# 这个脚本演示了如何使用metaSync包

# ===== 设置环境 =====
# 首先加载包（安装后）
# devtools::load_all("path/to/metaSync")
# 或
# library(metaSync)

# ===== 创建演示数据 =====
# 这里是创建虚拟Seurat对象用于演示的代码
# 实际使用时，你会用真实的单细胞对象

create_demo_object <- function(n_cells = 100, name = "demo") {
  set.seed(42)
  
  # 创建虚拟expression matrix
  expr_mat <- matrix(
    rpois(5000 * n_cells, lambda = 0.5),
    nrow = 5000,
    ncol = n_cells
  )
  rownames(expr_mat) <- paste0("Gene_", 1:5000)
  colnames(expr_mat) <- paste0(name, "_cell_", 1:n_cells)
  
  # 创建Seurat对象
  obj <- Seurat::CreateSeuratObject(counts = expr_mat)
  return(obj)
}

# ===== 演示用法 =====
demo_merge_metadata <- function() {
  # 创建主对象
  main_obj <- create_demo_object(n_cells = 300, name = "main")
  
  # 创建3个子对象，各100个细胞
  sub_obj1 <- create_demo_object(n_cells = 100, name = "main")
  sub_obj2 <- create_demo_object(n_cells = 100, name = "main")
  sub_obj3 <- create_demo_object(n_cells = 100, name = "main")
  
  # 添加celltype信息到子对象
  sub_obj1$celltype <- "Type_A"
  sub_obj2$celltype <- c(rep("Type_B", 80), rep("Type_C", 20))
  sub_obj3$celltype <- "Type_D"
  
  # 使用metaSync合并metadata
  merged_obj <- merge_metadata(
    main_obj = main_obj,
    sub_objects = list(sub_obj1, sub_obj2, sub_obj3),
    metadata_col = "celltype",
    na.fill = "Unknown",
    as.factor = TRUE,
    verbose = TRUE
  )
  
  return(merged_obj)
}

# ===== 真实使用示例 =====
# 这是用你的实际数据的方式：
# 
# library(metaSync)
# library(Seurat)
# 
# # 你的对象已经存在：
# # new_all_obj, C0_obj, C1_obj, ..., C5_obj
# 
# merged_obj <- merge_metadata(
#   main_obj = new_all_obj,
#   sub_objects = list(C0_obj, C1_obj, C2_obj, C3_obj, C4_obj, C5_obj),
#   metadata_col = "celltype",
#   na.fill = "Unclassified",
#   as.factor = TRUE,
#   verbose = TRUE
# )
# 
# # 检查结果
# table(merged_obj$celltype)
