# metaSync - 安装和使用指南

## 包结构

```
metaSync/
├── DESCRIPTION          # 包的元信息
├── NAMESPACE            # 导出的函数
├── LICENSE              # MIT许可证
├── README.md            # 项目说明
├── R/
│   └── merge_metadata.R # 核心函数
├── man/                 # 帮助文档（会自动生成）
└── data-raw/
    └── demo_script.R    # 演示脚本
```

## 安装方式

### 方式1：本地安装（开发中）

```r
# 从GitHub克隆或者从本地文件夹
devtools::load_all("path/to/metaSync")

# 或构建和安装
devtools::install("path/to/metaSync")
```

### 方式2：从GitHub安装（发布后）

```r
# 需要先将包发布到GitHub
devtools::install_github("leone/metaSync")
```

## 快速开始

### 基础用法

```r
library(metaSync)
library(Seurat)

# 你的对象（Seurat对象）
# new_all_obj: 主对象
# C0_obj, C1_obj, ..., C5_obj: 子对象，含有celltype标注

# 合并metadata
merged_obj <- merge_metadata(
  main_obj = new_all_obj,
  sub_objects = list(C0_obj, C1_obj, C2_obj, C3_obj, C4_obj, C5_obj),
  metadata_col = "celltype"
)

# 检查结果
table(merged_obj$celltype)
```

### 高级选项

```r
# 自定义填充值和其他参数
merged_obj <- merge_metadata(
  main_obj = new_all_obj,
  sub_objects = list(C0_obj, C1_obj, C2_obj, C3_obj, C4_obj, C5_obj),
  metadata_col = "celltype",
  by = "cell_names",              # 或 "rownames"
  na.fill = "Unclassified",       # 未匹配的细胞用什么标记
  as.factor = TRUE,               # 转换为因子？
  verbose = TRUE                  # 显示详细信息？
)
```

## 函数参数详解

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `main_obj` | 主对象（将接收metadata）| 必需 |
| `sub_objects` | 子对象列表 | 必需 |
| `metadata_col` | 要合并的metadata列名 | "celltype" |
| `by` | 匹配方式："cell_names" 或 "rownames" | "cell_names" |
| `na.fill` | 未匹配细胞的填充值 | "Unclassified" |
| `as.factor` | 结果转换为因子? | TRUE |
| `verbose` | 显示详细信息? | TRUE |

## 返回值

返回修改后的主对象，metadata列已更新。

## 工作流示例

```r
library(metaSync)
library(Seurat)

# 1. 加载或创建数据
load("your_data.RData")  # 包含new_all_obj, C0_obj, ..., C5_obj

# 2. 应用metaSync
new_all_obj <- merge_metadata(
  main_obj = new_all_obj,
  sub_objects = list(C0_obj, C1_obj, C2_obj, C3_obj, C4_obj, C5_obj),
  metadata_col = "celltype",
  na.fill = "Unclassified"
)

# 3. 验证结果
table(new_all_obj$celltype)

# 4. 继续下游分析
Ident(new_all_obj) <- "celltype"
DimPlot(new_all_obj, label = TRUE, repel = TRUE)

# 5. 保存
saveRDS(new_all_obj, "new_all_obj_annotated.rds")
```

## 故障排除

### 问题1：找不到某些细胞

**症状**：提示信息显示匹配数小于预期，有很多"Unclassified"

**原因**：细胞名称不匹配（大小写、空格、特殊字符差异）

**解决**：
```r
# 检查细胞名称
head(colnames(new_all_obj))
head(colnames(C0_obj))

# 如果名称不同，需要先标准化
```

### 问题2：错误："因子层次有错，产生了NA"

**原因**：metadata列已经是因子，不能赋予新值

**解决**：
```r
# 先转换为字符
new_all_obj$celltype <- as.character(new_all_obj$celltype)

# 然后使用metaSync
new_all_obj <- merge_metadata(new_all_obj, ...)
```

## 提计划和建议

如果你有以下想法，欢迎反馈：
- [ ] 同时合并多个metadata列
- [ ] 冲突检测（同一细胞来自不同对象有不同标注）
- [ ] 批量处理多个主对象
- [ ] 更灵活的匹配方式（如基于特征ID）
- [ ] 自动celltype标准化

## 许可证

MIT License
