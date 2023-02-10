# Subset cells in a seurat object using meta data
# Check the active Ident 

object_name@active.ident

# First need to change "Ident" to name of meta.data column name 

object_name <- SetIdent(object_name, value = "meta.data_column_name")
# Example: sce.seurat <- SetIdent(sce.seurat, value = "major_clust")

# Subset cells based on "idents" in the column name set as the new Ident

new_object <- subset(object_name, idents = "desired_column_name")
# Example: L2_3 <- subset(sce.seurat, idents = "L2/3_CUX2")

# View the amount of cells in new subsetted object for each of the "idents" 
# in the current meta.data column. All idents should have 0 cells except 
# the for the ident used to make the subsetted object. 

table(new_object@meta.data$meta.data_column_name)
# Example: table(L2_3@meta.data$major_clust)

# Save the subsetted object as a .rds 
saveRDS(new_object, paste0("new_object.Rds")) 