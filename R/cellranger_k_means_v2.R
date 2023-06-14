#' Finds a nice initialization for the K-means clustering
#'
#' @param x Vector of real values.
#'
#' @return Value roughly equal to 1/10 of the 99% percentile of values.
#'
simpler_ordmag <- function(x){
  # This function finds a nice initialization for the K-means clustering:
  # finds 99% percentile cut-off, and retains cell with count>99%cutoff/10;
  # then removes the top 1% of these;
  # if that removes all cells, then return the initial 99%cutoff/10, else, then compute the new 99% percentile cutoff
  # and return new99%cutoff/10
  # x: vector with counts per cell
  # output: count number that separates tentative empty droplets and cell-containing droplets

  p99_1 = stats::quantile(x, probs = c(0.99), na.rm=TRUE )
  ordmag = x > p99_1 / 10.0
  if (sum(ordmag) == 0){
    return(0)
  }
  top_1_ordmag = stats::quantile(x[ordmag], 0.99, na.rm=TRUE)
  ordmag2 = x < top_1_ordmag
  # use the previous threshold
  if (sum(ordmag2) == 0){
    return(p99_1 / 10.0)
  }
  p99_2 = stats::quantile(x[ordmag2], 0.99)
  #print(p99_2)

  return(p99_2 / 10.0)

}

#' Log-transform the data.
#'
#' @param data Dataframe with column atac_count and column rna_count
#'
#' @return Dataframe with with log-transformed column values
#'
log_transform <- function(data){
  # data: dataframe with column "atac_count" and "rna_count"
  # output: dataframe with log10 transformed columns
  data["atac_count"] <- log10(data["atac_count"]+0.1)
  data["rna_count"] <- log10(data["rna_count"]+0.1)
  return(data)
}




#' Assigns the value is_cell of the unique elements to the value is_cell of their corresponding duplicates.
#'
#' @param data Dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell' (bool), 'dup' (bool).
#'             data is sorted wrt 'atac_count' (int), 'rna_count' (int), 'dup' (bool).
#' @return Dataframe with altered is_cell column.
#'
transfer_labels_to_dedupl_v3 <- function(data){
  # data: dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell' (bool), 'dup' (bool)
  # data is sorted wrt 'atac_count' (int), 'rna_count' (int), 'dup' (bool)
  # take in the data dataframe and pass to the duplicates the labels you've calculated for the corresponding unduplicated barcode
  # return data

  print(colnames(data))

  uniques = which(data$dup==F)
  uniques = append(uniques, length(data$dup)+1)
  numb_of_duplicates = diff(uniques)-1

  i=1
  for (cell_no in uniques[seq(1, length(uniques)-1)]){
    if (numb_of_duplicates[i]!=0){
      data$is_cell[seq(cell_no+1, cell_no+numb_of_duplicates[i])]  = data$is_cell[cell_no]
    }
    i = i+1
  }
  
  return(data)
}



#' Call cells using k-means and store this information in is_cell.
#'
#' @param data Dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell'.
#' @param verbose if TRUE various intermediate steps and plots are printed.
#'
#' @return Dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell'.
#' @export
#'
call_cells <- function(data, verbose=TRUE){
  # data: dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell',
  #       e.g.   df <- data.frame( "rna_count"= seurat_024@meta.data[["nCount_RNA"]], "atac_count"= seurat_024@meta.data[["nCount_ATAC"]],
  #              "excluded"= rep(F,each=length(colnames(seurat_024))),  "is_cell"= rep(0,each=length(colnames(seurat_024)))    )
  # output data: dataframe with columns 'atac_count' (int), 'rna_count' (int), 'excluded' (bool), 'is_cell' (bool), "dup" (bool)
  # where the is_cell column has been calculated using K-means

  print("Call cells using cellranger-arc")
  data <- data[order(data$atac_count, data$rna_count, data$excluded),]
  exclude_rowname <- c("atac_count", "rna_count", "excluded", "is_cell")
  data["dup"] <- duplicated(data[, exclude_rowname])
  print("finished ordering")
  keep <- !data["dup"] & !data$excluded & data$atac_count>0 & data$rna_count>0   # !data$is_cell & !data["dup"]       #mask(data)
  print("calculated keep")


  threshold = vector(mode="numeric", length=2)
  names(threshold) <- c("atac", "rna")

  # 1: do ordmag calling and identify centroids
  for (assay in c("atac", "rna")){
    threshold[assay]= simpler_ordmag(data[keep, paste0(assay,"_count")])
    ordmag_filter <- data$atac_count > threshold["atac"] & data$rna_count > threshold["rna"]
  }

  print(paste0("thresholds are ", threshold))

  # if ordmag calls everything a cell, then do that
  if (length(unique(ordmag_filter))==1){
    data[keep, "is_cell"]=1

    print("all_called by ordmag")
    print(unique(ordmag_filter))

    data <- transfer_labels_to_dedupl_v3(data)
    print("finished label transfer")

    equation = c(0, 0)
  }else{
    set.seed(123)
    # Do scaling i.e. log transform
    df_for_km <- log_transform(data[,c("rna_count", "atac_count")])

    mean_atac_cell <- mean(df_for_km[keep & ordmag_filter, c("atac_count")])
    mean_rna_cell <- mean(df_for_km[keep & ordmag_filter, c("rna_count")])
    mean_atac_empty <- mean(df_for_km[keep & !ordmag_filter, c("atac_count")])
    mean_rna_empty <- mean(df_for_km[keep & !ordmag_filter, c("rna_count")])

    print("the means of the two clusters are: ")
    print(  data.frame( "rna_count"=c(mean_rna_empty ,mean_rna_cell), "atac_count"=c(mean_atac_empty, mean_atac_cell ) )   )

    km <- stats::kmeans( df_for_km[keep, c("rna_count", "atac_count")], centers=data.frame( "rna_count"=c(mean_rna_empty ,mean_rna_cell), "atac_count"=c(mean_atac_empty, mean_atac_cell ) )    )
    #print(paste0("number of cells in km object is ", sum(km[["cluster"]])-length(km[["cluster"]]) ) )

    # take center information and then exchange rna and atac (since atac needs to be on the x axis)
    center1 = c( km$centers[1,2], km$centers[1,1] )
    center2 = c( km$centers[2,2], km$centers[2,1] )

    # calculate the perpendicular bisector line to the segment defined by the two centers (i.e. the k-means line)
    equation = perp_bisector(center2, center1)
    print(equation)


    # if k-means calls everything a cell, do that
    if (length(unique(km[["cluster"]]))==1 ){
      data[keep, "is_cell"]=1
    }else{
      mean_1 <- mean(km[["centers"]][1])
      mean_2 <- mean(km[["centers"]][2])
      km[["cluster"]] <- km[["cluster"]]-1
      if (mean_1 > mean_2){
        km[["cluster"]] <- 1-km[["cluster"]]
      }
      data[keep, "is_cell"]= km[["cluster"]]

    }
  }

    
#   # plot cell calling before label transfer
#   cell_calling_plot(log10(data$atac_count + 0.1), log10(data$rna_count + 0.1), factor(data$is_cell ), paste0("cR-arc calling before label transfer: ", sum(data$is_cell), " cells;", equation))
#   graphics::abline(a=equation[1],b=equation[2], col="blue")

  # transfer labels
  data <- transfer_labels_to_dedupl_v3(data)
  print("finished label transfer")
    
#   # plot cell calling after label transfer
#   cell_calling_plot(log10(data$atac_count + 0.1), log10(data$rna_count + 0.1), factor(data$is_cell ), paste0("cR-arc calling after label transfer: ", sum(data$is_cell), " cells;", equation))
#   graphics::abline(a=equation[1],b=equation[2], col="blue")

  print("The number of cells called by k-means is ")
  print(paste0(sum(data$is_cell), " cells"))


  return(list(data, equation) )

}
