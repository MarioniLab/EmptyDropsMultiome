
x = 1
y = 1
line = c(5,-1)

line_new = calc_lower_line(x,y, line)

test_that("multiplication works", {
  expect_equal(line_new, c(2,-1))
})

line = c(0,-1)

test_that("expect correct intercept", {
  expect_equal(calc_intercept_of_parallel_line(line, 1.41421356), c(2,-1))
})


test_that("expect correct intercept", {
  expect_equal(dist_point_to_line(line, c(1,1)), 1.41421356)
})

set.seed(0)
count_matrix <- DropletUtils:::simCounts()
colnames(count_matrix) <- paste0("C", as.character(seq(dim(count_matrix)[2]) ) )
count_matrix_rna <- as.matrix(count_matrix)
set.seed(1)
count_matrix1 <- DropletUtils:::simCounts()
colnames(count_matrix1) <- paste0("C", as.character(seq(dim(count_matrix1)[2]) ) )
count_matrix_atac <- as.matrix(count_matrix1)

Total_RNA =Matrix::colSums(count_matrix_rna)
Total_chromatin=Matrix::colSums(count_matrix_atac)

df <- data.frame("Total_RNA"=Total_RNA, "Total_chromatin"=Total_chromatin)
equation_k_means <- c(1,-1)
equation_parallel <- calc_ambiguous_above(df, equation_k_means)

test_that("slopes match", {
  expect_equal(equation_parallel[2], equation_k_means[2])
})

plot(log10(df$Total_chromatin + 0.1), log10(df$Total_RNA + 0.1),
     #     xlim=c(0,1510) , ylim=c(0,8),
     pch=".",
     # cex=2,
     #col=factor(df$FDR_RNA <= 0.001 & !is.na(eD.out_multi$FDR_RNA) ),
     xlab="log10(atac_count)", ylab="log10(rna_count)",
     main=paste0()
)
abline(a=equation_k_means[1],b=equation_k_means[2], col="blue")
abline(a=equation_parallel[1],b=equation_parallel[2], col="blue")
abline(a=3.13+3.15,b=-1, col="blue")
print(dist_point_to_line(equation_parallel, c(3.13,3.15)))
print(dist_point_to_line(equation_k_means, c(3.13,3.15)))

# plotBarcodes <- function(ranks, totals, fitted=NULL, subset=NULL, ...) {
#   xlim <- range(ranks[ranks>0])
#   ylim <- range(totals[totals>0])
# 
#   # Dropping non-unique points to save space.
#   # Subsetting performed after range() to create comparable plots, if desired.
#   keep <- !duplicated(totals)
#   if (!is.null(subset)) {
#     alt.keep <- keep
#     keep[] <- FALSE
#     keep[subset] <- alt.keep[subset]
#   }
#   Rx <- ranks[keep]
#   Tx <- totals[keep]
# 
#   # Actually making the plot, also plotting the fitted values if requested.
#   plot(Rx, Tx, log="xy", xlab="Rank", ylab="Total count", xlim=xlim, ylim=ylim, ...)
#   if (!is.null(fitted)) {
#     Fx <- fitted[keep]
#     o <- order(Rx)
#     lines(Rx[o], Fx[o], col="red", lwd=2)
#   }
#   return(invisible(NULL))
# }
# 

# logtotals <- log10(Total_RNA+0.1) +  log10(Total_chromatin+0.1)
# o <- order(logtotals, decreasing=TRUE)
# ot <- logtotals[o]
# r <- seq_along(o)
# plotBarcodes(r, ot, pch=16, main="logtotals")
#
#
# prod <- Total_RNA  *  Total_chromatin^(-equation_k_means[2])
# o <- order(prod, decreasing=TRUE)
# ot <- prod[o]
# r <- seq_along(o)
# plotBarcodes(r, ot, pch=16, main="prod")


# line = c(5,  -1.2)
# point = c(2,1)
# 
# calc_lower_line(line, point[1], point[2])


