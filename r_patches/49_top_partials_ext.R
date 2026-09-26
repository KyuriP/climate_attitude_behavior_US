ut <- which(upper.tri(W_ext), arr.ind = TRUE)
d <- data.frame(a = rownames(W_ext)[ut[,1]], b = colnames(W_ext)[ut[,2]], w = W_ext[ut])
d <- d[order(-abs(d$w)), ]; d$rank <- seq_len(nrow(d)); print(head(transform(d, w = round(w, 4)), 14), row.names = FALSE)
cat("n retained ext:", sum(W_ext[upper.tri(W_ext)] != 0), "of", choose(nrow(W_ext),2), "\n49 DONE\n")
