library(msa)
mySequences <- readAAStringSet("cdr3.fasta")
aln <- msa(mySequences)
aln2 <- msaConvert(aln, type="seqinr::alignment")
library(seqinr)
d <- dist.alignment(aln2, "identity")
clusters <- hclust(d)
groups <- cutree(clusters)
summary(clusters)

# Amir's idea
# cls = cumsum(table(factor(ct, levels = unique(ct[clusters$order])))) / length(ct)
# abline(h = cls, v = cls)
# ct  = cutree(clusters, 40)
# cls = cumsum(table(factor(ct, levels = unique(ct[clusters$order])))) / length(ct)
# image(X[ clusters$order, clusters$order])
# abline(h = cls, v = cls)
# mtext(names(cls), side = 2, las = 2, at = rowMeans(cbind(c(0,cls[-length(cls)]),cls)))
# names(which(ct == 25))