load('junctions.rda')


# convert bedtools multicov outputs to CSV spreadsheets
bam.files <- read.delim('./BAM_files', header = FALSE, stringsAsFactors = FALSE)
ss5 <- read.delim('./introns/introns_known_5ss.bed.counts', header = FALSE, stringsAsFactors = FALSE)
colnames(ss5)[1:6] <- c('chrom', 'start', 'end', 'name', 'score', 'strand')
colnames(ss5)[7:ncol(ss5)] <- colnames(junc.all)[(ncol(junc.all) - nrow(bam.files) + 1):ncol(junc.all)]
write.csv(ss5, file = './introns/introns_known_5ss.bed.counts.csv', row.names = FALSE, quote = TRUE)
ss3 <- read.delim('./introns/introns_known_3ss.bed.counts', header = FALSE, stringsAsFactors = FALSE)
colnames(ss3) <- colnames(ss5)
write.csv(ss3, file = './introns/introns_known_3ss.bed.counts.csv', row.names = FALSE, quote = TRUE)


# generate diagnostic plots of 3'ss and  5'ss coverage
pseudocount <- 0.1
for (i in 1:nrow(bam.files)){
  pdf(paste('./images/', colnames(ss3)[i + 6], '_intron_coverage.pdf', sep = ''))
  plot(ss3[, i + 6] + pseudocount, ss5[, i + 6] + pseudocount,
       pch = 20, log = 'xy',
       xlab = "3' end counts", 
       ylab = "5' end counts", 
       main = paste(colnames(ss3)[i + 6], '(intron counts)', sep = '\n'),
       sub = paste('Pearson R=', round(cor(ss3[, i + 6] + pseudocount, ss5[, i + 6] + pseudocount), 3), sep = ''))
  abline(h = 0.5, col = 'red')
  abline(v = 0.5, col = 'red')
  plot((ss3[, i + 6] + pseudocount) / (ss5[, i + 6] + pseudocount),
       pch = 20, log = 'y',
       xlab = "genes", 
       ylab = "3' end counts / 5' end counts", 
       main = colnames(ss3)[i + 6])
  abline(h = 1, col = 'blue')
  dev.off()
}


# calculate splicing efficiency for each intron (splice site coverage / transread count)
samples <- colnames(ss5)[7:ncol(ss5)]
# all known introns
efficiency.5ss <- cbind(ss5[1:6], junc.known[, samples] / ss5[, samples])
colnames(efficiency.5ss) <- colnames(ss5)
write.csv(efficiency.5ss, file = './efficiency/splicing_efficiency_5ss.csv', row.names = FALSE, quote = TRUE)
efficiency.3ss <- cbind(ss3[1:6], junc.known[, samples] / ss3[, samples])
colnames(efficiency.5ss) <- colnames(ss5)
write.csv(efficiency.3ss, file = './efficiency/splicing_efficiency_3ss.csv', row.names = FALSE, quote = TRUE)
# filter out unreliable calculations for junctions with low read coverage
threshold.counts <- 5 # minimum number of transreads & reads at 3' / 5' intron end required to treat any calculations as valid
confidence.5ss <- matrix(TRUE, nrow(junc.known), length(samples))
colnames(confidence.5ss) <- samples
confidence.5ss[junc.known[, samples] < threshold.counts] <- FALSE
confidence.5ss[ss5[, samples] < threshold.counts] <- FALSE
confidence.3ss <- matrix(TRUE, nrow(junc.known), length(samples))
colnames(confidence.3ss) <- samples
confidence.3ss[junc.known[, samples] < threshold.counts] <- FALSE
confidence.3ss[ss3[, samples] < threshold.counts] <- FALSE
efficiency.5ss.conf <- efficiency.5ss
efficiency.3ss.conf <- efficiency.3ss
for (i in samples){
  efficiency.5ss.conf[confidence.5ss[, i] == FALSE, i] <- NA
  efficiency.3ss.conf[confidence.3ss[, i] == FALSE, i] <- NA
}
write.csv(efficiency.5ss.conf, file = './efficiency/splicing_efficiency_5ss_conf.csv', row.names = FALSE, quote = TRUE)
write.csv(efficiency.3ss.conf, file = './efficiency/splicing_efficiency_3ss_conf.csv', row.names = FALSE, quote = TRUE)


# compare mutant vs WT (relative splicing efficiency for each intron)
comparisons <- t(matrix(c('cwf10_1', 'WT_1'), 
                        nrow = 2, ncol = 1)) # prepare a table of sample pairs to be compared: c(mutant, WT)
# all known introns
efficiency.5ss.rel <- cbind(efficiency.5ss[, c('chrom', 'start', 'end', 'name', 'score', 'strand')], as.list(rep(NA, nrow(comparisons))))
efficiency.3ss.rel <- cbind(efficiency.3ss[, c('chrom', 'start', 'end', 'name', 'score', 'strand')], as.list(rep(NA, nrow(comparisons))))
for (i in 1:nrow(comparisons)){
  efficiency.5ss.rel[, i + 6] <- efficiency.5ss[, comparisons[i, 1]] / efficiency.5ss[, comparisons[i, 2]]
  colnames(efficiency.5ss.rel)[i + 6] <- paste(comparisons[i, 1], comparisons[i, 2], sep = '_')
  efficiency.3ss.rel[, i + 6] <- efficiency.3ss[, comparisons[i, 1]] / efficiency.3ss[, comparisons[i, 2]]
  colnames(efficiency.3ss.rel)[i + 6] <- paste(comparisons[i, 1], comparisons[i, 2], sep = '_')
}
write.csv(efficiency.5ss.rel, file = './efficiency/relative_splicing_efficiency_5ss.csv', row.names = FALSE, quote = TRUE)
write.csv(efficiency.3ss.rel, file = './efficiency/relative_splicing_efficiency_3ss.csv', row.names = FALSE, quote = TRUE)
# introns with coverage above threshold
efficiency.5ss.rel.conf <- cbind(efficiency.5ss.conf[, c('chrom', 'start', 'end', 'name', 'score', 'strand')], as.list(rep(NA, nrow(comparisons))))
efficiency.3ss.rel.conf <- cbind(efficiency.3ss.conf[, c('chrom', 'start', 'end', 'name', 'score', 'strand')], as.list(rep(NA, nrow(comparisons))))
for (i in 1:nrow(comparisons)){
  efficiency.5ss.rel.conf[, i + 6] <- efficiency.5ss.conf[, comparisons[i, 1]] / efficiency.5ss.conf[, comparisons[i, 2]]
  colnames(efficiency.5ss.rel.conf)[i + 6] <- paste(comparisons[i, 1], comparisons[i, 2], sep = '_')
  efficiency.3ss.rel.conf[, i + 6] <- efficiency.3ss.conf[, comparisons[i, 1]] / efficiency.3ss.conf[, comparisons[i, 2]]
  colnames(efficiency.3ss.rel.conf)[i + 6] <- paste(comparisons[i, 1], comparisons[i, 2], sep = '_')
}
write.csv(efficiency.5ss.rel.conf, file = './efficiency/relative_splicing_efficiency_5ss_conf.csv', row.names = FALSE, quote = TRUE)
write.csv(efficiency.3ss.rel.conf, file = './efficiency/relative_splicing_efficiency_3ss_conf.csv', row.names = FALSE, quote = TRUE)


# prepare figures for paper
for (i in 1:nrow(comparisons)){
  pdf(paste('./images/', comparisons[i, 1], '_vs_', comparisons[i, 2], '_splicing_efficiency.pdf', sep = ''))
  par(mar = c(5, 5, 5, 5), cex = 1.2, lwd = 2)
  limits.data <- cbind(efficiency.5ss[, comparisons[i, ]], efficiency.3ss[, comparisons[i, ]])
  limits.data[limits.data == 0 | limits.data == Inf] <- NA
  limits <- range(limits.data, na.rm = TRUE)
  # 5' splice site
  plot(log2(efficiency.5ss[, comparisons[i, 1]]), log2(efficiency.5ss[, comparisons[i, 2]]), 
       main = paste('Splicing efficiency for each intron\n(transreads / 5\' splice site reads)', sep = '\n'),
       ylab = paste(comparisons[i, 2], '\nlog2(splicing efficiency)', sep = ''),
       xlab = paste(comparisons[i, 1], '\nlog2(splicing efficiency)', sep = ''),
       pch = 21, col = 'darkgrey',
       xlim = log2(limits), ylim = log2(limits))
  points(log2(efficiency.5ss.conf[, comparisons[i, 1]]), log2(efficiency.5ss.conf[, comparisons[i, 2]]), 
         pch = 19, col = 'black')
  abline(a = 0, b = 1)
  legend('bottomright', legend = c(paste('>= ', threshold.counts, ' reads', sep = ''),
                                   paste('< ', threshold.counts, ' reads', sep = '')),
         title = 'coverage',
         bty = 'n', pch = c(19, 21), col = c('black', 'dark gray'))
  # 3' splice site
  plot(log2(efficiency.3ss[, comparisons[i, 1]]), log2(efficiency.3ss[, comparisons[i, 2]]),
       main = paste('Splicing efficiency for each intron\n(transreads / 3\' splice site reads)', sep = '\n'),
       ylab = paste(comparisons[i, 2], '\nlog2(splicing efficiency)', sep = ''),
       xlab = paste(comparisons[i, 1], '\nlog2(splicing efficiency)', sep = ''),
       pch = 21, col = 'darkgrey',
       xlim = log2(limits), ylim = log2(limits))
  points(log2(efficiency.3ss.conf[, comparisons[i, 1]]), log2(efficiency.3ss.conf[, comparisons[i, 2]]), 
         pch = 19, col = 'black')
  abline(a = 0, b = 1)
  legend('bottomright', legend = c(paste('>= ', threshold.counts, ' reads', sep = ''),
                                   paste('< ', threshold.counts, ' reads', sep = '')),
         title = 'coverage',
         bty = 'n', pch = c(19, 21), col = c('black', 'dark gray'))
  
  dev.off()
}