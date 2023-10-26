transread.folder <- './transreads/'
transread.file.extension <- '.trans.annotated'
transread.files <- list.files(transread.folder, pattern = transread.file.extension)


# compile a complete collection of splice junctions identified by regtools
first.run <- TRUE
for (i in transread.files) {
  junc.i <- read.delim(file = paste(transread.folder, i, sep = ''), header = TRUE, stringsAsFactors = FALSE)
  junc.i <- junc.i[, colnames(junc.i)[!colnames(junc.i) %in% c('name', 'score')]] # drop junction names and score
  if (first.run == TRUE) {
    junc.all <- junc.i
    first.run <- FALSE
  }
  else {
    junc.all.ID <- paste(junc.all$chrom, junc.all$start, junc.all$end, junc.all$strand)
    junc.i.ID <- paste(junc.i$chrom, junc.i$start, junc.i$end, junc.i$strand)
    junc.all <- rbind(junc.all, junc.i[!junc.i.ID %in% junc.all.ID, ])
  }
}


# extract numbers of transreads for all splice junctions from each sample
junc.all.ID <- paste(junc.all$chrom, junc.all$start, junc.all$end, junc.all$strand)
for (i in transread.files) {
  junc.all <- cbind(junc.all, 0)
  colnames(junc.all)[length(colnames(junc.all))] <- strsplit(i, split = '.', fixed = TRUE)[[1]][1]
  junc.i <- read.delim(file = paste(transread.folder, i, sep = ''), header = TRUE, stringsAsFactors = FALSE)
  junc.i.ID <- paste(junc.i$chrom, junc.i$start, junc.i$end, junc.i$strand)
  for (j in 1:length(junc.all.ID)){
    if (junc.all.ID[j] %in% junc.i.ID) {
      junc.all[j, ncol(junc.all)] <- junc.i[junc.i.ID == junc.all.ID[j], 'score']
    }
  }
}
junc.all <- junc.all[order(junc.all$chrom, junc.all$start, junc.all$end, junc.all$strand), ]
write.csv(junc.all, 
          file = paste(transread.folder, 'splice_junctions_coverage.all.csv', sep = ''),
          row.names = FALSE, quote = TRUE)
junc.known <- junc.all[junc.all$known_junction == 1, ]
write.csv(junc.known, 
          file = paste(transread.folder, 'splice_junctions_coverage.known.csv', sep = ''),
          row.names = FALSE, quote = TRUE)


# extract intron start/end coordinates (5' and 3' splice sites) and convert to BED format
junc.known.5ss <- as.data.frame(cbind(junc.known[, 'chrom'],
                                          junc.known[, 'start'], # BED conversion
                                          junc.known[, 'start'] + 1, # BED conversion
                                          junc.known[, 'transcripts'],
                                          '0',
                                          junc.known[, 'strand']),
                                    stringsAsFactors = FALSE)
colnames(junc.known.5ss) <- c('chrom', 'start', 'end', 'name', 'score', 'strand')
junc.known.3ss <- as.data.frame(cbind(junc.known[, 'chrom'],
                                          junc.known[, 'end'] - 2, # BED conversion
                                          junc.known[, 'end'] - 1, # BED conversion
                                          junc.known[, 'transcripts'],
                                          '0',
                                          junc.known[, 'strand']),
                                    stringsAsFactors = FALSE)
colnames(junc.known.3ss) <- c('chrom', 'start', 'end', 'name', 'score', 'strand')


# swap "start" and "end" for genes on "-" strand to account for gene orientation
temp.5ss <- junc.known.5ss
junc.known.5ss[junc.known.5ss$strand == '-', c('start', 'end')] <- junc.known.3ss[junc.known.3ss$strand == '-', c('start', 'end')]
junc.known.3ss[junc.known.3ss$strand == '-', c('start', 'end')] <- temp.5ss[temp.5ss$strand == '-', c('start', 'end')]
rm(temp.5ss)
write.table(junc.known.5ss, file = './genome/introns_known_5ss.bed', 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
write.table(junc.known.3ss, file = './genome/introns_known_3ss.bed',
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')
save(junc.all, junc.known, junc.known.5ss, junc.known.3ss, file = 'junctions.rda')