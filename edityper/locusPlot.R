#!/usr/bin/env Rscript

library(utils)
library(graphics)

# The locusplot function
locusplot <- function(counts, fastq.name, num.reads, max.reads = NULL) {
  if (is.null(x = max.reads)) {
    max.reads <- num.reads
  }
  locus.colors <- c(
    'Coverage' = 'grey',
    'Insertions' = 'red',
    'Deletions' = 'green',
    'Mismatches' = 'blue'
  )
  par(xpd = TRUE, mar = c(4, 5, 4, 5))
  # Make the plot
  barplot(
    height = counts,
    col = locus.colors[rownames(x = counts)],
    border = NA,
    ylim = c(0, max.reads),
    main = fastq.name,
    axes = FALSE
  )
  # Add xaxis and label
  plot.max <- max(barplot(height = counts, plot = FALSE))
  axis(
    side = 1,
    at = as.integer(x = seq.int(from = 0, to = plot.max, length.out = 10)),
    las = 1
  )
  mtext(text = 'Position', side = 1, line = 2)
  # Add yaxis and label
  y.pos <- seq.int(from = 0, to = max.reads, length.out = 10)
  axis(
    side = 2,
    at = as.integer(x = y.pos),
    las = 1
  )
  mtext(text = 'Number of Reads', side = 2, line = 3.5)
  # Add yaxis2 and label
  axis(
    side = 4,
    at = y.pos,
    # labels = paste0(round(x = y.pos / num.reads * 100, digits = 2), '%'),
    labels = round(x = y.pos / num.reads * 100, digits = 2),
    las = 1
  )
  mtext(text = 'Percent', side = 4, line = 3.5, srt = 180)
  # Add a legend
  legend(
    'topleft',
    fill = locus.colors[rownames(x = counts)],
    legend = names(x = locus.colors[rownames(x = counts)])
  )
}

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(x = args) != 2) {
  stop("Usage: locusPlot.R counts.file num.reads")
}
counts.file <- args[1]
num.reads <- as.integer(x = args[2])

if (!file.exists(counts.file)) {
  stop(paste("Cannot find counts file", counts.file))
}

# Get some information
fastq.name <- basename(path = counts.file)
fastq.name <- strsplit(x = fastq.name, split = '_locus')
fastq.name <- unlist(x = fastq.name)[1]

# Read in and order the counts table
counts <- read.table(file = counts.file, header = TRUE, as.is = TRUE)
counts <- t(x = as.matrix(x = counts))
counts <- counts[order(apply(X = counts, MARGIN = 1, FUN = max), decreasing = FALSE), ]

# Make our plots
pdf(file = file.path(dirname(path = counts.file), paste0(fastq.name, '_locus.pdf')))
# Overarching locus plot
locusplot(counts = counts, fastq.name = fastq.name, num.reads = num.reads)
# # Locus plot zoomed in to only events
# locusplot(counts = counts, fastq.name = fastq.name, num.reads = num.reads, max.reads = max(counts[1:3, ]))
# Locus plot per event type
for (i in rownames(x = counts)[1:3]) {
  locusplot(
    counts = counts[c(i, 'Coverage'), ],
    fastq.name = fastq.name,
    num.reads = num.reads,
    max.reads = max(counts[i, ])
  )
}
dev.off()
