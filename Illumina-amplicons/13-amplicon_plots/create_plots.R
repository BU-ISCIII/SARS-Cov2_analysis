library(ggplot2)
library(tidyr)
library(stringr)


## Variables
# Directory with htseq results
directory <- "./"
# File names present in directory
sampleFiles <- grep("_coverage.txt",list.files(directory,recursive=TRUE,include.dirs=FALSE),value=TRUE)
# Get only sampleNames
sampleNames <- sub("(.*)/(.*)_coverage.txt","\\2",sampleFiles)

####MERGE DATA###### 
Final_table = NULL

for (index in 1:length(sampleNames)) {
  sample_name <- sampleNames[index]
  sample_file <- sampleFiles[index]
  sample_table <- read.table(file = sample_file, header = F, sep = "\t")
  if (index==1) {
    sample_table <- sample_table[,c(1,2,3,4,5,7,8)]
    header <- c("CHROMOSOME", "START", "END", "AMPLICON", "BASES", "AMPLICON_SIZE", "X1.COVERAGE")
    colnames(sample_table) <- header
    sample_table$DEPTH <- round(((sample_table$BASES * 150)/sample_table$AMPLICON_SIZE), digits = 3)
    sample_table$COVERAGE <- round((sample_table$X1.COVERAGE)*100, digits = 0)
    sample_table <- sample_table[,c(1,2,3,4,6,9,8)]
    header <- c("CHROMOSOME", "START", "END", "AMPLICON", "AMPLICON SIZE", paste("COVERAGE",sample_name, sep = "_"), paste("DEPTH",sample_name, sep = "_"))
    colnames(sample_table) <- header
    Final_table <- sample_table
  } else {
    sample_table <- sample_table[,c(4,5,7,8)]
    header <- c("AMPLICON", "BASES", "AMPLICON_SIZE", "X1.COVERAGE")
    colnames(sample_table) <- header
    sample_table$DEPTH <- round(((sample_table$BASES * 150)/sample_table$AMPLICON_SIZE), digits = 3)
    sample_table$COVERAGE <- round((sample_table$X1.COVERAGE)*100, digits = 0)
    sample_table <- sample_table[,c(1,6,5)]
    header <- c("AMPLICON", paste("COVERAGE",sample_name, sep = "_"), paste("DEPTH",sample_name, sep = "_"))
    colnames(sample_table) <- header
    Final_table <- merge(x = Final_table, y = sample_table, by.x = "AMPLICON", by.y = "AMPLICON", all = T)
  }   
}

write.table(x = Final_table, file = "./per_amplican_all_samples_coverage_depth.txt", col.names = T, row.names = F, quote = F, sep = "\t")

####PER SAMPLE PLOT COVERAGE###### 
for (index in 1:length(sampleNames)) {
  sample_name <- sampleNames[index]
  sample_file <- sampleFiles[index]
  sample_table <- read.table(file = sample_file, header = F, sep = "\t")
  sample_table <- sample_table[,c(2,4,8)]
  header <- c("START", "AMPLICON", "COVERAGE")
  colnames(sample_table) <- header
  plot_title <- paste(sample_name, "per amplicon coverage", sep = " ")
  sample_table$AMPLICON <- factor(sample_table$AMPLICON, levels = sample_table$AMPLICON[order(sample_table$START)])
  bp<-ggplot(data=sample_table, aes(x=AMPLICON, y=COVERAGE)) +
    geom_bar(stat="identity") + ggtitle(plot_title) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
  output_file_name_2 <- paste(sample_name, "coverage.pdf", sep = "_")
  output_file_name <- paste(sample_name, output_file_name_2, sep = "/")
  pdf(file=output_file_name)
  plot(bp)
  dev.off()
}

####PER SAMPLE PLOT DEPTH###### 
for (index in 1:length(sampleNames)) {
  sample_name <- sampleNames[index]
  sample_file <- sampleFiles[index]
  sample_table <- read.table(file = sample_file, header = F, sep = "\t")
  sample_table$DEPTH <- (sample_table$V5 * 150)/sample_table$V7
  sample_table <- sample_table[,c(2,4,9)]
  header <- c("START", "AMPLICON", "DEPTH")
  colnames(sample_table) <- header
  plot_title <- paste(sample_name, "per amplicon depth", sep = " ")
  sample_table$AMPLICON <- factor(sample_table$AMPLICON, levels = sample_table$AMPLICON[order(sample_table$START)])
  bp<-ggplot(data=sample_table, aes(x=AMPLICON, y=DEPTH)) +
    geom_bar(stat="identity") + ggtitle(plot_title) + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
  output_file_name_2 <- paste(sample_name, "depth.pdf", sep = "_")
  output_file_name <- paste(sample_name, output_file_name_2, sep = "/")
  pdf(file=output_file_name)
  plot(bp)
  dev.off()
}

####ALL SAMPLES PLOT COVERAGE###### 
merged_table <- Final_table
merged_table$POS <- round((merged_table$START + merged_table$END)/2)

rest_columns <- seq(6,ncol(merged_table)-2,2)
amplicon_coverage <- merged_table[,c(ncol(merged_table),rest_columns)]

data_long <- gather(amplicon_coverage, sample, coverage, 2:ncol(amplicon_coverage), factor_key=TRUE)
samples_names <- data_long$sample
samples_names <- str_remove(string = samples_names, pattern  = "COVERAGE_")
data_long$sample <- samples_names

plt <- ggplot(data = data_long, mapping = aes(x = POS, y = coverage, color = sample)) + 
  geom_line()+
  labs(title="Coverage per sample per position", y="Coverage", x="Postion", color="Sample name") +  # title and caption
  theme(panel.grid.minor = element_blank()) # turn off minor grid
pdf(file="all_samples_coverage.pdf")
plot(plt)
dev.off()

####ALL SAMPLES PLOT DEPTH###### 
rest_columns <- seq(7,ncol(merged_table)-1,2)
amplicon_depth <- merged_table[,c(ncol(merged_table),rest_columns)]

data_long <- gather(amplicon_depth, sample, depth, 2:ncol(amplicon_depth), factor_key=TRUE)
samples_names <- data_long$sample
samples_names <- str_remove(string = samples_names, pattern  = "DEPTH_")
data_long$sample <- samples_names

plt <- ggplot(data = data_long, mapping = aes(x = POS, y = depth, color = sample)) + 
  geom_line()+
  labs(title="Depth per sample per position", y="Depth", x="Postion", color="Sample name") +  # title and caption
  theme(panel.grid.minor = element_blank()) # turn off minor grid
pdf(file="all_samples_depth.pdf")
plot(plt)
dev.off()

####PER AMPLICON PLOT COVERAGE######
rest_columns <- seq(6,ncol(merged_table)-2,2)
per_amplicon_coverage <- merged_table[,c(1,ncol(merged_table),rest_columns)]
per_amplicon_coverage <- per_amplicon_coverage[order(per_amplicon_coverage$POS),]
per_amplicon_coverage$AMPLICON <- factor(per_amplicon_coverage$AMPLICON, levels = per_amplicon_coverage$AMPLICON[order(per_amplicon_coverage$POS)])
per_amplicon_coverage$POS <- NULL
pdf(file="all_amplicons_coverage.pdf")
for (row_num in 1:nrow(per_amplicon_coverage)) {
  new_table <- as.data.frame(per_amplicon_coverage[c(row_num),])
  amplicon_name <- as.character(per_amplicon_coverage[row_num,1])
  new_table$AMPLICON <- NULL
  new_table <- rbind(colnames(new_table),new_table)
  new_table <- t(new_table)
  new_table <- as.data.frame(new_table)
  plot_title <- paste(amplicon_name, "per sample coverage", sep = " ")
  colnames(new_table) <- c("Sample_name", "COVERAGE")
  new_table$COVERAGE <- as.numeric(as.character(new_table$COVERAGE))
  amplcn <- ggplot(data=new_table, aes(x=Sample_name, y=COVERAGE, group=1)) +
    geom_line() + ggtitle(plot_title) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
    ylim(0,100)
  plot(amplcn)
}
dev.off()

####PER AMPLICON PLOT DEPTH######
rest_columns <- seq(7,ncol(merged_table)-1,2)
per_amplicon_depth <- merged_table[,c(1,ncol(merged_table),rest_columns)]
per_amplicon_depth <- per_amplicon_depth[order(per_amplicon_depth$POS),]
per_amplicon_depth$AMPLICON <- factor(per_amplicon_depth$AMPLICON, levels = per_amplicon_depth$AMPLICON[order(per_amplicon_depth$POS)])
per_amplicon_depth$POS <- NULL
pdf(file="all_amplicons_depth.pdf")
for (row_num in 1:nrow(per_amplicon_depth)) {
  new_table <- as.data.frame(per_amplicon_depth[c(row_num),])
  amplicon_name <- as.character(per_amplicon_depth[row_num,1])
  new_table$AMPLICON <- NULL
  new_table <- rbind(colnames(new_table),new_table)
  new_table <- t(new_table)
  new_table <- as.data.frame(new_table)
  plot_title <- paste(amplicon_name, "per sample depth", sep = " ")
  colnames(new_table) <- c("Sample_name", "DEPTH")
  new_table$DEPTH <- as.numeric(as.character(new_table$DEPTH))
  amplcn <- ggplot(data=new_table, aes(x=Sample_name, y=DEPTH, group=1)) +
    geom_line() + ggtitle(plot_title) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5)) +
    ylim(0,max(per_amplicon_depth[,c(2:ncol(per_amplicon_depth))]))
  plot(amplcn)
}
dev.off()

save.image()
  