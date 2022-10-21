file_names <- list.files(pattern = "\\coverage$")
files <- lapply(file_names, function(x) read.table(x, header = FALSE)[,c(1,3)])
files <- suppressWarnings(Reduce(function(x, y) merge(x, y, by = "V1"), files))

row.names(files) <- files$V1; files <- files[,-1]
colnames(files) <- c(file_names)

write.csv(files, file="genome_length.csv")

files <- lapply(file_names, function(x) read.table(x, header = FALSE)[,c(1,5)])
files <- suppressWarnings(Reduce(function(x, y) merge(x, y, by = "V1"), files))

row.names(files) <- files$V1; files <- files[,-1]
colnames(files) <- c(file_names)

write.csv(files, file="cov_bases.csv")


