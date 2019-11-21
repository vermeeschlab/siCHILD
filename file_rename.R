
# file_rename(from = "C:/Users/msc2/Desktop/rabata.txt",to = "C:/Users/msc2/Desktop/Halwa/BADMASHI/SCOP/rabata.txt")

file_rename <- function(from, to) {
	todir <- dirname(to)
	if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
	file.rename(from = from,  to = to)
}



