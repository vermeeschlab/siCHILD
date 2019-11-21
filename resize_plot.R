rm(list=ls())
unlink(".RData")


# /uz/data/hydra/shared_app/apps/R/3.4.0/bin/R --save --restore --args /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/ resize_plot_list.txt < ~/git_gcap/sichild/resize_plot.R 

args <- commandArgs(TRUE)
wd <- args[1]	# /uz/data/hydra/genomicscore/raw/HiSeqComputed/new/diagnostics/uz/cme_pgd/gcpu/samples/
dir_list <- args[2]		# resize_plot_list.txt


dl <- read.table(paste0(wd,dir_list),stringsAsFactors=F)
cmd <- c()
for (i in 1 : nrow(dl)) {
	fn <- list.files(paste0(wd,dl[i,]),pattern = ".jpg$", recursive = TRUE)
	if(!identical(fn, character(0)) ){
		for (j in 1:length(fn)){
			fn_new <- paste0(fn[j],".resize.jpg")
			cmd1 <- paste0("/uz/data/hydra/shared_app/apps/imagemagick/7.0.7-27/bin/convert -flatten -density 250 -quality 100 -resize 10% ",wd,dl[i,],"/",fn[j]," ",wd,dl[i,],"/",fn_new)
			cmd2 <- paste0("rm ",wd,dl[i,],"/",fn[j])
			cmd <- c(cmd1,cmd2,cmd)
		}
	}
}


out_fn <- paste0(wd,dir_list,".cmd.txt")
write.table(cmd,file = out_fn, quote=F,col.names=F, row.names=F)


