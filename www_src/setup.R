
library("pkg2html")
library("markdown")

source("config.R")

if (!file.exists(dest))
    dir.create(dest)

stopifnot(file.exists(dest)) 

template <- system.file("template", package = "pkg2html")

system(paste("cp -ra", file.path(template, "*"), dest, sep = " "))

wd <- setwd(file.path(dest, "_data"))
R2yaml(pkg)
writeLines(bib2yaml(file.path(wd, bib), 
           c("Hothorn:2006:JCGS", "Hothorn:2006:Biostatistics:16344280",
             "Strobl:2006:BMC-Bioinformatics", "Zeileis+Hothorn+Hornik:2008",
             "Strobl_Boulesteix_Kneib_2008", "Strobl_Malley_Tutz_2009")), 
           con = "cites.yml")

setwd(wd)
setwd(file.path(dest, "_posts"))
NEWS2md(pkg)

setwd(wd)

Rmd <- list.files(pattern = "Rmd$")

for (f in Rmd)
    writeLines(Rmd2html(f), con = file.path(dest, gsub("Rmd$", "html", f)))

x <- readLines(file.path(dest, "_data", "pkg.yml"))
x <- c(x, paste("headpic: /img/", pkg, ".png", sep = ""))
writeLines(x, con = file.path(dest, "_data", "pkg.yml"))

file.copy("_config.yml", dest, overwrite = TRUE)

yml <- list.files(pattern = "yml$")
yml <- yml[-grep("^_", yml)]
sapply(yml, function(f) file.copy(f, file.path(dest, "_data"), overwrite = TRUE))

system(paste("cat site.css >> ", file.path(dest, "css", "main.css")))

system(paste("cp ", file.path(publish, "img", "*"), file.path(dest, "img")))

