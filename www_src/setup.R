
library("pkg2html")
library("markdown")

pkg <- "party"
download.file("http://user.math.uzh.ch/hothorn/TH.bib", dest = "TH.bib")
dest <- "html"

if (!file.exists(dest))
    dir.create(dest)

stopifnot(file.exists(dest)) 
system("rm -rf www/_posts/*")

template <- system.file("template", package = "pkg2html")

system(paste("cp -ra", file.path(template, "*"), dest, sep = " "))

wd <- setwd(file.path(dest, "_data"))
R2yaml(pkg)
writeLines(bib2yaml(file.path(wd, "TH.bib"), 
           c("Hothorn:2006:JCGS", "Strobl:2006:BMC-Bioinformatics")), con = "cites.yml")

setwd(wd)
setwd(file.path(dest, "_posts"))
NEWS2md(pkg)

setwd(wd)

Rmd <- list.files(pattern = "Rmd$")

for (f in Rmd)
    writeLines(Rmd2html(f), con = file.path(dest, gsub("Rmd$", "html", f)))

file.remove("TH.bib")

setwd(dest)
