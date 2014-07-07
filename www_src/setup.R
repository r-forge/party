
library("pkg2html")
library("markdown")

pkg <- "party"
download.file("http://user.math.uzh.ch/hothorn/TH.bib", dest = "TH.bib")
system("cat party.bib >> TH.bib")
dest <- "../www"

if (!file.exists(dest))
    dir.create(dest)

stopifnot(file.exists(dest)) 

template <- system.file("template", package = "pkg2html")

system(paste("cp -ra", file.path(template, "*"), dest, sep = " "))

wd <- setwd(file.path(dest, "_data"))
R2yaml(pkg)
writeLines(bib2yaml(file.path(wd, "TH.bib"), 
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

file.remove("TH.bib")

x <- readLines(file.path(dest, "_data", "pkg.yml"))
x <- c(x, "headpic: /img/party.png")
writeLines(x, con = file.path(dest, "_data", "pkg.yml"))

yml <- list.files(pattern = "yml$")
sapply(yml, function(f) file.copy(f, file.path(dest, "_data"), overwrite = TRUE))

system(paste("cat site.css >> ", file.path(dest, "html/css/main.css")))


