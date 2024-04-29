
pkg <- "party"
dest <- "html"
publish <- "../www"

download.file("http://user.math.uzh.ch/hothorn/bib/TH.bib", dest = "TH.bib")
system("cat party.bib >> TH.bib")
bib <- "TH.bib"