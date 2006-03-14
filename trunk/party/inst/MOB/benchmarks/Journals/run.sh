
#!/bin/bash

R --no-save < JournalsSetup.R

R CMD BATCH Journals_ctree.R &

R CMD BATCH Journals_lm.R &

R CMD BATCH Journals_M5P.R &

R CMD BATCH Journals_mob.R & 

R CMD BATCH Journals_rpart.R &

R CMD BATCH GUIDE/Journals_GUIDE.R &
