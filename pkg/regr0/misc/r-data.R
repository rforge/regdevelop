##- d.blast <- d.blast[,-3]
##- write.table(d.blast, file="~/data/blast.dat")
d.blast <- read.table("~/data/blast.dat",header=T)
tit(d.blast) <- "Blasting for Tunnel Excavation"
doc(d.blast) <-
  c("Blasting for tunnel excavation: measurements of tremor in nearby house basements",
 "location: Schaffhausen, Switzerland","source: Basler and Hofmann, Zurich")

d.surveyenvir <- d.umweltumf
names(d.surveyenvir) <- c("age","sex","education","location","townsize",
 "party","disturbance","gov","responsability","weight")
levels(d.surveyenvir$sex) <- c("m","f")
levels(d.surveyenvir$education) <- c("no.training","apprentiship","no.degree",
 "college","uni")
levels(d.surveyenvir$disturbance) <- c("none","slight","clear","severe")
d.surveyenvir$disturbance <- ordered(d.surveyenvir$disturbance)
levels(d.surveyenvir$gov) <- c("does.enough","not.enough")
levels(d.surveyenvir$responsability) <- c("individuals","government","both")
tit(d.surveyenvir) <- "Survey on Environment"
doc(d.surveyenvir) <- c("Survey on attitudes towards environmental problems",
    "source: Umweltschutz im Privatbereich.",
    "Erhebung des EMNID, Zentralarchiv fuer empirische Sozialforschung der Universitaet Koeln")

##- package.skeleton("a",list=c("d.blast","d.surveyenvir","d.rehab"),path="..",
##-                  force=T)
## 
## --------------------------------------
source("~/data/d.rehab.dat")
dd <- d.rehab
levels(dd$SPITAL) <- LETTERS[1:length(levels(dd$SPITAL))]
levels(dd$ARZT) <- c(LETTERS[1:(length(levels(dd$ARZT))-1)],"Other")
dv <- read.table("~/data/d.rehab.txt", sep="&", header=F, as.is=1:3)
t.j <- pmatch(names(dd),dv[,1])
t.j[2] <- 2
d.rehab <- dd[,!is.na(t.j)][,order(nna(t.j))]
attr(d.rehab,"variables") <- dv
doc(d.rehab) <- c("Data on rehabilitation of rheumatic arthritis patients")
tit(d.rehab) <- "Rehabilitation of Rheumatic Arthritis"

## --------------------------------------
# dump("d.fossiles","/u/stahel/data/d-fossiles.R")
source("/u/stahel/data/d-fossiles.R")
tit(d.fossiles) <- "Shapes of Gephyrocapsa Shells"
doc(d.fossiles) <-
  c("Shapes of Gephyrocapsa shells and environmental characteristics",
    "Source: Joerg Bollmann, Jorijntje Henderiks and Bernhard Brabec",
    "Geological Institute, ETH Zurich",
    "Marine Micropaleontology 29, 319-350 (1997)")

package.skeleton("a",list=c("d.blast"),path="..", # ,"d.surveyenvir","d.rehab"
                 force=T)

