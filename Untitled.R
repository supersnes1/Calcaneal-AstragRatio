# lab 1b  names <- read.csv("/Users/emdoughty/Downloads/LS7B-1_1bRoster.csv")
# lab 1i  names <- read.csv("/Users/emdoughty/Downloads/LS7B-1_1iRoster.csv")
# lab 1L  
names <- read.csv("/Users/emdoughty/Downloads/LS7B-1_1LRoster.csv")

names
gardens <- c("Hawaiian", "North Conifer", "Pollination", "Fern", "Mediterranean", " Desert", "Cycad")

groupAll <- seq(1,nrow(names),1)
groupList <- list()

for(xx in 1:6)
{
  if(length(groupAll)-4 < 0) {print(names[groupAll,]); break}
  sel <- sample(groupAll, 4)

  groupList[[xx]] <- names[sel,]

  groupAll <- groupAll[!groupAll %in% sel]
}

sample(gardens, 6)

