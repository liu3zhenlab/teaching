######################################################
# Project presentation ordering
######################################################
#setwd("/Users/liu3zhen/Documents/GitHub/teaching/PLPTH813Bioinformatis/2025/5_project")
# talks for randomization
st <- paste0("S", 1:7)
dt <- paste0("D", 1:6)
tt <- paste0("T", 1:2)
talks <- c(st, dt, tt)

# prerequested talks on May 1st
preset_talks <- c("D7", "T3", "T4")

# randomization
talk_random <- sample(talks, length(talks))
may1_talks <- sample(c(preset_talks, talk_random[1]))
talk_order <- c(may1_talks, talk_random[-1])

talk_time <- c(rep("May1am", 4),
               rep("May6am", 4),
               rep("May8am", 4),
               rep("May8pm", 6))

talks_final <- data.frame(Order = 1:length(talk_final),
                          time = talk_time,
                          talk = talk_order)

# output 
talks_list <- read.delim("talk.list.txt")
talks_out <- merge(talks_final, talks_list, by="talk")
talks_out <- talks_out[order(talks_out$Order), ]
write.table(talks_out, "BA25.project.presentation.txt", row.names=F, quote=F, sep="\t")
