#fix for recomb rate window estimates
#forces windows to be non-overlapping

ev.dir < -file.path(home.dir, "evs")

recomb <- read.table(file = file.path(ev.dir, "recomb_rate.txt"),header=TRUE)

recomb$pos1 <- recomb$pos1 + 1
recomb$pos2 <- recomb$pos2 - 1

recomb <- recomb[complete.cases(recomb), ]

write.table(recomb, file = "recomb_rate.txt", row.names = FALSE, quote = FALSE)
