library('ggplot2')
genome.repeats = read.table('../read-mapping/mtb-selected-region-repeats.txt')
colnames(genome.repeats) <- c('start1', 'start2', 'length')

genome.repeats$end1 = genome.repeats$start1 + genome.repeats$length
genome.repeats$end2 = genome.repeats$start2 + genome.repeats$length
genome.repeats$middle1 = (genome.repeats$start1 + genome.repeats$end1) / 2
genome.repeats$middle2 = (genome.repeats$start2 + genome.repeats$end2) / 2

p <- ggplot(genome.repeats) + 
  geom_point(aes(x = start1, y = 1)) +
  geom_point(aes(x = start2, y = 1)) +
  geom_point(aes(x = end1, y = 1)) +
  geom_point(aes(x = end2, y = 1)) +
  geom_curve(aes(x = start1, y = 1, xend = end1, yend = 1, colour = "start 1 - end 1")) +
  geom_curve(aes(x = start2, y = 1, xend = end2, yend = 1, colour = "start 2 - end 2")) +
  geom_curve(aes(x = middle1, y = 1, xend = middle2, yend = 1, colour = "middle 1 - middle 2"))
p
