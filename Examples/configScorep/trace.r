library(tidyverse)
library(dplyr)
library(magrittr)

#!/usr/bin/env Rscript
argas <- commandArgs(TRUE)

df <- read_csv("traces.csv", col_names=FALSE, col_types=cols()) %>%
        rename(Rank = X1, Start = X2, End = X3, Duration = X4, Value = X5) %>%
        select(-contains("X")) %>%
        mutate(Duration = End - Start) %>%
        mutate(Value = gsub("^P", "", Value)) %>%
        mutate(Rank = as.integer(gsub("rank", "", Rank))) %>%
        print
        
df %>%
        ggplot(aes(xmin = Start, xmax = End, ymin = Rank, ymax = Rank+0.9, fill = Value)) +
        theme_bw(base_size=16) +
        geom_rect() +
        scale_fill_brewer(palette = "Set1") +
        xlab("Tempo [segundos]") +
        ylab("Rank")
        
ggsave(argas[1], width=16, height=9)

rd <- df %>%
	group_by(Rank) %>%
	summarize(Rank.Duration = sum(Duration))
	write.csv(rd,"./Rank.Duration.csv", row.names = FALSE)
	
vd <- df %>%
	group_by(Value) %>%
	summarize(Value.Duration = sum(Duration))
	write.csv(vd,"./Value.Duration.csv", row.names = FALSE)
	
df %>% pull(End) %>% max -> makespan

rp <- df %>%
        group_by(Rank) %>%
        summarize(Rank.Percentage = sum(Duration)/makespan * 100)
        write.csv(rp,"./Rank.Percentage.csv", row.names = FALSE)

vp <- df %>%
	group_by(Value, Rank) %>%
	summarize(Value.Percentage = sum(Duration)/makespan * 100)
	write.csv(vp,"./Value.Percentage.csv", row.names = FALSE)
	
vp %>%
	ggplot(aes(x = Value.Percentage, y = Rank, fill = Value)) +
	theme_bw(base_size=16) +
        geom_bar(position="stack", stat="identity") +
        scale_fill_brewer(palette = "Set1") +
        xlab("Duração/tempo total [%]") +
        ylab("Rank")
ggsave("value.percentage.pdf", width=16, height=9)

