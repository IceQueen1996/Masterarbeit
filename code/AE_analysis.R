Sys.setenv(JAVA_HOME="C:\\Program Files\\Java\\jdk1.8.0_191\\jre")
library(xlsx)
library(ggplot2)
library(dplyr)


# read AE data
ae <- read.xlsx("C:\\Users\\Elly\\Uni\\Masterarbeit\\data\\xlsx\\AE.xlsx", 1)
enrol <- read.xlsx("C:\\Users\\Elly\\Uni\\Masterarbeit\\data\\xlsx\\ENROL.xlsx", 1)

merged <- merge(ae, enrol, by="subjid")

######################################## AE Types Analysis #############################################

data <- data.frame(merged$aetype, merged$sex, merged$age)
colnames(data) <- c("aetype", "sex", "age")

ggplot(data, aes(x=aetype, y=age, fill=sex)) +
  geom_boxplot() +
  labs(x = "Type of Adverse Event") +
  theme_bw(base_size = 20)

###################################### TOP 10 occuring AE's ############################################
data <- as.data.frame(sort(table(merged$aecode), decreasing = TRUE)[1:10])
colnames(data) <- c("group", "value")

# Compute the position of labels
data <- data %>%
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  geom_text(aes(y = ypos, label = value), color = "black", size=6) +
  scale_fill_brewer(palette="Spectral") +
  labs(fill='Adverse Event', title = "Top 10 Occuring Adverse Events") +
  theme_void() +
  theme(legend.text=element_text(size=13), legend.title =element_text(size=15), plot.title = element_text(size = 16,hjust=0.5))

############################## TLR + MI Analysis #####################################################
data <- merged %>% filter(aecode == "TLR" | aecode == "MI")
data <- data.frame(data$aecode, data$sex, data$age)
colnames(data) <- c("aecode", "sex", "age")

ggplot(data, aes(x=aecode, y=age, fill=sex)) +
  geom_boxplot() +
  labs(x = "Adverse Event") +
  theme_bw(base_size = 20)

