
#script to visualize gene expression data (GSE183947)

#load libraries
library(tidyverse)
library(ggplot2)

#data
#dat.long to be used 

#basic format for ggplot
#pipe operator ' + '

#ggplot(data, aes(x=variable1, y= variable2)) +
#  geom_bar()

#barplot
#to see the dataset you loaded onto your environment use the following functions-
#GSE183947_long_format %>%
#  head()

#use the 'filter' function to use the only variables 
#you need to visualize the data
#use the 'fill' function to color the variables i.e., - tissue type 
#GSE183947_long_format %>%
#  filter(gene == 'BRCA1') %>%
#  head()

#barplot
GSE183947_long_format %>%
   filter(gene == 'BRCA1') %>%
   ggplot(.,aes(x = samples, y = FPKM, fill = tissue)) +
   geom_col()
# to see the overlapping parts of the graph,
# we can change opacity parameter by using 'alpha' function

#density
GSE183947_long_format %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3)

##geom_boxplot()

GSE183947_long_format %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  geom_violin()

#to see the expression difference between two diff genes
#we use a scatter plot!

# we can filter out the two genes we want to compare-

#  dat.long %>%
#  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
 
# to use the wide format/ convert format 
# use the 'spread' function

GSE183947_long_format %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  head
  
#to plot this data frame use 'geom_point' for plotting as scatter plot
# to fit a straight line use 'geom_smooth' function
# to not to use any 'confidendnese intervals'  use
# method ='lm', se = FALSE

GSE183947_long_format %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)


# heatmaps are better choice to see the diff expressions among  more than two genes
# use the  %in%  function to use just the selected genes 
# that are included in the variable called 'genes.of.interest'

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
GSE183947_long_format %>%
  filter(gene %in% genes.of.interest) %>%
  head()

# to create a heatmap use 'geom_tile()' function

GSE183947_long_format %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile()

# to change the color of the heatmap use 'scale_fill_gradient' function

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

GSE183947_long_format %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')

#to save any plot you created from the gene exp analysis:

# to do that first 
# rename the variable, for example 'p'
# run the the whole script under the new variable name
# use 'ggsave' function 
# to save as pdf use following linescript -

genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
p<-GSE183947_long_format %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')


ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8)

##plotting multiple graph onto the same image using 'ggpbur'
install.packages("ggpubr")
library(ggpubr)

p<-GSE183947_long_format %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')
q<-GSE183947_long_format %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')
r<-GSE183947_long_format %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)
s<-GSE183947_long_format %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  geom_violin()

ggarrange(p, q, r, s ,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("ggpubrsave.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)




