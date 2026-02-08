library(tidyverse)
library(rstatix)
library(xlsx)
dfData = read.delim("F:/project/Code-Data-available/git/gene.txt")
dfClass = read.delim("F:/project/Code-Data-available/git/class.txt")

df = dfData %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "Sample",values_to = "value") %>%  
  left_join(dfClass,by=c("Sample" = "Sample"))
df

dfFC = df %>%
  group_by(Genes,Class) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%              
  pivot_wider(names_from = Class,values_from = mean) %>%  
  summarise(FC = VT/control)      # 计算Foldchange                    
dfFC

dfP = df %>%
  group_by(Genes) %>%
  t_test(value ~ Class,var.equal=T) # t_test 方法源于rstatix包; var.equal=T等方差
dfP

dfP_FDR = dfP %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  # p.adjust对P值FDR校正,算法选择BH
dfP_FDR


dfData %>%
  as_tibble() %>%
  left_join(dfFC) %>%
  left_join(dfP_FDR)


write.xlsx(dfFC,"F:/project/Code-Data-available/git/dfFC.xlsx")
write.xlsx(dfP_FDR,"F:/project/Code-Data-available/git/dfP_FDR.xlsx")

library(ggplot2)
library(ggrepel)

data = read.delim("F:/project/Code-Data-available/git/dfFC.xlsx",header = T)

FC = 1.5 
PValue = 0.05

data$PValue < - as.numeric(data$PValue)
data$FC < - as.numeric(data$FC)

data$sig[(-1*log10(data$PValue) < -1*log10(PValue)|data$PValue=="#N/A")|(log2(data$FC) < log2(FC))& log2(data$FC) > -log2(FC)] <- "NotSig"
data$sig[-1*log10(data$PValue) >= -1*log10(PValue) & log2(data$FC) >= log2(FC)] <- "Up"
data$sig[-1*log10(data$PValue) >= -1*log10(PValue) & log2(data$FC) <= -log2(FC)] <- "Down"

data$label=ifelse(data$Marker == 1, as.character(data$Name), '')

#绘图
ggplot(data,aes(log2(data$FC),-1*log10(data$PValue))) +    
  geom_point(aes(color = sig)) +                          
  labs(title="volcanoplot",                                
       x="log[2](FC)", 
       y="-log[10](PValue)") + 
  # scale_color_manual(values = c("red","green","blue")) + 
  geom_hline(yintercept=-log10(PValue),linetype=2)+        
  geom_vline(xintercept=c(-log2(FC),log2(FC)),linetype=2)+ 
  geom_text_repel(aes(x = log2(data$FC),                  
                      y = -1*log10(data$PValue),          
                      label=label),                       
                  max.overlaps = 10000,                    
                  size=3,                                  
                  box.padding=unit(0.5,'lines'),          
                  point.padding=unit(0.1, 'lines'), 
                  segment.color='black',      

#火山图              
install.packages("readxl")
library(readxl)

df <- read_excel("F:/project/Code-Data-available/git/rvo.xlsx", sheet = 1)
head(df)

library(ggplot2)
library(ggrepel)

ggplot(df, aes(log2FoldChange, -log10(pvalue))) + geom_point(alpha=0.6, size=2)

#数据分类
df$group<-as.factor(ifelse(df$pvalue < 0.05 & abs(df$log2FoldChange) >= 1, 
                           ifelse(df$log2FoldChange>= 1 ,'up','down'),'NS')) 
p<-ggplot(df, aes(log2FoldChange, -log10(pvalue))) + 
  geom_point(aes(color = group),alpha=0.6, size=2)+
  scale_color_manual(values = c( 'blue','grey','red'))#颜色
p

p1<-p+geom_vline(xintercept = c(-1, 1), lty=3,color = 'black', lwd=0.5) + #竖直辅助线
  geom_hline(yintercept = -log10(0.05), lty=3,color = 'black', lwd=0.5)+#垂直辅助线
  theme_bw()+
  theme(legend.title = element_blank(),
        panel.grid = element_blank())+
  labs(title="volcanoplot",
       x = 'log2 fold change',
       y = '-log10 pvalue')+ylim(0, 5)+xlim(-3, 3)
p1

df$label=ifelse(df$Marker == 1, as.character(df$gene), '')

df$label<-ifelse(df$pvalue<0.05&abs(df$log2FoldChange)>=4,"Y","N")
df$label<-ifelse(df$label == 'Y', as.character(df$gene), '')
#使用ggrepel包完成
library(ggrepel)
p1+
  geom_text_repel(aes(x = log2FoldChange,
                      y = -log10(pvalue),          
                      label=label),                       
                  max.overlaps = 10000,
                  size=3,
                  box.padding=unit(0.8,'lines'),
                  point.padding=unit(0.8, 'lines'),
                  segment.color='black',
                  segment.size=1,
                  show.legend=FALSE,
                  nudge_y=1,
                  nudge_x=-0.2
                  )

library(xlsx)
write.xlsx(df,"F:/project/Code-Data-available/git/diff.xlsx")
write.xlsx(dfP_FDR,"F:/project/Code-Data-available/git/dfP_FDR.xlsx")

library(ggplot2)
library(dplyr)
test$gene <- factor(test$gene)
head(test)
test$gene <- factor(test$gene,levels=c("Gad1","Gad2","Slc6a1","Slc6a13","Slc6a11",
                                       "Glud1","Gls","Slc17a7","Slc17a6","Slc17a8",
                                       "Vamp2","Rab3a","Snap25","Nrxn1","Pclo"
                                       ))
ggplot(test,aes(x=gene,y=max,fill=lg))+theme_classic()+
  geom_bar(stat = "identity",width=1)+scale_fill_gradient2(high="red",low="blue", midpoint=0.4,limits=c(0,1),na.value ="red" )
