########################################################################
###                   RNAseq manuscript high definition plots
#########################################################################

# Complementary script to produce high quality images 
# Requires cairo, ggarrange, 

# Cairo can be used to produce high quality images 
library(Cairo)
library(ggpubr)

publication_plots <-  list()

#########################################################################
### Fig 1: volcano and MA plot

# must be run after DESeq2 script so all data sets are properly established


tiff(filename = "Fig1.tiff",
     width = 300,
     height = 150,
     units = "mm",
     pointsize = 22, # size of text
     res = 300, # desired dpi
     type = "cairo", # specify use of Cairo
     compression = "lzw" # to reduce size
  ) 
  
# ------------ code to produce MA plot
publication_plots$ma_plot<- 
  temp_objects$maplot_combat %>% 
  arrange(padj) %>% 
  # select top n differentially expressed 
  mutate(diff = ifelse(padj < 0.3, "yes", "no"),
         top_genes = case_when(padj < 0.3 & row_number() <= 15 ~ as.character(SYMBOL),
                               TRUE ~ NA_character_)) %>%
  ggplot(aes(x = log2(baseMean),
             y = log2FoldChange,
             label = top_genes,
             color = factor(diff))) +
  geom_point(cex = 1.5) +
  scale_color_manual(values = c("grey70", "blue4")) +
  geom_label_repel(box.padding = 0.5, max.overlaps = Inf, color="blue") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.text = element_text(size = 14),
        plot.subtitle = element_text(color = "red4", face = "italic"),
        text = element_text(family = "Arial")) +
  labs(title = "MA-plot batch corrected counts",
       x = "Normalized mean counts",
       y = "Shrunk Log2 Fold Change")

# ------------ code to produce Volcano plot
publication_plots$volcano_plot<- 
  temp_objects$combat_volcano %>%
  # plot log2foldchange(x) vs -Log10pvalue (y)
  ggplot(data = .,
         mapping = aes(
           x = log2FoldChange,
           y = -log10(pvalue),
           col = diff_exp)) +
  geom_point() +
  theme_classic() +
  # red=up, green=down, grey=NA
  scale_color_manual(values = c("green", "grey50", "red")) +
  labs(title = "Volcano plot of DEG by NTM-PD status",
       col = "NTM disease") +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.001), col = "red", linetype = "dashed") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.text = element_text(size = 14),
        plot.subtitle = element_text(face = "italic", color = "red4")) +
  coord_cartesian(ylim = c(0, 7), xlim = c(-6, 6))

# code to join both plots
ggarrange(publication_plots$ma_plot, 
          publication_plots$volcano_plot,
          ncol = 2, nrow = 1, legend = "none", labels = c("A.","B."))

dev.off()

#########################################################################
### S1 Fig: Baseline CBC counts

# run it after running the exploratory_CBC.rmd

{
  tiff(filename = "S1Figure.tiff",
       width = 300,
       height = 250,
       units = "mm",
       pointsize = 22, # size of text
       res = 300, # desired dpi
       type = "cairo", # specify use of Cairo
       compression = "lzw" # to reduce size
  ) 
  ggarrange(CBCplots$wbc, CBCplots$neut,
            CBCplots$lym, CBCplots$mono,
            labels = c("A.","B.","C.","D."), 
            common.legend = T, legend = "bottom"
  ) %>% print()
  dev.off()
}

#########################################################################
### S3 Fig: Cybersort concordance with CBC counts

# must be run after CBC script
# load Cairo if necessary
library(Cairo)

tiff(filename = "S3 Fig.tiff",
     width = 400,
     height = 150,
     units = "mm",
     pointsize = 14, # size of text
     res = 300, # desired dpi
     type = "cairo", # specify use of Cairo
     compression = "lzw" # to reduce size
)
# plot of concordance of calculated percentages of cell populations
ggarrange(plotlist = comparison_plots_list[c("lym_perc","neu_perc","mon_perc")], 
          # selects only percentage plots
          ncol=3,
          nrow=1,
          labels = c("A","B","C"))

dev.off()

#########################################################################
### Figs S4: PCA plot by RNA to NTM growth window

# run it after DESeq2 script
  

tiff(filename = "S4_Fig.tiff",
     width = 200,
     height = 170,
     units = "mm",
     pointsize = 14, # size of text
     res = 300, # desired dpi
     type = "cairo", # specify use of Cairo
     compression = "lzw" # to reduce size
     )
  df <- plotPCA(temp_objects$vst_combat, intgroup="rna_growth_window", returnData=T)
  percentVar <- round(100 * attr(df, "percentVar"))
  ggplot(df, aes(PC1, PC2, color = rna_growth_window)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() +
    theme_classic() +
    labs(color = "Interval between NTM\n growth and RNA sample") +
    theme(panel.grid.major = element_line(color = "grey76",
                                          size = 0.5,
                                          linetype = 2),
          legend.position = "bottom",
          legend.title = element_text(hjust = 0.5)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
dev.off()

#########################################################################
### Figs S5: Mean expression of genes found in Cowman cohort

### ----------------- Top 15 DEG in Cowman et al. (13/15)
# after DESeq2 script

tiff(filename = "S4Fig.tiff",
     width = 300,
     height = 250,
     units = "mm",
     pointsize = 14, # size of text
     res = 300, # desired dpi
     type = "cairo", # specify use of Cairo
     compression = "lzw" # to reduce size
)
# must be run after DESeq2 script
ggarrange(plotlist = plot_candidates[c(1:14)], 
          ncol = 4, nrow = 4, legend = "none") %>% 
  annotate_figure(left = text_grob("Normalized counts", rot = 90))

dev.off()

