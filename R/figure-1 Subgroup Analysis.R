#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
# figure-1 Subgroup Analysis.R
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

library(tidyverse)
library(ggplot2)
library(reshape2)     # melt()
library(RColorBrewer) # color palettes
library(grid)         # arrange plots
library(gridExtra)    # arrange plots
library(cowplot)      # common legend

setwd('UCSF-R/ipd-ma-cd2 GitHub')

#------------------------------------------------------------------------------#

# Load data (pre-calculated drug preferences
# SE calculated using 95% CI boostrapped with n=10,000 sim)
crohns <- read.csv("data/treatment_preferences-n10000.csv") 

#------------------------------------------------------------------------------#
# Subgroup Populations
#------------------------------------------------------------------------------#

# TNFI > (IL12, INTG)
prefer_tnfi <- crohns %>% 
  filter(p12_ohe == 1 & drug1 %in% c('tnfi')) 

# (TNFI, IL12) > INTG
prefer_tnfi_il12 <- crohns %>% 
  filter(p12_ohe == 0 & p23_ohe == 1 
         & drug1 %in% c('tnfi','il12') 
         & drug2 %in% c('tnfi','il12'))

# IL12 > (TNFI, INTG)
prefer_il12 <- crohns %>% 
  filter(p12_ohe == 1 & drug1 %in% c('il12'))

drugClassSummaries <- function(subpopulation) {
  # returns N, placebo reduction average, and drug reduction average 
  # for given subpopulation - results used for bar plots
  
  # Anti-TNF summary
  tnfi_summary <- subpopulation %>% summarise(
    N = n(), 
    subgroup = 'Anti-TNF',
    placebo.avg = mean(plac.attrib), 
    drug.avg = mean(tnfi.attrib))    # change here
  
  # Anti-IL12/23
  il12_summary <- subpopulation %>% summarise(
    N = n(), 
    subgroup = 'Anti-IL-12-23',
    placebo.avg = mean(plac.attrib), 
    drug.avg = mean(il12.attrib))    # change here
  
  # Anti-Integrin
  intg_summary <- subpopulation %>% summarise(
    N = n(), 
    subgroup = 'Anti-Integrin',
    placebo.avg = mean(plac.attrib), 
    drug.avg = mean(intg.attrib))    # change here
  
  results <- rbind(data.frame(), tnfi_summary, il12_summary, intg_summary) %>% 
    # pivot longer to format for bar plots
    pivot_longer(cols = ends_with("avg"), names_to = "group", values_to = "avg")
  
  return( results )
}

# summarize drug class results - for bar plot
subgroup1_tnfi <- drugClassSummaries(prefer_tnfi) 
subgroup2_both <- drugClassSummaries(prefer_tnfi_il12) 
subgroup3_il12 <- drugClassSummaries(prefer_il12)

#------------------------------------------------------------------------------#
# Plots
#------------------------------------------------------------------------------#

# plot titles
titles <- c(
  sprintf("Prefer Anti-TNF Only (N = %d)",             nrow(prefer_tnfi)),
  sprintf("Prefer Anti-TNF or Anti-IL-12/23 (N = %d)", nrow(prefer_tnfi_il12)),
  sprintf("Prefer Anti-IL-12/23 Only (N = %d)",        nrow(prefer_il12))
)

# subgroup labels
labels <- c('Prefer Anti-TNF', 
            'Prefer Anti-TNF or Anti-IL-12/23', 
            'Prefer Anti-IL-12/23')

# subgroup colors
colors <- brewer.pal(3, "Set1")

#------------------------------------------------------------------------------#
# Row 1 - Subgroup Bar Plots
#------------------------------------------------------------------------------#

# bar plot function
subgroupBarPlot <- function(subgroup_df, title = '', ylab = '', xlab = '', drug_color = 'black') {
  
  # add placebo.avg + drug.avg - print total at top of bar plot
  totals <- subgroup_df %>% 
    group_by(subgroup) %>% 
    summarize(total = round(sum(avg),1))
  
  # bar plot
  p <- subgroup_df %>% 
    group_by(subgroup) %>% 
    
    # add vars for stacked bar label (P, D) and height (middle of bar) 
    mutate(bar_label = ifelse(group == 'placebo.avg', 'P' ,'D')) %>%
    mutate(label_y = cumsum(avg) - 0.5 * avg) %>% 
    
    ggplot(aes(x = reorder(subgroup, -avg), y = avg, fill = group)) + 
    # bar plot
    geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
    
    # add total average cdai reduction to top of each bar
    geom_text(data = totals, aes(x = subgroup, y = total, label = total, fill = NULL), nudge_y = 10) + 
    
    # add bar label (P, D) in middle of stacked bar
    geom_text(aes(y = label_y, label = bar_label), size = 4, color = 'grey') +
    
    # add title, labels
    ggtitle(title) + xlab(xlab) + ylab(ylab) +
    
    # manually fill bar colors 
    scale_fill_manual(
      name = NULL,
      labels = c("Drug Attributable", "Placebo Attributable"),
      values = c(drug_color,"white")) +
    
    # scale y axis - constant throughout all bar plots
    scale_y_continuous(limits = c(0,150)) + 
    
    # theme specifications
    theme_bw() + 
    theme(
      plot.title = element_text(size = 12),
      panel.grid.major.x = element_blank()) + 
    theme(legend.position = "none")
  
  p
}

# row 1 (r1) bar plots
r1.tnfi <- subgroupBarPlot(subgroup1_tnfi, 
                           title      = titles[1],
                           drug_color = colors[1]) + 
  geom_vline(xintercept = 1.5, linewidth = 1.5) # separate preferred vs not-preferred

r1.both <- subgroupBarPlot(subgroup2_both, 
                           title      = titles[2],
                           drug_color = colors[2]) + 
  geom_vline(xintercept = 2.5, linewidth = 1.5) # separate preferred vs not-preferred

r1.il12 <- subgroupBarPlot(subgroup3_il12, 
                           title      = titles[3],
                           drug_color = colors[3]) + 
  geom_vline(xintercept = 1.5, linewidth = 1.5) # separate preferred vs not-preferred

#------------------------------------------------------------------------------#
# Row 2 - Binary Covariates
#------------------------------------------------------------------------------#

# combine subgroups, add label column
combined_df <- rbind( 
  prefer_tnfi %>% mutate(label=labels[1]),
  prefer_tnfi_il12 %>% mutate(label=labels[2]),
  prefer_il12 %>% mutate(label=labels[3])
) %>% 
  # select covariates (continuous, binary) + label
  dplyr::select(CDAI_baseline:CRP, HxOfTNFi:Ileal, label) %>% 
  mutate(label = factor(label, levels = labels))

r2.binary <- combined_df %>% 
  
  # select binary covariates, label
  dplyr::select(HxOfTNFi:Ileal, label) %>% 
  
  # format (melt) data for bar plots
  group_by(label) %>%
  summarize('Female'                  = 100 * ( 1 - sum(Sex_Male)/n() ),
            "Prior Anti-TNF Use"      = 100 * ( sum(HxOfTNFi)/n() ),
            "Immunomodulator Use"     = 100 * ( sum(ImmUse)/n() ),
            "Steroid Use"             = 100 * ( sum(SteroidUse)/n() ),
            "Disease Location: Ileum" = 100 * ( sum(Ileal)/n() ) ) %>% 
  reshape2::melt() %>% 
  
  # ggplot
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  
  # axis labels
  xlab("") + ylab("Percentage (%)") +
  
  # theme specifications  
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank()) + 
  
  # dashed line at 50%
  geom_hline(yintercept=50, linetype = "dashed")

#------------------------------------------------------------------------------#
# Row 3 - Continuous Covariates
#------------------------------------------------------------------------------#

createViolinPlot <- function(combined_df, var = '', 
                             title='', xlab='', ylab=''){
  
  violin <- combined_df %>% 
    # ggplot
    ggplot(aes(x = label, y = .data[[var]], fill = label)) + 
    
    # violin plot
    geom_violin(alpha = 0.8) +
    
    # box-quantile plot (within violin plot)
    geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
    
    # title, labels
    ggtitle(title) + xlab(xlab) + ylab(ylab) + 
    
    # colors
    scale_fill_manual(values = colors) +
    
    # theme specifications
    theme_bw() + 
    theme(
      legend.position    = "none", 
      plot.title         = element_text(size = 12),
      axis.text.x        = element_blank(), 
      panel.grid.major.x = element_blank()) + 
    
    # rotate x axis labels 90 degrees
    scale_x_discrete(guide = guide_axis(angle = 90))
  
  violin
  
}


r3.cdai <- createViolinPlot(combined_df, 
                            var   = 'CDAI_baseline', 
                            title = 'CDAI Baseline')

r3.age <- createViolinPlot(combined_df, 
                           var   = 'Age', 
                           title = 'Age (yrs)')

r3.bmi <- createViolinPlot(combined_df, 
                           var   = 'BMI', 
                           title = 'BMI (kg/m2)')

r3.crp <- createViolinPlot(combined_df, 
                           var   = 'CRP', 
                           title = 'CRP (mg/L)')

#------------------------------------------------------------------------------#
# Combine Plots
#------------------------------------------------------------------------------#

#----------------------------- 
# Legend
#-----------------------------

# dummy plot - just need legend details
p.legend <- data.frame(a = seq(1, length(labels)),
                       b = factor(labels, levels = labels) ) %>%
  ggplot(aes(y = a, fill = b)) +
  geom_bar() +
  scale_fill_manual(
    name = NULL,
    values = c(colors) ) + 
  theme(legend.position = "bottom")

# extract legend
overall.legend <- cowplot::get_legend(p.legend)

#----------------------------- 
# Title
#-----------------------------

overall.title <- ggdraw() + 
  draw_label(
    "Subgroup Analysis",
    fontface = 'bold',
    x = 0.01, 
    hjust = 0,
    size = 20,
  ) +
  theme(plot.margin = margin(0, 0, 0, 0))

#----------------------------- 
# Plots by row
#-----------------------------

# row 1
overall.r1 <- grid.arrange(
  # concat plots into row
  cowplot::plot_grid(r1.tnfi, r1.both, r1.il12, ncol = 3),
  # add title, axis titles
  # top = textGrob("Comparing Treatment Recommendations by Subgroups", 
  #               gp=gpar(fontface='bold',  cex = 1.5)),
  left = textGrob("Average CDAI Reduction", gp=gpar(fontsize=11), rot=90, vjust = 1.5), 
  bottom = textGrob("Treatment Received", gp=gpar(fontsize=11), vjust = -1) )

# row 2
overall.r2 <- grid.arrange(
  # concat plots into row
  cowplot::plot_grid(r2.binary, ncol=1))
  # add title
  # top = textGrob("Binary Covariate Percentages", gp=gpar(fontface='bold',  cex = 1.5)) )

# row 3
overall.r3 <- grid.arrange(
  # concat plots into row
  cowplot::plot_grid(r3.cdai, r3.age, r3.bmi, r3.crp, ncol=4),
  # add title, axis title
  #top = textGrob("Continuous Covariate Distributions", gp=gpar(fontface='bold',  cex = 1.5)), 
  bottom = textGrob("Subgroups", gp=gpar(fontsize=11), vjust = -1) )

#----------------------------- 
# Final Plot
#-----------------------------

# combined plot
subgroup.final <- cowplot::plot_grid(
  overall.r1, overall.r2, overall.r3, overall.legend, 
  ncol = 1, nrow = 4, rel_heights = c(1, 1, 1, 0.1),
  labels = c('a','b')
)

#------------------------------------------------------------------------------#

# save 
ggsave(filename = 'images/figure-2.pdf', 
       plot = subgroup.final, 
       width = 12, height = 9, units='in', dpi=300)

ggsave(filename = 'images/figure-2.jpeg', 
       plot = subgroup.final, 
       width = 12, height = 9, units='in', dpi=300)

ggsave(filename = 'images/figure-2.png', 
       plot = subgroup.final, 
       width = 12, height = 9, units='in', dpi=300)

ggsave(filename = 'images/figure-2.tiff', 
       plot = subgroup.final, 
       width = 12, height = 9, units='in', dpi=300)

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#