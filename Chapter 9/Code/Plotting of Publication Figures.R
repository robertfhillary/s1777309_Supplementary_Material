setwd("./Epigenetic Clocks WHO")

## Installing Requisite Packages

if(!require(ggplot2)){
  install.packages("ggplot2")
}

if(!require(cowplot)){
  install.packages("cowplot")
}

library(ggplot2)
library(cowplot)


## Create function for element textbox - to give individual plots in Figure 1 headers 

element_textbox <- function(...) {
  el <- element_text(...)
  class(el) <- c("element_textbox", class(el))
  el
}

element_grob.element_textbox <- function(element, ...) {
  text_grob <- NextMethod()
  rect_grob <- element_grob(calc_element("strip.background", theme_bw()))
  
  ggplot2:::absoluteGrob(
    grid::gList(
      element_grob(calc_element("strip.background", theme_bw())),
      text_grob
    ),
    height = grid::grobHeight(text_grob), 
    width = grid::unit(1, "npc")
  )
}

## Composite Plot for Figure 1 - Continuous, Categorical and All-Cause Mortality merged into one plot

x = read.csv("Fig1.csv")

## Define order of variables by magnitude of effect sizes 

x$Clock2 <- x$Clock
levels(x$Clock2) = c(" ", "  ", "   ", "    ")
x$TraitVar <- paste0(x$Trait,x$Clock2)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$test))]))
x$facet = factor(x$Group, levels = c("Continuous", "Disease", "All-Cause Mortality"))

## Continuous variable portion of plot 

graph_legend <- ggplot(x1,aes(y=test, x=TraitVar, colour = Clock)) + 
  geom_point() + theme(legend.title=element_text(hjust=0.5)) + guides(colour = guide_legend(override.aes = list(shape = 17, size = 3))) 


x1 <- x[which(x$Group %in% "Continuous"), ]
graph1 <- ggplot(x1,aes(y=test, x=TraitVar, colour = Clock)) + 
  geom_point(position=position_dodge(width=0.5), shape = 17, size = 3)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Standardised Beta [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 0, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() + labs(title = "Continuous") + 
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  )


## Categorical variable portion of plot 

x2 = x[which(x$Group %in% "Disease"),]

graph2 <- ggplot(x2,aes(y=test, x=TraitVar, colour = Clock)) + 
  geom_point(position=position_dodge(width=0.5), shape = 17, size = 3)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Odds Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+ scale_y_continuous(limits = c(0,6.5), breaks = c(0,1,2,3,4,5,6)) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() +  labs(title = "Disease") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  )

graph2

## All-cause mortality portion of plot 

x3 = x[which(x$Group %in% "All-Cause Mortality"),]

graph3 <- ggplot(x3,aes(y=test, x=TraitVar, colour = Clock)) + 
  geom_point(position=position_dodge(width=0.5), shape = 17, size = 3)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+ scale_y_continuous(limits = c(0,6.5), breaks = c(0,1,2,3,4,5,6)) +
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "none",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() +  labs(title = "Mortality") +
  theme(
    plot.title = element_textbox(
      hjust = 0.5, margin = margin(t = 5, b = 5)
    )
  )

graph3

## Create combined legend for three individual plots in Figure 1 

legend <- get_legend(graph_legend)

## Assemble all components of plot together

plot = ggdraw(plot_grid(graph1, plot_grid(graph2, graph3, nrow = 2, rel_widths = c(1/4,1/4)), legend, ncol = 3, rel_widths = c(0.75,1/2,0.3)))

plot

## Plot for Figure 2 - Split by Epigenetic Clock 

x = read.csv("Figure_2.csv")

## Define order of variables by magnitude of effect sizes 

x$Clock2 <- x$Clock
levels(x$Clock2) = c(" ", "  ", "    ")
x$TraitVar <- paste0(x$Trait,x$Clock2)
x$TraitVar = factor(x$TraitVar, levels=unique(x$TraitVar[rev(order(x$test))]))


graph1 <- ggplot(x,aes(y=test, x=TraitVar, colour = Clock)) + 
  geom_point(position=position_dodge(width=0.5), shape = 17, size = 3)+
  geom_errorbar(aes(ymin = LCI, ymax = UCI),
                position = position_dodge(0.5), width = 0.05,
                colour = "black")+
  ylab("Hazard Ratio [95% Confidence Interval]")+ xlab ("")+
  geom_hline(yintercept = 1, linetype = "dotted")+
  theme(axis.text.x = element_text(size = 12, vjust = 0.5), legend.position = "right",
        plot.title = element_text(size = 12))+ theme(legend.title = element_text(hjust = 0.5)) +
  coord_flip() 

graph1 + facet_wrap(vars(Clock), scales = "free_y")

          