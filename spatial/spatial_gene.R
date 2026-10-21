rm(list = ls())
setwd("\\\\EMTAB13530")

library(Seurat)
library(GSVA)
library(tidyverse)
library(scales)

st_obj = Load10X_Spatial(data.dir = "\\EMTAB13530\\P10_T1\\")
st_obj <- NormalizeData(st_obj)

p1 <- SpatialFeaturePlot(st_obj, features = "MAOB") +
  scale_fill_viridis_c(
    labels = number_format(accuracy = 0.1)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("P10_T1MAOB")
p1
ggsave(filename = "/P10MAOB.pdf", p1, height = 4, width = 4.5, dpi = 600)

p2 <- SpatialFeaturePlot(st_obj, features = "SLC47A1") +
  scale_fill_viridis_c(
    labels = number_format(accuracy = 0.1)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("P10_T1SLC47A1")
p2
ggsave(filename = "/P10SLC47A1.pdf", p2, height = 4, width = 4.5, dpi = 600)

p3 <- SpatialFeaturePlot(st_obj, features = "AKAP12") +
  scale_fill_viridis_c(
    labels = number_format(accuracy = 0.1)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("P10_T1AKAP12")
p3
ggsave(filename = "/P10AKAP12.pdf", p3, height = 4, width = 4.5, dpi = 600)

p4 <- SpatialFeaturePlot(st_obj, features = "SLCO1B3") +
  scale_fill_viridis_c(
    labels = number_format(accuracy = 0.1)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("P10_T1SLCO1B3")
p4
ggsave(filename = "/P10SLCO1B3.pdf", p4, height = 4, width = 4.5, dpi = 600)

p5 <- SpatialFeaturePlot(st_obj, features = "LDHA") +
  scale_fill_viridis_c(
    labels = number_format(accuracy = 0.1)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("P10_T1LDHA")
p5
ggsave(filename = "/P10LDHA.pdf", p5, height = 4, width = 4.5, dpi = 600)

p6 <- SpatialFeaturePlot(st_obj, features = "ALDOA") +
  scale_fill_viridis_c(
    labels = number_format(accuracy = 0.1)
  ) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.key.size = unit(0.4, "cm")
  ) +
  ggtitle("P10_T1ALDOA")
p6
ggsave(filename = "/P10ALDOA.pdf", p6, height = 4, width = 4.5, dpi = 600)
