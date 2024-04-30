#! /usr/bin/R

source("renv/activate.R")

library(ComplexHeatmap)
library(cowplot)
library(data.table)
library(GAMBLR)
library(ggalluvial)
library(ggbeeswarm)
library(ggpubr)
library(glue)
library(rstatix)
library(tidyverse)

colours <- list(
    pathology = c(
        "BL" = "#926CAD",
        "FL" = "#EA8368", 
        "HGBCL-DH-BCL2" = "#7A1616",
        "DLBCL" = "#479450",
        "HGBCL-DH-BCL6" = "#002642",
        "HGBCL-NOS" = "#294936",
        "COMFL" = "#8BBC98",
        "PBL" = "#E058C0"
    ), 
    ICC_class = c(
        "BL" = "#926CAD",
        "FL" = "#EA8368", 
        "HGBCL-DH-BCL2" = "#7A1616",
        "DLBCL" = "#479450",
        "HGBCL-DH-BCL6" = "#002642",
        "HGBCL-NOS" = "#294936",
        "COMFL" = "#8BBC98",
        "PBL" = "#E058C0"
    ),
    morphology = c(
        "High-grade" = "#7A1616",
        "DLBCL" = "#479450",
        "BL" = "#926CAD",
        "FL" = "#EA8368",
        "COMFL" = "#8BBC98",
        "PBL" = "#E058C0"
    ),
    myc_partner_group = c(
        "IGH" = "#DA2C38",
        "IGK/L" = "#ec7c83",
        "recurrent non-IG" = "#226F54",
        "non-IG" = "#87C38F",
        "no break" = "#F9F7DC",
        "not found" = "black"
    ),
    myc_sv = c(POS = "#02a628", NEG = "black"),
    bcl2_sv = c(POS = "#02a628", NEG = "black"),
    bcl6_sv = c(POS = "#02a628", NEG = "black"),
    cpr_cytology = c(
        "DLBCL" = "#479450",
        "INT" = "#7a1648",
        "BLASTOID" = "#165e7a",
        "UNC" = "#7A1616"
    ),
    group = c(
        "BL" = "#926CAD",
        "FL" = "#EA8368",
        "HGBCL-DH-BCL2" = "#7A1616",
        "HGBCL-DH-BCL2 DHITsig+" = "#9c0202",
        "MYC-DHITsig+" = "#086375",
        "BCL2-DHITsig+" = "#1DD3B0",
        "DHITsig+" = "#D62828",
        "ABC-DLBCL" = "#05ACEF",
        "GCB-DLBCL" = "#F58F20",
        "PBL" = "#E058C0"
    ),
    ighc = c(
        "Variable" = "#595959",
        "Emu" = "#c41230",
        "IGHM" = "#F9BD1F",
        "IGHD" = "#B76D29",
        "IGHG3" = "#E90C8B",
        "IGHG1" = "#2CACE3",
        "IGHA1" = "#046852",
        "3'RR1" = "#c41230",
        "IGHG2" = "#5C266C",
        "IGHG4" = "#39B54B",
        "IGHE" = "#FE9003",
        "IGHA2" = "#115284",
        "3'RR2" = "#c41230",
        "NULL" = "darkgrey"
    )
)

colours$group_with_ebv <- c(
        "BL EBV-" = unname(colours$group["BL"]),
        "BL EBV+" = "#B398C6",
        "HGBCL-DH-BCL2" = unname(colours$ICC_class["HGBCL-DH-BCL2"]),
        colours$ICC_class[c(
            "FL", 
            "DLBCL"
    )])
