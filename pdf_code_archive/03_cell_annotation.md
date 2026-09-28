# Cell-type annotation and UMAP

PDF transcription only. The PDF may wrap code or omit glyphs;
use the curated scripts for executable analysis.

## Source PDF page 29

```text
Idents(merged) <- factor(Idents(merged), levels = as.character(0:14))

genes <- c(
"TraesCS7A02G088700","TraesCS4A02G388000","TraesCS7D02G084100",
"TraesCS4B02G240900","TraesCS4D02G240700","TraesCS4A02G063800",
"TraesCS4A02G016400","TraesCS4B02G287800","TraesCS4D02G286500",
"TraesCS6B02G050700","TraesCS6D02G041700","TraesCS6A02G036100",
"TraesCS4D02G032000","TraesCS4B02G033300","TraesCS4A02G279200","TraesCS4D02G031800",
"TraesCS4D02G226100",
"TraesCS4B02G297500","TraesCS4D02G296400","TraesCS4A02G007400",
"TraesCS2A02G474000","TraesCS2D02G473700",
"TraesCS1A02G341300","TraesCS1D02G343400",
"TraesCS2A02G350700","TraesCS2D02G348800","TraesCS6D02G199500","TraesCS6A02G216100",
"TraesCS1D02G194000","TraesCS1A02G189900","TraesCS1B02G192000",
"TraesCS4A02G322200","TraesCS5D02G552200",
"TraesCS2A02G066800","TraesCS5A02G165400","TraesCS5D02G169900",
"TraesCS7A02G276400","TraesCS1B02G317500","TraesCS1A02G403300",
"TraesCS5D02G329200","TraesCS5B02G322900","TraesCS5A02G322500",
"TraesCS7B02G354800","TraesCS3A02G230000","TraesCS3B02G259300",
"TraesCS6D02G314300","TraesCS6A02G334800","TraesCS6B02G365200"
)

gene_symbols <- c(
"SULTR3;4","SULTR3;4","SULTR3;4",
"TaGSr","TaGSr","TaGSr",
"TaSUT1","TaSUT1","TaSUT1",
"CPIII","CPIII","CPIII",
"gl-OXO","gl-OXO","gl-OXO","gl-OXO",
"HIC",
"FDH","FDH","FDH",
"ATML1","ATML1",
"DCR","DCR",
"EP3","EP3","EP3","EP3",
"ALMT12","ALMT12","ALMT12",
"MYB60","MYB60",
"RBCS","RBCS","RBCS",
"CAB3","CAB3","CAB3",
"LHCB2.1","LHCB2.1","LHCB2.1",
"CA1","CA1","CA1",
"AOC2","AOC2","AOC2"
)
```

## Source PDF page 30

```text

symbol2cell <- c(
"SULTR3;4" = "Vascular Bundle",
"TaGSr" = "Vascular Bundle",
"TaSUT1" = "Phloem",
"CPIII" = "Procambium",
"gl-OXO" = "Inner Sheath",
"HIC"   = "Guard Cells/Epidermis",
"FDH"   = "Epidermis",
"ATML1" = "Epidermis",
"DCR"   = "Epidermis",
"EP3"   = "Epidermis",
"ALMT12" = "Guard Cells",
"MYB60" = "Guard Cells",
"RBCS"  = "Mesophyll",
"CAB3"  = "Mesophyll",
"LHCB2.1" = "Mesophyll",
"CA1"   = "Mesophyll",
"AOC2"  = "Mesophyll"
)

dp <- DotPlot(merged, features = genes)
d <- dp$data
d$feat <- as.character(d$features.plot)
pos_map <- setNames(seq_along(genes), genes)
d$ypos <- pos_map[d$feat]
d$label <- gene_symbols[match(d$feat, genes)]

group_df <- data.frame(symbol = unique(gene_symbols), stringsAsFactors= FALSE)
group_df$celltype <- symbol2cell[group_df$symbol]
group_df$minpos <- sapply(group_df$symbol, function(s) min(which(gene_symbols == s)))
group_df$maxpos <- sapply(group_df$symbol, function(s) max(which(gene_symbols == s)))
group_df$mid <- (group_df$minpos + group_df$maxpos)/2
d$idx <- as.numeric(as.character(d$id))

p <- ggplot(d, aes(x = idx, y = ypos, size = pct.exp, color = avg.exp)) +
geom_point(shape = 16) +
scale_y_reverse(
breaks = NULL,
labels = NULL,
expand = c(0.01, 0)
) +
```

## Source PDF page 31

```text
scale_x_continuous(
breaks = 0:14,
labels = 0:14,
expand = c(0, 0)
) +
scale_size(range = c(1, 9), guide = guide_legend(reverse = TRUE)) +
scale_color_distiller(palette = "PuRd", direction = 1) +
theme_minimal(base_size = 12) +
theme(
axis.title = element_blank(),
axis.text.x = element_text(size = 11),
axis.text.y = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.margin = unit(c(0.5, 3, 0.5, 7), "cm"),
axis.line.y = element_line(color = "black"),
axis.line.x = element_line(color = "black"),
axis.ticks.y = element_blank()
) +
coord_cartesian(
xlim = c(-0.5, 14.5),
clip = "off"
) +
geom_segment(
data = group_df,
aes(x = -2.5, xend = -2.5, y = minpos - 0.4, yend = maxpos + 0.4),
inherit.aes = FALSE,
linewidth = 0.6
) +
geom_segment(
data = group_df,
aes(x = -2.55, xend = -2.45, y = minpos - 0.4, yend = minpos - 0.4),
inherit.aes = FALSE,
linewidth = 0.6
) +
geom_segment(
data = group_df,
aes(x = -2.55, xend = -2.45, y = maxpos + 0.4, yend = maxpos + 0.4),
inherit.aes = FALSE,
linewidth = 0.6
) +
geom_text(
data = group_df,
aes(x = -2.7, y = mid, label = celltype),
```

## Source PDF page 32

```text
hjust = 1,
vjust = 0.5,
size = 4,
inherit.aes = FALSE
) +
geom_text(
data = data.frame(y = seq_along(genes), symbol = gene_symbols),
aes(x = -0.8, y = y, label = symbol),
hjust = 1,
vjust = 0.5,
size = 3.5,
inherit.aes = FALSE
)

ggsave("marker_dotplot_merged_coleoptile.pdf", p, width = 12, height =15)
```

## Source PDF page 33

```text
markers_by_group <- list(
Mesophyll_cells = c("TraesCS2D02G065200","TraesCS2A02G187200","TraesCS7D02G276300","TraesCS2A02G206200"),
Epidermal_cell = c("TraesCS3B02G167000","TraesCS5D02G440700","TraesCS5B02G145900","TraesCS1B02G292200"),
Vascular_cell = c("TraesCS4B02G237300","TraesCS7D02G358600","TraesCS6A02G218800","TraesCS7D02G126000","TraesCS1A02G212300","TraesCS2D02G431500"),
Proliferating_cell = c("TraesCS1D02G350900","TraesCS4B02G268400"),
Hydathode_cells = c("TraesCS5A02G017900")
)

genes <- unlist(markers_by_group, use.names = FALSE)
gene_symbols <- genes
symbol2cell <- unlist(lapply(names(markers_by_group), function(g) {
s <- markers_by_group[[g]]
setNames(rep(g, length(s)), s)
}))
dp <- DotPlot(merged, features = genes, assay = "RNA")
d <- dp$data
present_genes <- intersect(genes, unique(as.character(d$features.plot)))
if(length(present_genes) < length(genes)){
warning("       Seurat                :\n", paste(setdiff(genes, present_genes), collapse = ", "))
}
genes <- genes[genes %in% present_genes]
gene_symbols <- gene_symbols[genes %in% present_genes]
pos_map <- setNames(seq_along(genes), genes)
d$feat <- as.character(d$features.plot)
d$ypos <- pos_map[d$feat]
d$label <- gene_symbols[match(d$feat, genes)]
celltypes <- unique(unname(symbol2cell[genes]))
celltypes <- names(markers_by_group)
group_ranges <- lapply(celltypes, function(ct) {
genes_in_ct <- intersect(markers_by_group[[ct]], genes)
if(length(genes_in_ct) == 0) return(NULL)
idxs <- pos_map[genes_in_ct]
data.frame(celltype = ct,
minpos = min(idxs),
maxpos = max(idxs),
mid = (min(idxs) + max(idxs))/2,
stringsAsFactors = FALSE)
```

## Source PDF page 34

```text
})
group_df_unique <- do.call(rbind, group_ranges)
d$idx <- as.numeric(as.character(d$id))
if(any(is.na(d$idx))) stop("cluster id           cluster labels0,1,...        ")

p <- ggplot(d, aes(x = idx, y = ypos, size = pct.exp, color = avg.exp)) +
geom_point(shape = 16) +
scale_y_continuous(
breaks = NULL,
labels = NULL,
expand = c(0.01, 0)
) +
scale_x_continuous(
breaks = 0:14,
labels = 0:14,
expand = c(0, 0)
) +
scale_size(range = c(1, 9), guide = guide_legend(reverse = TRUE)) +
scale_color_distiller(palette = "PuRd", direction = 1) +
theme_minimal(base_size = 12) +
theme(
axis.title = element_blank(),
axis.text.x = element_text(size = 11),
axis.text.y = element_blank(),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
plot.margin = unit(c(0.5, 3, 0.5, 10), "cm"),
axis.line.y = element_line(color = "black"),
axis.line.x = element_line(color = "black"),
axis.ticks.y = element_blank()
) +
coord_cartesian(
xlim = c(-0.5, 14.5),
clip = "off"
) +
geom_segment(
data = group_df_unique,
aes(x = -8, xend = -8, y = minpos - 0.4, yend = maxpos + 0.4),
inherit.aes = FALSE,
linewidth = 0.6
) +
geom_segment(
data = group_df_unique,
aes(x = -8.05, xend = -7.95, y = minpos - 0.4, yend = minpos - 0.4),
```

## Source PDF page 35

```text
inherit.aes = FALSE,
linewidth = 0.6
) +
geom_segment(
data = group_df_unique,
aes(x = -8.05, xend = -7.95, y = maxpos + 0.4, yend = maxpos + 0.4),
inherit.aes = FALSE,
linewidth = 0.6
) +
geom_text(
data = group_df_unique,
aes(x = -8.2, y = mid, label = celltype),
inherit.aes = FALSE,
hjust = 1,
vjust = 0.5,
size = 4
) +
geom_text(
data = data.frame(y = seq_along(genes), symbol = gene_symbols, stringsAsFactors = FALSE),
aes(x = -1, y = y, label = symbol),
inherit.aes = FALSE,
hjust = 1,
vjust = 0.5,
size = 3.5
)

ggsave("marker_dotplot_merged_leaf.pdf", p, width = 10, height = max(6, length(genes)*0.35), limitsize = FALSE)
df <- merged.markers %>%
dplyr::select(cluster, p_val_adj, gene) %>%
dplyr::arrange(p_val_adj)

rownames(df) <- NULL

write.csv(df, file = "merged_markers.csv", row.names = FALSE)
```

## Source PDF page 37

```text
geneset_list <- list(
Epidermal_Cells = c(
"TraesCS2D02G473700","TraesCS2B02G497500","TraesCS2A02G474000","TraesCS7B02G208600",
"TraesCS7A02G308400","TraesCS4A02G007400","TraesCS1B02G052800","TraesCS3D02G052100",
"TraesCS7B02G220000","TraesCS1B02G436300","TraesCS3D02G272400","TraesCS3B02G266100",
"TraesCS3B02G062800","TraesCS7D02G210500","TraesCS4A02G407300","TraesCS3D02G361300",
"TraesCS3A02G368400","TraesCS3D02G112000","TraesCS2B02G583800","TraesCS7D02G350500",
"TraesCS3B02G133700","TraesCS1D02G414600","TraesCS2B02G347200","TraesCS2D02G294500",
"TraesCS5A02G238300","TraesCS4A02G021300","TraesCS3D02G338300","TraesCS6B02G253500",
"TraesCS5D02G234800","TraesCS1D02G148900","TraesCS7B02G408900","TraesCS3B02G306600",
"TraesCS3B02G063000","TraesCS3B02G209100","TraesCS6B02G241000","TraesCS5A02G093000",
"TraesCS2D02G324800","TraesCS1D02G126600","TraesCS7D02G381200","TraesCS3D02G074900",
"TraesCS6D02G141500","TraesCS1B02G093600","TraesCS7B02G159400","TraesCS1A02G341300",
"TraesCSU02G166000","TraesCS3A02G132600","TraesCS2D02G115800","TraesCS2D02G597400",
"TraesCS1D02G161500","TraesCS1B02G168500","TraesCS2D02G379500","TraesCS2D02G268400",
"TraesCS7D02G093200","TraesCS1A02G418800","TraesCS5B02G145900","TraesCS2D02G503600",
"TraesCS2D02G349000","TraesCS1A02G166100","TraesCS2D02G348900","TraesCS7B02G273400",
"TraesCS3D02G238300","TraesCS7B02G354800","TraesCS4D02G263500"
),

Proliferating_Cells_G2M = c(
"TraesCS7D02G181300","TraesCS7B02G085000","TraesCS7A02G179800","TraesCS7A02G062700",
"TraesCS7D02G057400","TraesCS4A02G056100","TraesCS2A02G521300","TraesCS2A02G202400",
"TraesCS7B02G142600","TraesCS3B02G297900","TraesCS1D02G283300","TraesCS2D02G227600",
"TraesCS5B02G199400","TraesCS7A02G239600","TraesCS2A02G412600","Tra
```

## Source PDF page 38

```text
"TraesCS4A02G101200","TraesCS5B02G037500","TraesCS6A02G161200","TraesCS7A02G335100",
"TraesCS5B02G144100","TraesCS3B02G126500","TraesCS4A02G117700","TraesCS6A02G268600",
"TraesCS4A02G206100","TraesCS2B02G250000","TraesCS6A02G103200","TraesCS7A02G346300",
"TraesCS5A02G401900","TraesCS1D02G095500","TraesCS5D02G313800","TraesCS5D02G231200",
"TraesCS1A02G230600","TraesCS5D02G206500","TraesCS3A02G277100","TraesCS3B02G246600",
"TraesCS5A02G032500","TraesCS1A02G054400","TraesCS3A02G123800","TraesCS5B02G114300",
"TraesCS2D02G154600","TraesCS4B02G248800","TraesCS7A02G008900","TraesCS2B02G422800",
"TraesCS7A02G279000","TraesCS4D02G188000","TraesCS2D02G401700","TraesCSU02G082100",
"TraesCS3B02G196600","TraesCS4B02G138000","TraesCS5D02G421200","TraesCS1D02G386100",
"TraesCS3D02G190400","TraesCS1D02G302100","TraesCS3A02G172000","TraesCS5B02G395000",
"TraesCS3A02G153700","TraesCS5D02G044700","TraesCS4A02G114600","TraesCS2B02G363300",
"TraesCS2D02G093700","TraesCS1A02G170100","TraesCS6A02G379000","TraesCS2A02G445000",
"TraesCS7B02G472800","TraesCS4B02G289700","TraesCS2A02G275700","TraesCS3B02G183400",
"TraesCS3A02G275300","TraesCS2D02G065700","TraesCS5D02G386000","TraesCS6A02G325300",
"TraesCS7A02G549500","TraesCS7D02G053600","TraesCS1A02G284300","TraesCS3B02G363200",
"TraesCS7D02G247700","TraesCS5A02G415900","TraesCS2B02G263800","TraesCS2B02G415500",
"TraesCS1B02G210600","TraesCS1B02G290400","TraesCS3A02G121700","TraesCS4B02G137400",
"TraesCS1B02G147800","TraesCS4B02G077500","TraesCS4A02G128400","TraesCS6A02G314400",
"TraesCS2D02G292800","TraesCS6D02G403700","TraesCS1D02G309300","TraesCS5B02G415800",
"TraesCS2A02G165900","TraesCS6B02G195800","TraesCS5D02G078500","TraesCS1D02G309400",
"TraesCSU02G115700","TraesCS7B02G395100"
),

Guard_Cells = c(
"TraesCS5B02G236800","TraesCS5A02G238300","TraesCS6B02G253500","Tra
```

## Source PDF page 39

```text
"TraesCS1B02G162500","TraesCS2D02G528900","TraesCS2A02G526100","TraesCS7B02G232000",
"TraesCS2A02G120700","TraesCS3B02G199100","TraesCS6B02G007700","TraesCS2D02G195800",
"TraesCS3A02G073200","TraesCS4A02G013600","TraesCS2A02G343000","TraesCS6A02G172800",
"TraesCS7A02G355700","TraesCS1D02G289400","TraesCS1A02G037400","TraesCS1A02G332300",
"TraesCS1D02G043900","TraesCS2A02G357100","TraesCS2D02G395700","TraesCS4B02G236700"
),

Mesophyll_Cells = c(
"TraesCS3A02G230000","TraesCS3D02G223300","TraesCS3B02G259300","TraesCS3D02G516800",
"TraesCS3B02G577800","TraesCS5A02G457500","TraesCS3B02G398300","TraesCS4D02G003600",
"TraesCS1D02G428200","TraesCS2B02G264100","TraesCS2A02G247300","TraesCS5A02G322500",
"TraesCS7B02G388800","TraesCS2B02G312600","TraesCS2A02G235000","TraesCS6A02G025700",
"TraesCS2B02G369600","TraesCS5B02G256400","TraesCS2B02G342200","TraesCS4D02G168800",
"TraesCS4D02G134900","TraesCS2D02G569700","TraesCS3A02G170700","TraesCS4A02G246100",
"TraesCS5A02G482800","TraesCS7B02G148000","TraesCS5D02G219100","TraesCS2D02G214700",
"TraesCS3D02G081300","TraesCS4B02G100700","TraesCS6D02G287300","TraesCS5D02G464700",
"TraesCS2D02G220800","TraesCS6B02G350500","TraesCS6D02G162200","TraesCS2B02G079500",
"TraesCS2B02G273900","TraesCS6D02G273900","TraesCS1B02G420100","TraesCS2B02G342000",
"TraesCS7D02G016800","TraesCS2D02G097500","TraesCS4D02G205400","TraesCS7B02G283000",
"TraesCS5B02G083300","TraesCS6A02G303800","TraesCS5A02G535400","TraesCS5A02G301100",
"TraesCS5B02G199200","TraesCS7A02G406800","TraesCS6B02G015000","TraesCS2A02G075000",
"TraesCS1B02G186300","TraesCS4B02G109900","TraesCS6D02G356100","TraesCS4B02G170800",
"TraesCS4D02G107400","TraesCS2B02G240000","TraesCS7B02G186500","TraesCS7B02G084000",
"TraesCS1D02G292500","TraesCS5D02G178200","TraesCS2B02G233400","TraesCS3B02G453600",
```

## Source PDF page 40

```text
"TraesCS7D02G373100","TraesCS2D02G200700","TraesCS5B02G094300","TraesCS5D02G044400",
"TraesCS2A02G297300","TraesCS7D02G209000","TraesCS6A02G290600","TraesCS1A02G273900",
"TraesCS7D02G359900","TraesCS4D02G212000","TraesCS7A02G314100","TraesCS7B02G129300",
"TraesCS7B02G247200","TraesCS2D02G253200","TraesCS2A02G114400","TraesCS7D02G000300",
"TraesCS3B02G242400","TraesCS4A02G499500","TraesCS7A02G276400","TraesCS1B02G432700",
"TraesCS3B02G344200","TraesCS2B02G272300","TraesCS6A02G159800","TraesCS7D02G276300",
"TraesCS2A02G252600"
),

Proliferating_Cells_S = c(
"TraesCS7D02G477300","TraesCS7A02G491100","TraesCS7B02G394900","TraesCS1D02G396600",
"TraesCS7D02G477400","TraesCS2D02G220400","TraesCS1A02G306500","TraesCS7D02G381300",
"TraesCS4D02G076900","TraesCS6B02G049300","TraesCS1B02G192500","TraesCS5D02G392100",
"TraesCS1A02G005700","TraesCS6B02G081800","TraesCS2B02G109900","TraesCS6B02G341900",
"TraesCS4B02G052900","TraesCS1A02G295400","TraesCS1B02G416600","TraesCS5D02G400300",
"TraesCS1D02G002400","TraesCS6B02G376900","TraesCS4D02G246300","TraesCS6B02G145200",
"TraesCS2D02G381600","TraesCS1D02G331700","TraesCS4D02G234700","TraesCS5D02G474800",
"TraesCS1B02G357600","TraesCS7B02G230000","TraesCS5B02G024300","TraesCS7D02G299600",
"TraesCS5B02G493100","TraesCS1A02G025700","TraesCS5B02G448400","TraesCS5A02G142600",
"TraesCS7A02G204800","TraesCS2D02G241800"
),

Phloem_Cells = c(
"TraesCS3A02G002000","TraesCS3B02G001500","TraesCSU02G056100","TraesCS3B02G032900",
"TraesCS5B02G182600","TraesCS4B02G206800","TraesCS5D02G161000","TraesCS3D02G256300",
"TraesCS6A02G309300","TraesCS2D02G225000","TraesCS4B02G218800","TraesCS4A02G130000",
"TraesCS1B02G266700","TraesCS2D02G423000","TraesCS5A02G519300","TraesCS6A02G274000",
```

## Source PDF page 41

```text
"TraesCS3A02G202000","TraesCS6D02G303600","TraesCS3D02G109400","TraesCS4D02G176500",
"TraesCS6B02G188000","TraesCS6B02G296800","TraesCS4B02G140500","TraesCS5A02G082200",
"TraesCS3B02G288000","TraesCS5D02G335500","TraesCS3D02G318700","TraesCS4A02G230800",
"TraesCS5A02G078100","TraesCS5D02G090400","TraesCS6D02G009300","TraesCS6B02G157700",
"TraesCS1B02G261100","TraesCS6B02G342600","TraesCS5A02G103500","TraesCS3A02G283300",
"TraesCS7B02G302900","TraesCS5D02G360700"
),

Xylem_Cells = c(
"TraesCS7D02G191600","TraesCS7A02G190600","TraesCS7B02G095500","TraesCS6B02G337300",
"TraesCS6D02G287800","TraesCS7D02G302400","TraesCS3B02G120700","TraesCS1B02G098000",
"TraesCS4A02G087000","TraesCS3A02G409200","TraesCS4A02G178900","TraesCS7A02G331100",
"TraesCS4A02G099100","TraesCS4A02G004100","TraesCS5D02G385300","TraesCS2A02G444000",
"TraesCS6D02G378300","TraesCS7D02G230700","TraesCS4B02G320300","TraesCS2A02G050800",
"TraesCS4A02G044500","TraesCS2B02G334600","TraesCS5D02G404300","TraesCS7D02G094300",
"TraesCS4D02G133800","TraesCS6A02G084500","TraesCS6A02G257900","TraesCS2A02G206500",
"TraesCS4A02G067000","TraesCS7B02G148100","TraesCS3A02G246600","TraesCS4D02G265500",
"TraesCS6A02G028700","TraesCS4D02G217700","TraesCS5A02G138700","TraesCS2D02G018000",
"TraesCS2A02G198000","TraesCS2D02G467000","TraesCS6D02G220600","TraesCS5D02G255700",
"TraesCS3D02G289600","TraesCS4A02G314800","TraesCS7A02G276300","TraesCS2A02G292600",
"TraesCS2A02G541800","TraesCS3B02G539600","TraesCS1A02G254800","TraesCS7A02G008800"
),

Parenchyma_Cells = c(
"TraesCS3B02G441300","TraesCS3D02G402200","TraesCS3B02G441100","TraesCS3D02G402300",
"TraesCS2A02G268200","TraesCS2B02G530500","TraesCS5B02G535900","TraesCS3A02G077700",
"TraesCS7A02G251300","TraesCS5D02G552900","TraesCS7D02G196700","Tra
```

## Source PDF page 42

```text
"TraesCS4A02G380200","TraesCS5B02G319300","TraesCS1A02G051800","TraesCS4D02G095300",
"TraesCS5B02G155900","TraesCS5B02G011300","TraesCS2B02G399800","TraesCS3A02G324000",
"TraesCSU02G025500","TraesCS3D02G357700","TraesCS6B02G028400","TraesCS5A02G262700",
"TraesCS5B02G193100","TraesCS6D02G242100","TraesCS6A02G340100","TraesCS5D02G200800",
"TraesCS4B02G090800","TraesCS1A02G311100","TraesCS5A02G225700","TraesCS5A02G205200",
"TraesCS7D02G412400","TraesCS7D02G491400","TraesCS1B02G248700","TraesCS4A02G133100",
"TraesCS2A02G484700","TraesCS7D02G067300"
),

Hydathode_Cells = c(
"TraesCS3B02G293200","TraesCS5B02G181500","TraesCS3A02G483000","TraesCS3A02G260100",
"TraesCS7D02G161200","TraesCS5A02G004800","TraesCS3A02G260200","TraesCS7A02G558500",
"TraesCS4A02G296300","TraesCS5A02G183300"
),

Vascular_Cells = c(
"TraesCS6A02G404500","TraesCS5A02G337400","TraesCS5D02G342000","TraesCS2A02G511000",
"TraesCS5B02G336300","TraesCS4B02G237300","TraesCS7D02G355900","TraesCS4A02G061700",
"TraesCS7A02G366900","TraesCS6D02G314300"
),

Proliferating_Cells = c(
"TraesCS3B02G131300","TraesCS3B02G130500","TraesCS3B02G131000","TraesCS4D02G212300",
"TraesCS2B02G047400","TraesCS4D02G145500","TraesCS3D02G045500","TraesCS4D02G213300",
"TraesCS4D02G213100","TraesCS6D02G322300"
),

Companion_Cells = c(
"TraesCS3B02G441300","TraesCS3D02G402200","TraesCS3B02G441100","TraesCS3D02G402300",
"TraesCS3A02G406800","TraesCS3D02G357900","TraesCS5D02G423800","TraesCS2A02G488200",
"TraesCS4B02G348900","TraesCS4B02G292700","TraesCSU02G198800","Trae
```

## Source PDF page 43

```text
"TraesCS3B02G225200","TraesCS1D02G426700","TraesCS2A02G386500","TraesCS4D02G176500",
"TraesCS7A02G403500","TraesCS3B02G153400","TraesCS6A02G155400","TraesCS3D02G125300",
"TraesCS5B02G222800","TraesCS1B02G351600","TraesCS5D02G498500","TraesCS5A02G092000"
)
)


merged_celltype <- c(
"Epidermal_Cells", "Proliferating_Cells_G2M", "Guard_Cells",
"Mesophyll_Cells", "Proliferating_Cells_S", "Phloem_Cells",
"Xylem_Cells", "Parenchyma_Cells", "Hydathode_Cells",
"Vascular_Cells", "Proliferating_Cells", "Companion_Cells"
)

for(celltype in merged_celltype) {

genes <- geneset_list[[celltype]]
genes = unique(genes)
geneSets <- GeneSet(genes, setName="Scores")
assay <- GetAssay(merged, "RNA")
exprMatrix <- GetAssayData(assay, "counts")
cells_AUC <- AUCell_run(exprMatrix, geneSets)
merged$score <- getAUC(cells_AUC)["Scores",]

p <- FeaturePlot(merged, "score") +
ggtitle(celltype) +
theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(celltype, "_merged_featureplot_ALL.pdf"), p, width =8, height = 8)

print(paste(celltype, " "))
}
```

## Source PDF page 44

```text
geneset_list <- list(
Epidermal_Cells = c(
"TraesCS1A02G418800","TraesCS5B02G145900","TraesCS2D02G503600",
"TraesCS2D02G349000","TraesCS1A02G166100","TraesCS2D02G348900","TraesCS7B02G273400",
"TraesCS3D02G238300","TraesCS7B02G354800","TraesCS4D02G263500"
),

Mesophyll_Cells = c(
"TraesCS4A02G499500","TraesCS7A02G276400","TraesCS5A02G322500","TraesCS1B02G432700",
"TraesCS2D02G253200","TraesCS3B02G344200","TraesCS2B02G272300","TraesCS6A02G159800",
"TraesCS7D02G276300","TraesCS2A02G252600"
),

Hydathode_Cells = c(
"TraesCS3B02G293200","TraesCS5B02G181500","TraesCS3A02G483000","TraesCS3A02G260100",
"TraesCS7D02G161200","TraesCS5A02G004800","TraesCS3A02G260200","TraesCS7A02G558500",
"TraesCS4A02G296300","TraesCS5A02G183300"
),

Vascular_Cells = c(
"TraesCS6A02G404500","TraesCS5A02G337400","TraesCS5D02G342000","TraesCS2A02G511000",
"TraesCS5B02G336300","TraesCS4B02G237300","TraesCS7D02G355900","TraesCS4A02G061700",
"TraesCS7A02G366900","TraesCS6D02G314300"
),

Proliferating_Cells = c(
"TraesCS3B02G131300","TraesCS3B02G130500","TraesCS3B02G131000","TraesCS4D02G212300",
"TraesCS2B02G047400","TraesCS4D02G145500","TraesCS3D02G045500","TraesCS4D02G213300",
"TraesCS4D02G213100","TraesCS6D02G322300"
)
)


merged_celltype <- c(
```

## Source PDF page 45

```text
"Epidermal_Cells", "Mesophyll_Cells", "Hydathode_Cells",
"Vascular_Cells", "Proliferating_Cells"
)

for(celltype in merged_celltype) {

genes <- geneset_list[[celltype]]
genes = unique(genes)
geneSets <- GeneSet(genes, setName="Scores")
assay <- GetAssay(merged, "RNA")
exprMatrix <- GetAssayData(assay, "counts")
cells_AUC <- AUCell_run(exprMatrix, geneSets)
merged$score <- getAUC(cells_AUC)["Scores",]

p <- FeaturePlot(merged, "score") +
ggtitle(celltype) +
theme(plot.title = element_text(hjust = 0.5))

ggsave(paste0(celltype, "_merged_featureplot_top10_degs.pdf"), p, width = 8, height = 8)

print(paste(celltype, " "))
}
```

## Source PDF page 46

```text
all_clusters <- sort(unique(merged.markers$cluster))

for (i in all_clusters) {
message("Processing Cluster: ", i)
cluster_data <- merged.markers %>%
filter(cluster == i & p_val_adj < 0.05 & avg_log2FC > 0)
deg_genes <- cluster_data$gene
deg_genes <- deg_genes[deg_genes != "" & !is.na(deg_genes)]
if (length(deg_genes) == 0) {
message(" Skipping Cluster ", i, ": No significant marker genes found.")
next
}
ego_BP <- tryCatch({
enrichGO(
gene        = deg_genes,
OrgDb       = org.Ta.eg.db,
keyType     = "GENEID",
ont         = "BP",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05
)
}, error = function(e) NULL)
ego_CC <- tryCatch({
enrichGO(
gene        = deg_genes,
OrgDb       = org.Ta.eg.db,
keyType     = "GENEID",
ont         = "CC",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
qvalueCutoff = 0.05
)
}, error = function(e) NULL)
ego_MF <- tryCatch({
enrichGO(
gene        = deg_genes,
OrgDb       = org.Ta.eg.db,
keyType     = "GENEID",
ont         = "MF",
pAdjustMethod = "BH",
pvalueCutoff = 0.05,
```

## Source PDF page 47

```text
qvalueCutoff = 0.05
)
}, error = function(e) NULL)

df_BP <- if (!is.null(ego_BP)) as.data.frame(ego_BP) %>% mutate(type= "biological process") else data.frame()
df_CC <- if (!is.null(ego_CC)) as.data.frame(ego_CC) %>% mutate(type= "cellular component") else data.frame()
df_MF <- if (!is.null(ego_MF)) as.data.frame(ego_MF) %>% mutate(type= "molecular function") else data.frame()

final_df <- bind_rows(df_BP, df_CC, df_MF)

if (nrow(final_df) > 0) {
file_name <- paste0("merged_cluster_", i, "_GO.csv")

write.csv(final_df, file = file_name, row.names = FALSE)
message(" Saved: ", file_name)
} else {
message(" No GO terms enriched for Cluster ", i)
}
}
```

## Source PDF page 48

```text
library(Seurat)
library(ggplot2)
library(ggrepel)
library(grid)
library(dplyr)

#=====================================================
# 1.    cluster
#=====================================================

merged <- RenameIdents(
merged,
"0" = "Defense-detoxification cells",
"1" = "Mesophyll cells I",
"2" = "Wounding-responsive cells",
"3" = "Defense-thickening cells I",
"4" = "Protein-synthesis cells",
"5" = "JA-mediated defense cells",
"6" = "Mesophyll cells II",
"7" = "Defense-secondary metabolism cells",
"8" = "Epidermis cells I",
"9" = "Defense-thickening cells II",
"10" = "Trichome cells",
"11" = "Proliferating-S cells",
"12" = "Epidermis cells II",
"13" = "Phloem cells",
"14" = "Guard cells"
)

merged$celltype <- Idents(merged)
merged@meta.data$celltype <- Idents(merged)

#=====================================================
# 2.
#=====================================================

cluster_levels <- as.character(0:14)

cluster_names <- c(
"Defense-detoxification cells",
"Mesophyll cells I",
"Wounding-responsive cells",
"Defense-thickening cells I",
```

## Source PDF page 49

```text
"Protein-synthesis cells",
"JA-mediated defense cells",
"Mesophyll cells II",
"Sclerenchyma cells",
"Epidermis cells I",
"Defense-thickening cells II",
"Trichome cells",
"Proliferating-S cells",
"Epidermis cells II",
"Phloem cells",
"Guard cells"
)

names(cluster_names) <- cluster_levels

mycols_named <- setNames(mycols, cluster_levels)

#=====================================================
# 3.
#=====================================================

draw_key_circled_num <- function(data, params, size){

idx <- names(mycols_named)[match(data$colour,
unname(mycols_named))]

grobTree(
circleGrob(
r = 0.42,
gp = gpar(fill = data$colour, col = NA)
),
textGrob(
idx,
x = 0.5,
y = 0.5,
gp = gpar(
col = "white",
fontsize = 8,
fontface = "bold"
)
)
)
}

#=====================================================
# 4.   UMAP
```

## Source PDF page 50

```text
#=====================================================

umap_df <- as.data.frame(Embeddings(merged, "umap"))
colnames(umap_df) <- c("UMAP_1","UMAP_2")

umap_df$cluster <- as.character(merged$seurat_clusters)

#=====================================================
# 5.   cluster
#=====================================================

center <- umap_df %>%
group_by(cluster) %>%
summarise(
cx = mean(UMAP_1),
cy = mean(UMAP_2),
.groups = "drop"
)

#=====================================================
# 6.
#=====================================================

celltypepos <- umap_df %>%
left_join(center, by = "cluster") %>%
mutate(
dist = (UMAP_1 - cx)^2 + (UMAP_2 - cy)^2
) %>%
group_by(cluster) %>%
slice_min(dist, n = 1, with_ties = FALSE) %>%
ungroup()

#=====================================================
# 7.   12 cluster
#=====================================================

celltypepos[celltypepos$cluster=="12","UMAP_1"] <-
celltypepos[celltypepos$cluster=="12","UMAP_1"] - 0.8

celltypepos[celltypepos$cluster=="12","UMAP_2"] <-
celltypepos[celltypepos$cluster=="12","UMAP_2"] + 0.25

#=====================================================
# 8.
#=====================================================

```

## Source PDF page 51

```text
p <- DimPlot(
merged,
group.by = "seurat_clusters",
reduction = "umap",
label = FALSE,
pt.size = 0.3,
cols = mycols_named
) +

scale_color_manual(
breaks = cluster_levels,
values = mycols_named,
labels = cluster_names,
guide = guide_legend(
override.aes = list(size = 5)
)
) +

geom_label_repel(

data = celltypepos,

aes(
x = UMAP_1,
y = UMAP_2,
label = cluster,
colour = cluster
),

inherit.aes = FALSE,

fontface = "bold",

size = 4,

fill = "white",

label.size = 0.25,

box.padding = 0.35,

point.padding = 0.2,

force = 0.5,

max.overlaps = Inf,
```

## Source PDF page 52

```text

segment.color = NA,  #

show.legend = FALSE,

seed = 54
) +

labs(title = NULL) +

theme_dr(
xlength = 0.2,
ylength = 0.2,
arrow = arrow(
length = unit(0.2, "inches"),
type = "closed"
)
) +

theme(
panel.grid = element_blank(),
axis.title = element_text(
face = 2,
hjust = 0.03
),
legend.title = element_blank(),
legend.text = element_text(
size = 10,
margin = margin(l = 5)
),
legend.key = element_blank(),
legend.background = element_blank(),
legend.key.height = unit(0.8, "cm"),
legend.key.width = unit(1.1, "cm")
)

#=====================================================
# 9.
#=====================================================

p$layers[[1]]$geom$draw_key <- draw_key_circled_num

#=====================================================
# 10.
#=====================================================

```
