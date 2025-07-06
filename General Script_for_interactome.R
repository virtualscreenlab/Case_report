# --- Library Loading ---
if (!requireNamespace("HGNChelper", quietly = TRUE)) install.packages("HGNChelper")
if (!requireNamespace("STRINGdb", quietly = TRUE)) install.packages("STRINGdb")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")

library(HGNChelper)
library(STRINGdb)
library(igraph)
library(openxlsx)
library(dplyr)
library(tibble)
library(ggraph)
library(ggplot2)


# Set working directory
output_dir <- "path/to/your/output/folder"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# Load data
data <- read.delim(
  "path/to/your/input/folder/Genes_K.csv",
  sep = ";",
  header = TRUE,
  col.names = c("pvalue", "logFC", "gene"),
  stringsAsFactors = FALSE
) %>%
  mutate(gene = trimws(gene))

# Gene symbol validation and correction
corrected <- checkGeneSymbols(unique(data$gene), species = "human")
data <- data %>%
  mutate(
    original_gene = gene,
    gene = corrected$Suggested.Symbol[match(gene, corrected$x)]
  ) %>%
  filter(!is.na(gene)) %>%
  select(-original_gene)

# Connect to STRING database
string_db <- STRINGdb$new(
  species = 9606,
  version = "12.0",
  score_threshold = 400,
  network_type = "full"
)

# Map data to STRING identifiers
mapped_data <- string_db$map(data, "gene", removeUnmappedRows = TRUE)

# Retrieve protein-protein interaction network
ppi_network <- string_db$get_subnetwork(mapped_data$STRING_id)

# Convert to tidygraph format
ppi_tbl <- as_tbl_graph(ppi_network) %>%
  left_join(mapped_data, by = c("name" = "STRING_id")) %>%
  mutate(community = as.factor(group_louvain()))

# --- Visualization in PNG and PDF formats ---
p <- ggraph(ppi_tbl, layout = "fr") +
  geom_edge_link(color = "grey80", alpha = 0.5, width = 0.3) +
  geom_node_point(aes(
    color = community,
    size = ifelse(is.na(logFC), 1, abs(logFC))
  )) +
  geom_node_text(
    aes(label = ifelse(abs(logFC) > 1.5, gene, "")),
    repel = TRUE,
    size = 2.5,
    color = "black",
    max.overlaps = 100
  ) +
  scale_color_viridis_d() +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  ) +
  ggtitle("STRING-like Protein Interaction Network")

ggsave("PPI_STRING_Style_HQ.png", plot = p, width = 16, height = 16, dpi = 600, bg = "white")
ggsave("PPI_STRING_Style_HQ.pdf", plot = p, width = 16, height = 16)

# --- Generate interactive STRING link (max 300 nodes) ---
string_ids <- unique(mapped_data$STRING_id)
string_ids <- string_ids[1:min(300, length(string_ids))]  # Limit to 300 nodes
string_link <- string_db$get_link(string_ids)
cat("\n🔗 Interactive STRING network link:\n", string_link, "\n")
browseURL(string_link)

# --- Top 30 proteins by number of interactions ---
top30 <- degree(ppi_network, mode = "all") %>%
  sort(decreasing = TRUE) %>%
  head(30) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("STRING_id") %>%
  rename(connections = ".")

top30 <- top30 %>%
  left_join(mapped_data, by = "STRING_id") %>%
  select(STRING_id, connections, gene, pvalue)

# --- Excel table with top 30 proteins ---
wb <- createWorkbook()
addWorksheet(wb, "Top30_Proteins")
writeData(wb, "Top30_Proteins", top30)
setColWidths(wb, "Top30_Proteins", cols = 1:4, widths = "auto")
saveWorkbook(wb, "Top30_Proteins.xlsx", overwrite = TRUE)

# --- Node degree distribution analysis ---

# Calculate degrees for all nodes
node_degrees <- degree(ppi_network, mode = "all")

# Count distribution: how many nodes have each possible degree
degree_table <- as.data.frame(table(node_degrees))
colnames(degree_table) <- c("n_neighbors", "n_nodes")
degree_table$n_neighbors <- as.numeric(as.character(degree_table$n_neighbors))
degree_table$n_nodes <- as.numeric(degree_table$n_nodes)

# Remove zero degrees (if any)
degree_table <- degree_table[degree_table$n_neighbors > 0, ]

# --- Power-law approximation k = C * n^a ---
# Logarithmic transformation for linear regression
degree_table$log_n_neighbors <- log10(degree_table$n_neighbors)
degree_table$log_n_nodes <- log10(degree_table$n_nodes)

fit <- lm(log_n_nodes ~ log_n_neighbors, data = degree_table)
a <- coef(fit)[2] # exponent
C <- 10^coef(fit)[1] # coefficient
r2 <- summary(fit)$r.squared

# Add approximated values
degree_table$predicted_n_nodes <- C * degree_table$n_neighbors^a

# --- Save data to Excel ---
addWorksheet(wb, "Degree_Distribution")
writeData(wb, "Degree_Distribution", degree_table)
writeData(wb, "Degree_Distribution", data.frame(
  Equation = sprintf("k = %.2f * n^%.2f", C, a),
  R2 = r2
), startRow = nrow(degree_table) + 3)
saveWorkbook(wb, "Top30_Proteins.xlsx", overwrite = TRUE)

# --- Degree distribution visualization ---
library(scales)
p_deg <- ggplot(degree_table, aes(x = n_neighbors, y = n_nodes)) +
  geom_point(size = 2, color = "#21908CFF") +
  geom_line(aes(y = predicted_n_nodes), color = "red", linetype = "dashed") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  labs(
    title = "Degree Distribution (log-log)",
    x = "Number of neighbors (n)",
    y = "Number of nodes (k)",
    subtitle = sprintf("Approximation: k = %.2f * n^%.2f, R² = %.3f", C, a, r2)
  ) +
  theme_minimal(base_size = 14)

ggsave("Degree_Distribution_loglog.png", plot = p_deg, width = 8, height = 6, dpi = 300)
ggsave("Degree_Distribution_loglog.pdf", plot = p_deg, width = 8, height = 6)

cat("\n✅ Node degree distribution analysis completed!\n",
    "• Plot: Degree_Distribution_loglog.png\n",
    "• PDF: Degree_Distribution_loglog.pdf\n",
    "• Excel (2nd sheet): Top30_Proteins.xlsx\n",
    sprintf("• Approximation: k = %.2f * n^%.2f, R² = %.3f\n", C, a, r2)
)


# --- Summary ---
cat("\n✅ Files saved:\n",
    "• PNG plot: PPI_STRING_Style_HQ.png\n",
    "• PDF plot: PPI_STRING_Style_HQ.pdf\n",
    "• Excel table: Top30_Proteins.xlsx\n",
    "• STRING link: ", string_link, "\n")
