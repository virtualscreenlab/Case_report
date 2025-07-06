# --- Загрузка библиотек ---
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

# --- Директории ---
output_dir <- "C:/PPI/PatientK/Result"
readme_dir <- file.path(output_dir, "READMEs")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(readme_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# --- Загрузка и коррекция данных ---
data <- read.delim(
  "C:/PPI/PatientK/input/Genes_K.csv",
  sep = ";",
  header = TRUE,
  col.names = c("pvalue", "logFC", "gene"),
  stringsAsFactors = FALSE
)
corrected <- checkGeneSymbols(unique(data$gene), species = "human")
data$gene <- corrected$Suggested.Symbol[match(data$gene, corrected$x)]
data <- data[!is.na(data$gene), ]

# --- STRING ---
string_db <- STRINGdb$new(
  species = 9606,
  version = "12.0",
  score_threshold = 400,
  network_type = "full"
)
mapped_data <- string_db$map(data, "gene", removeUnmappedRows = TRUE)
ppi_network <- string_db$get_subnetwork(mapped_data$STRING_id)
ppi_igraph <- ppi_network

# --- Глобальные метрики ---
density <- edge_density(ppi_igraph)
assort <- assortativity_degree(ppi_igraph, directed = FALSE)
modularity_louvain <- modularity(cluster_louvain(ppi_igraph))
diam <- diameter(ppi_igraph)
avg_path <- mean_distance(ppi_igraph, directed = FALSE)
deg <- degree(ppi_igraph, mode = "all")
num_remove <- round(0.05 * vcount(ppi_igraph))
hubs <- names(sort(deg, decreasing = TRUE)[1:num_remove])
ppi_frag <- delete_vertices(ppi_igraph, hubs)
frag_largest_comp <- max(components(ppi_frag)$csize)
frag_ratio <- frag_largest_comp / vcount(ppi_igraph)

# --- Единичные метрики в один Excel ---
global_metrics <- tibble(
  density = density,
  assortativity = assort,
  modularity_louvain = modularity_louvain,
  diameter = diam,
  avg_shortest_path = avg_path,
  hubs_removed = num_remove,
  largest_comp_after_removal = frag_largest_comp,
  fraction_remaining = frag_ratio
)
write.xlsx(global_metrics, file = file.path(output_dir, "Global_Metrics.xlsx"))
writeLines(
  "Файл содержит глобальные метрики PPI-сети:
- density: плотность сети
- assortativity: коэффициент ассоциативности степеней
- modularity_louvain: модульность по Louvain
- diameter: диаметр сети
- avg_shortest_path: средняя длина кратчайших путей
- hubs_removed: число удалённых хабов (5% от всех)
- largest_comp_after_removal: размер наибольшей компоненты после удаления хабов
- fraction_remaining: доля вершин в наибольшей компоненте после удаления хабов",
  con = file.path(readme_dir, "Global_Metrics_README.txt"))

# --- Топ-20 хабов ---
top20 <- sort(deg, decreasing = TRUE)[1:20]
top20_df <- tibble(
  STRING_id = names(top20),
  degree = as.integer(top20)
) %>%
  left_join(mapped_data, by = "STRING_id") %>%
  select(STRING_id, gene, degree, pvalue, logFC)
write.xlsx(top20_df, file = file.path(output_dir, "Top20_Hubs.xlsx"))
writeLines(
  "Файл содержит топ-20 белков-хабов по числу взаимодействий (degree).
Degree рассчитывается функцией degree() из igraph.",
  con = file.path(readme_dir, "Top20_Hubs_README.txt"))

# --- Кластеры Louvain ---
communities_louvain <- cluster_louvain(ppi_igraph)
modules_df <- tibble(
  STRING_id = names(membership(communities_louvain)),
  gene = mapped_data$gene[match(names(membership(communities_louvain)), mapped_data$STRING_id)],
  module = as.integer(membership(communities_louvain))
)
write.xlsx(modules_df, file = file.path(output_dir, "Functional_Modules_Louvain.xlsx"))
writeLines(
  "Файл содержит принадлежность белков к функциональным модулям (кластеризация Louvain, igraph::cluster_louvain).
module — номер модуля.",
  con = file.path(readme_dir, "Functional_Modules_Louvain_README.txt"))

# --- Кластеры FastGreedy (альтернативная кластеризация) ---
if (is.connected(ppi_igraph) && !is.directed(ppi_igraph)) {
  communities_fg <- cluster_fast_greedy(ppi_igraph)
  modules_fg_df <- tibble(
    STRING_id = names(membership(communities_fg)),
    gene = mapped_data$gene[match(names(membership(communities_fg)), mapped_data$STRING_id)],
    module = as.integer(membership(communities_fg))
  )
  write.xlsx(modules_fg_df, file = file.path(output_dir, "Functional_Modules_FastGreedy.xlsx"))
  writeLines(
    "Файл содержит принадлежность белков к модулям по FastGreedy (igraph::cluster_fast_greedy).
module — номер модуля.",
    con = file.path(readme_dir, "Functional_Modules_FastGreedy_README.txt"))
}

# --- Ключевые регуляторы (betweenness) ---
betw <- betweenness(ppi_igraph, normalized = TRUE)
top20_betw <- sort(betw, decreasing = TRUE)[1:20]
top20_betw_df <- tibble(
  STRING_id = names(top20_betw),
  gene = mapped_data$gene[match(names(top20_betw), mapped_data$STRING_id)],
  betweenness = as.numeric(top20_betw)
)
write.xlsx(top20_betw_df, file = file.path(output_dir, "Top20_Betweenness_Regulators.xlsx"))
writeLines(
  "Файл содержит топ-20 белков по посреднической центральности (betweenness, igraph::betweenness).
Высокое значение указывает на роль белка как 'посредника' в передаче сигналов.",
  con = file.path(readme_dir, "Top20_Betweenness_Regulators_README.txt"))

# --- Эксцентриситет (eccentricity) ---
ecc <- eccentricity(ppi_igraph)
ecc_df <- tibble(
  STRING_id = names(ecc),
  gene = mapped_data$gene[match(names(ecc), mapped_data$STRING_id)],
  eccentricity = as.integer(ecc)
)
write.xlsx(ecc_df, file = file.path(output_dir, "Eccentricity.xlsx"))
writeLines(
  "Эксцентриситет — максимальное расстояние от узла до любого другого (igraph::eccentricity).
Характеризует 'отдалённость' белка в сети.",
  con = file.path(readme_dir, "Eccentricity_README.txt"))

# --- Локальный коэффициент кластеризации ---
clust <- transitivity(ppi_igraph, type = "local", isolates = "zero")
clust_df <- tibble(
  STRING_id = names(clust),
  gene = mapped_data$gene[match(names(clust), mapped_data$STRING_id)],
  clustering = clust
)
write.xlsx(clust_df, file = file.path(output_dir, "Node_Clustering_Coefficient.xlsx"))
writeLines(
  "Локальный коэффициент кластеризации (igraph::transitivity, type='local').
Показывает плотность связей вокруг каждого белка.",
  con = file.path(readme_dir, "Node_Clustering_Coefficient_README.txt"))

# --- Межмодульные связи ---
edge_df <- igraph::as_data_frame(ppi_igraph, what = "edges")
edge_df$from_comm <- membership(communities_louvain)[edge_df$from]
edge_df$to_comm <- membership(communities_louvain)[edge_df$to]
intermodular_edges <- sum(edge_df$from_comm != edge_df$to_comm)
intermodular_df <- tibble(intermodular_edges = intermodular_edges)
write.xlsx(intermodular_df, file = file.path(output_dir, "Intermodular_Edges.xlsx"))
writeLines(
  "Число межмодульных рёбер — количество рёбер между разными модулями (по Louvain).
Рассчитано через сравнение принадлежности концов рёбер.",
  con = file.path(readme_dir, "Intermodular_Edges_README.txt"))

# --- Распределение степеней ---
deg_dist <- as.data.frame(table(deg))
colnames(deg_dist) <- c("degree", "frequency")
write.xlsx(deg_dist, file = file.path(output_dir, "Degree_Distribution.xlsx"))
writeLines(
  "Файл содержит распределение степеней (degree distribution) — сколько белков имеют каждое значение degree.",
  con = file.path(readme_dir, "Degree_Distribution_README.txt"))

# --- Мультицентральность ---
closeness <- closeness(ppi_igraph, normalized = TRUE)
multi_centrality <- scale(deg) + scale(betw) + scale(closeness)
multi_centrality_df <- tibble(
  STRING_id = V(ppi_igraph)$name,
  gene = mapped_data$gene[match(V(ppi_igraph)$name, mapped_data$STRING_id)],
  multi_centrality = as.numeric(multi_centrality)
)
write.xlsx(multi_centrality_df, file = file.path(output_dir, "MultiCentrality.xlsx"))
writeLines(
  "Интегральная мультицентральность: сумма нормированных degree, betweenness и closeness для поиска ключевых белков.",
  con = file.path(readme_dir, "MultiCentrality_README.txt"))

# --- Shortest Paths Distribution (гистограмма) ---
sp_matrix <- distances(ppi_igraph)
sp_dist <- as.vector(sp_matrix[upper.tri(sp_matrix)])
sp_dist <- sp_dist[is.finite(sp_dist)]
sp_hist <- as.data.frame(table(sp_dist))
colnames(sp_hist) <- c("Shortest_Path_Length", "Frequency")
write.xlsx(sp_hist, file = file.path(output_dir, "Shortest_Paths_Distribution.xlsx"))
writeLines(
  "Гистограмма длин кратчайших путей между всеми парами белков (igraph::distances).
Shortest_Path_Length — длина пути, Frequency — число пар с такой длиной.",
  con = file.path(readme_dir, "Shortest_Paths_Distribution_README.txt"))

# --- Кластеризация: визуализация модулей Louvain и FastGreedy ---
# Louvain
png(file.path(output_dir, "PPI_Louvain_Clusters.png"), width=1800, height=1800, res=200)
ggraph(ppi_igraph, layout = "fr") +
  geom_edge_link(color = "grey80", alpha = 0.5, width = 0.3) +
  geom_node_point(aes(color = as.factor(membership(communities_louvain))), size=2) +
  theme_void() + ggtitle("Louvain Clusters")
dev.off()
writeLines(
  "PNG-файл: визуализация кластеров по Louvain (цвет = модуль)",
  con = file.path(readme_dir, "PPI_Louvain_Clusters_README.txt"))

# FastGreedy (если применимо)
if (exists("communities_fg")) {
  png(file.path(output_dir, "PPI_FastGreedy_Clusters.png"), width=1800, height=1800, res=200)
  ggraph(ppi_igraph, layout = "fr") +
    geom_edge_link(color = "grey80", alpha = 0.5, width = 0.3) +
    geom_node_point(aes(color = as.factor(membership(communities_fg))), size=2) +
    theme_void() + ggtitle("FastGreedy Clusters")
  dev.off()
  writeLines(
    "PNG-файл: визуализация кластеров по FastGreedy (цвет = модуль)",
    con = file.path(readme_dir, "PPI_FastGreedy_Clusters_README.txt"))
}

# --- Генерация ссылок на интерактом и кластеры STRING ---
organism_id <- 9606 # Homo sapiens
string_interactome_url <- paste0("https://string-db.org/cgi/network?species=", organism_id)
string_ids <- unique(mapped_data$STRING_id)
string_ids <- string_ids[1:min(300, length(string_ids))]
string_link <- string_db$get_link(string_ids)
string_clusters_url <- paste0(string_link, "&clustering_method=mcl")

links_text <- paste(
  "Ссылки на интерактивные ресурсы STRING для вашего анализа:",
  "",
  "1. Целый интерактом (все PPI для организма):",
  string_interactome_url,
  "",
  "2. Интерактивная сеть для вашего поднабора белков:",
  string_link,
  "",
  "3. Кластеризация (MCL) для вашего поднабора белков (можно выбрать и KMEANS в интерфейсе):",
  string_clusters_url,
  "",
  "Пояснения:",
  "- Интерактом — это полный граф белок-белковых взаимодействий для организма (например, человека).",
  "- Ваша интерактивная сеть — это подграф, построенный по вашему списку белков, с учётом порога доверия.",
  "- Кластеризация (MCL/KMEANS) выделяет функциональные модули, которые можно визуализировать и скачать прямо в STRING.",
  "- Подробнее о кластеризации: https://string-db.org/help/interactive_network/",
  sep = "\n"
)
writeLines(links_text, file.path(output_dir, "STRING_Links_and_Clusters.txt"))

cat("✅ Все файлы сохранены в папке", output_dir, "\nREADME-файлы — в папке", readme_dir, "\n")


