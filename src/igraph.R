library(igraph)

edges <- core_ped %>%
  select(child = id, father = pid, mother = mid) %>%
  pivot_longer(cols = c(father, mother), names_to = "parent_type", values_to = "parent_id") %>%
  filter(!is.na(parent_id)) %>%
  select(parent_id, child)

g <- graph_from_data_frame(edges, directed = TRUE, vertices = core_ped$id)

png("igraph.png", width=800, height=600)  # You can adjust width and height as needed
plot(g, 
     vertex.label = V(g)$name,
     vertex.color = ifelse(V(g)$Y_ALS == 1, "red", "lightgray"),
     vertex.size = 30,
     edge.arrow.size = 0.5,
     main = "Pedigree Graph: All Individuals")
dev.off()

png("subgraph.png", width=800, height=600)  # You can adjust width and height as needed
plot(subg,
     vertex.label = V(subg)$name,
     vertex.color = "red",
     vertex.size = 30,
     main = "Affected Individuals Subgraph")
dev.off()

paths <- sapply(1:ncol(affected_pairs), function(i) {
  pair <- affected_pairs[, i]
  sp <- shortest_paths(g, from = pair[1], to = pair[2])$vpath[[1]]
  length(sp) - 1  # number of edges in the path
})


png("kinship_graph.png", width=800, height=600)  # You can adjust width and height as needed
plot(kinship_graph, 
     vertex.label = V(g)$name,
     vertex.color = ifelse(V(g)$Y_ALS == 1, "red", "lightgray"),
     vertex.size = 30,
     edge.arrow.size = 0.5,
     main = "Pedigree Graph: All Individuals")
dev.off()
