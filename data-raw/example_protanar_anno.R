## code to prepare `example_protanar_anno` dataset goes here

example_protanar_anno <- data.frame(
  label = c(
    "g1_1_1", "g1_1_2", "g1_2_1", "g1_2_2", "g1_3_1", "g1_3_2",
    "g2_1_1", "g2_1_2", "g2_2_1", "g2_2_2", "g2_3_1", "g2_3_2"
  ),
  group = c(
    "g1", "g1", "g1", "g1", "g1", "g1",
    "g2", "g2", "g2", "g2", "g2", "g2"
  ),
  biol_repl = c(
    1, 1, 2, 2, 3, 3,
    1, 1, 2, 2, 3, 3
  ),
  tech_repl = c(
    1, 2, 1, 2, 1, 2,
    1, 2, 1, 2, 1, 2
  )
)

usethis::use_data(example_protanar_anno, overwrite = TRUE)
