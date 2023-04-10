## code to prepare `example_protanar_1` dataset goes here

set.seed(42)

example_protanar_1 <- data.frame(
	id = c(
		'p1', 'p1', 'p1', 'p1', 'p1', 'p1', 
		'p1', 'p1', 'p1', 'p1', 'p1', 'p1', 
		'p2', 'p2', 'p2', 'p2', 'p2', 'p2', 
		'p2', 'p2', 'p2', 'p2', 'p2', 'p2', 
		'p3', 'p3', 'p3', 'p3', 'p3', 'p3',
		'p3', 'p3', 'p3', 'p3', 'p3', 'p3' 
		),
	label = c(
		'g1_1_1', 'g1_1_2', 'g1_2_1', 'g1_2_2', 'g1_3_1', 'g1_3_2',
		'g2_1_1', 'g2_1_2', 'g2_2_1', 'g2_2_2', 'g2_3_1', 'g2_3_2',
		'g1_1_1', 'g1_1_2', 'g1_2_1', 'g1_2_2', 'g1_3_1', 'g1_3_2',
		'g2_1_1', 'g2_1_2', 'g2_2_1', 'g2_2_2', 'g2_3_1', 'g2_3_2',
		'g1_1_1', 'g1_1_2', 'g1_2_1', 'g1_2_2', 'g1_3_1', 'g1_3_2',
		'g2_1_1', 'g2_1_2', 'g2_2_1', 'g2_2_2', 'g2_3_1', 'g2_3_2'
		),
	LFQ = c(
		stats::rnorm(6, 10, 1),
		stats::rnorm(6, 10, 1),
		stats::rnorm(6, 15, 1),
		stats::rnorm(6, 10, 1),
		stats::rnorm(6, 10, 1),
		stats::rnorm(6, 15, 1)
		)
)

usethis::use_data(example_protanar_1, overwrite = TRUE)
