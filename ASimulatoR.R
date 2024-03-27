library(ASimulatoR)

input_dir = "/home/miyokawa/ドキュメント/ngs/araport/RSim"
num_reps = c(6,6)
outdir = "/media/miyokawa/8TB-Data3/RSim_out_sp"
max_genes = NULL
probs_as_freq = T
event_freq = setNames(rep(1/4, 4), c('es', 'ir', 'a3', 'a5'))

seed = 110
ncores = 4
readlen = 150
seq_depth = 6e07f
error_rate = 0.001
             
simulate_alternative_splicing(seed = seed,
                              ncores = ncores,
                              readlen = readlen,
                              seq_depth = seq_depth,
                              input_dir = input_dir,
                              outdir = outdir, 
                              event_probs = event_freq,
                              probs_as_freq = probs_as_freq, 
                              max_genes = max_genes,
                              error_rate = error_rate,
                              num_reps = num_reps,
                              verbose = TRUE)
