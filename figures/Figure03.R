
# SNPs -> CE (exposure) -> Meniere (outcome) --------------------------------------------------------------

load("GCST90025818R11H8MENIERE_1e6_MAGMArefclump.RData")

mr_leaveoneout_plot(exp_out_loo)[[1]] +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 6))

#forest
p <- ggforestplot::forestplot(
  df = exp_out_res,
  name = outcome,
  estimate = b,
  se = se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR (95% CI)",
  colour = method,
  shape = method,
  logodds = T) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))


# SNPs -> Meniere (exposure) -> CE (outcome) ---------------------------------
load("R11H8MENIEREGCST90025818_1e6_MAGMArefclump.RData")

mr_leaveoneout_plot(exp_out_loo)[[1]] +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 6))

#forest
p <- ggforestplot::forestplot(
  df = exp_out_res,
  name = outcome,
  estimate = b,
  se = se,
  pvalue = pval,
  psignif = 0.05,
  xlab = "OR (95% CI)",
  colour = method,
  shape = method,
  logodds = T) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 18))
