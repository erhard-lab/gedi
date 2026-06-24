#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(ggplot2)
	library(cowplot)
})

# Subtracted Target ECDF
get_F_T_adj <- function(S_T, S_D_CDS, grid) {
  if (length(S_D_CDS) == 0 || length(S_T) <= length(S_D_CDS)) {
    return(ecdf(S_T)(grid))
  }
  n_t <- length(S_T)
  n_d <- length(S_D_CDS)
  
  F_raw <- (n_t * ecdf(S_T)(grid) - n_d * ecdf(S_D_CDS)(grid)) / (n_t - n_d)
  
  # Enforce monotonicity and bounds [0, 1]
  F_mono <- pmin(1, pmax(0, cummax(F_raw)))
  F_mono <- F_mono / max(F_mono[length(F_mono)], 1e-12)
  return(pmin(1, F_mono))
}



t=read.csv(file,check.names=FALSE)

t$Len <- nchar(t$Sequence)

if (!"UniqueGenome" %in% names(t)) t$UniqueGenome <- "All"

t$Group <- paste(t$UniqueGenome, t$Len, sep = "_")

idx_T     <- t$Annotation == "CDS" & t$Decoy != "D"
idx_D_CDS <- t$Annotation == "CDS" & t$Decoy == "D"
idx_D     <- t$Decoy == "D"



stats <- aggregate(Score ~ Group + UniqueGenome + Len, data = t, FUN = length)
names(stats)[names(stats) == "Score"] <- "N_Total"

counts_T <- aggregate(Score ~ Group, data = t[idx_T, ], FUN = length)
counts_D <- aggregate(Score ~ Group, data = t[idx_D, ], FUN = length)

stats <- merge(stats, counts_T, by = "Group", all.x = TRUE)
names(stats)[names(stats) == "Score"] <- "N_T"
stats <- merge(stats, counts_D, by = "Group", all.x = TRUE)
names(stats)[names(stats) == "Score"] <- "N_D"

stats$N_T[is.na(stats$N_T)] <- 0
stats$N_D[is.na(stats$N_D)] <- 0

stats$valid <- stats$N_T >= 100 & stats$N_D >= 100

valid_stats <- stats[stats$valid, ]
if (nrow(valid_stats) == 0) stop("No valid groups (>= 100 targets and >= 100 decoys).")

stats$mapped_Group <- stats$Group
for (i in seq_len(nrow(stats))) {
  if (!stats$valid[i]) {
    cur_gen <- stats$UniqueGenome[i]
    cur_len <- stats$Len[i]
    
    cands <- valid_stats[valid_stats$UniqueGenome == cur_gen, ]
    if (nrow(cands) == 0) cands <- valid_stats
    
    min_dist <- min(abs(cands$Len - cur_len))
    cands <- cands[abs(cands$Len - cur_len) == min_dist, ]
    best_cand <- cands[which.max(cands$N_Total)[1], ]
    
    stats$mapped_Group[i] <- best_cand$Group
  }
}

fallback_map <- setNames(stats$mapped_Group, stats$Group)



# ECDF Mixture Modeling
S_T_list     <- split(t$Score[idx_T], t$Group[idx_T])
S_D_CDS_list <- split(t$Score[idx_D_CDS], t$Group[idx_D_CDS])
S_D_list     <- split(t$Score[idx_D], t$Group[idx_D])

t$pi <- NA_real_
t$PEP <- NA_real_
t$Q <- NA_real_

groups <- split(seq_len(nrow(t)), list(t$Group, t$Annotation), drop = TRUE)

for (idx in groups) {
  grp_str <- t$Group[idx[1]]
  S_M <- t$Score[idx]
  
  if (length(S_M) == 0) next
  
  mapped_grp <- fallback_map[grp_str]
  S_T     <- S_T_list[[mapped_grp]]
  S_D_CDS <- S_D_CDS_list[[mapped_grp]]
  S_D     <- S_D_list[[mapped_grp]]
  
  grid <- seq(0, 1, length.out = 100)
  F_M <- ecdf(S_M)(grid)
  F_T <- get_F_T_adj(S_T, S_D_CDS, grid)
  F_D <- ecdf(S_D)(grid)
  
  loss <- function(p) sum((F_M - (p * F_T + (1 - p) * F_D))^2)
  pi_hat <- optimize(loss, interval = c(0, 1))$minimum
  
  tail_M <- sapply(S_M, function(s) mean(S_M >= s))
  tail_D <- sapply(S_M, function(s) mean(S_D >= s))
  fdr <- pmin(1, (1 - pi_hat) * tail_D / pmax(tail_M, 1e-12))
  
  ord <- order(S_M, decreasing = FALSE)
  Q <- numeric(length(S_M))
  Q[ord] <- cummin(fdr[ord])
  
  eps <- 0.05
  mass_M <- sapply(S_M, function(s) mean(abs(S_M - s) <= eps))
  mass_D <- sapply(S_M, function(s) mean(abs(S_D - s) <= eps))
  pep <- pmin(1, (1 - pi_hat) * mass_D / pmax(mass_M, 1e-12))
  
  t$pi[idx] <- pi_hat
  t$PEP[idx] <- pep
  t$Q[idx] <- Q
}

pdf(paste0(prefix,".pep.mixmod.pdf"), width = 8, height = 6)

group_meta <- unique(t[, c("UniqueGenome", "Len", "Group")])
unique_groups <- group_meta$Group[order(group_meta$UniqueGenome, group_meta$Len)]

for (g in unique_groups) {
  sub_t <- t[t$Group == g, ]
  
  mapped_g <- fallback_map[g]
  is_fallback <- g != mapped_g
  
  if (!is_fallback) {
  
	  S_T     <- S_T_list[[mapped_g]]
	  S_D_CDS <- S_D_CDS_list[[mapped_g]]
	  S_D     <- S_D_list[[mapped_g]]
	  
	  if (length(S_T) < 100 || length(S_D) < 100) next
	  
	  gen_val <- sub_t$UniqueGenome[1]
	  len_val <- sub_t$Len[1]
	  
	  title_prefix <- sprintf("Genome: %s | Length: %s", gen_val, len_val)
	  
	  grid_plot <- seq(0, 1, length.out = 500)
	  F_T_plot <- get_F_T_adj(S_T, S_D_CDS, grid_plot)
	  F_D_plot <- ecdf(S_D)(grid_plot)
	  
	  df_td <- data.frame(
	    Score = rep(grid_plot, 2),
	    Fx = c(F_T_plot, F_D_plot),
	    Class = rep(c("Target (Subtracted)", "Decoy"), each = 500)
	  )
	  
	  p1 <- ggplot(df_td, aes(x = Score, y = Fx, color = Class)) +
	    geom_step(linewidth = 1) +
	    scale_color_manual(values = c("Target (Subtracted)" = "forestgreen", "Decoy" = "firebrick")) +
	    labs(title = paste(title_prefix, "\nTarget vs Decoy ECDFs"), y = "F(x)") +
	    theme_cowplot()
	  print(p1)
	  
	  cats <- unique(sub_t$Annotation)
	  
	  for (c in cats) {
	    idx <- sub_t$Annotation == c
	    S_M <- sub_t$Score[idx]
	    pi_hat <- sub_t$pi[idx][1]
	    
	    if (is.na(pi_hat) || length(S_M) == 0) next
	    
	    qvals <- sub_t$Q[idx]
	    fdr_pass <- S_M[which(qvals <= 0.10)]
	    cutoff <- if (length(fdr_pass) > 0) min(fdr_pass) else NA
	    
	    F_fit <- pi_hat * F_T_plot + (1 - pi_hat) * F_D_plot
	    
	    df_fit <- data.frame(Score = grid_plot, Fx = F_fit,T=F_T_plot,D=F_D_plot)
	    df_obs <- data.frame(Score = S_M)
	    
	    p2 <- ggplot() +
	      stat_ecdf(data = df_obs, aes(x = Score, color = "Observed"), geom = "step", linewidth = 1) +
	      geom_line(data = df_fit, aes(x = Score, y = Fx, color = "Fitted Mixture"), linewidth = 1, linetype = "dashed") +
	      geom_line(data = df_fit, aes(x = Score, y = T, color = "Target"), linewidth = 1, linetype = "dotted") +
	      geom_line(data = df_fit, aes(x = Score, y = D, color = "Decoy"), linewidth = 1, linetype = "dotted") +
	      scale_color_manual(
	        name = "Distribution",
	        values = c("Observed" = "black", "Fitted Mixture" = "dodgerblue", "Target" = "forestgreen", "Decoy" = "firebrick")
	      ) +
	      labs(
	        title = paste(title_prefix, "| Annotation:", c),
	        subtitle = sprintf("Estimated targets: n=%.0f (%.1f%%)", pi_hat*sum(idx), pi_hat*100),
	        y = "F(x)"
	      ) +
	      theme_cowplot()
	    
	    if (!is.na(cutoff)) {
	      p2 <- p2 + 
	        geom_vline(xintercept = cutoff, linetype = "solid", color = "darkorange", linewidth = 0.5) +
	        annotate("text", x = cutoff, y = 0.05, label = "10% FDR", color = "darkorange", 
	                 angle = 90, vjust = -0.5)
	    }
	    print(p2)
	  }
	}
}
invisible(dev.off())




pdf(paste0(prefix,".pep.mixcomp.pdf"), width = 8, height = 6)
idx_T     <- t$Annotation == "CDS" & t$Decoy != "D"
idx_D     <- t$Decoy == "D"

counts_T <- aggregate(Score ~ Group, data = t[idx_T, ], FUN = length)
counts_D <- aggregate(Score ~ Group, data = t[idx_D, ], FUN = length)

valid_groups <- intersect(
  counts_T$Group[counts_T$Score >= 100], 
  counts_D$Group[counts_D$Score >= 100]
)

plot_data <- list()
grid_plot <- seq(0, 1, length.out = 500)

for (g in valid_groups) {
  sub_t <- t[t$Group == g, ]
  
  gen_val <- sub_t$UniqueGenome[1]
  len_val <- as.numeric(sub_t$Len[1]) # Numeric for continuous gradient
  
  S_T     <- sub_t$Score[sub_t$Annotation == "CDS" & sub_t$Decoy != "D"]
  S_D_CDS <- sub_t$Score[sub_t$Annotation == "CDS" & sub_t$Decoy == "D"]
  S_D     <- sub_t$Score[sub_t$Decoy == "D"]
  
  # Target (Subtracted)
  F_T_plot <- get_F_T_adj(S_T, S_D_CDS, grid_plot)
  plot_data[[length(plot_data) + 1]] <- data.frame(
    UniqueGenome = gen_val, Len = len_val, Class = "Target (Subtracted)",
    Score = grid_plot, Fx = F_T_plot
  )
  
  # Decoy
  F_D_plot <- ecdf(S_D)(grid_plot)
  plot_data[[length(plot_data) + 1]] <- data.frame(
    UniqueGenome = gen_val, Len = len_val, Class = "Decoy",
    Score = grid_plot, Fx = F_D_plot
  )
}

df_plot <- do.call(rbind, plot_data)
genomes <- sort(unique(df_plot$UniqueGenome))
classes <- c("Target (Subtracted)", "Decoy")

for (g in genomes) {
  for (cls in classes) {
    sub_df <- df_plot[df_plot$UniqueGenome == g & df_plot$Class == cls, ]
    
    if (nrow(sub_df) == 0) next
    
    p <- ggplot(sub_df, aes(x = Score, y = Fx, color = Len, group = Len)) +
      geom_step(linewidth = 0.8) +
      scale_color_viridis_c(option = "plasma") + 
      labs(
        title = sprintf("Genome: %s | Class: %s", g, cls),
        x = "Score", 
        y = "F(x)", 
        color = "Length"
      ) +
      theme_cowplot()
    
    print(p)
  }
}

invisible(dev.off())

t$Len=NULL
t$Group=NULL
t$pi=NULL
write.table(t,file,row.names=FALSE,col.names=TRUE,sep=",",quote=FALSE)

