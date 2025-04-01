library(ggplot2) 

### Custom ggplot theme
theme_set(theme_bw() +
            theme(plot.title = element_text(hjust = 0.5, size = 24),
                  plot.subtitle = element_text(hjust = 0.5, size = 18),
                  axis.title = element_text(size = 20),
                  strip.text = element_text(size = 12),
                  plot.caption = element_text(size = 10),
                  legend.text = element_text(size = 12),
                  legend.position = "bottom"))

### Custom table 1 functions to get , after the thousands
render.continuous <- function(x, ...) {
  with(stats.default(x, ...), c("",
                                "Mean (SD)"         = sprintf("%s (%s)",
                                                              signif_pad(MEAN,   3, big.mark=","),
                                                              signif_pad(SD,     3, big.mark=","))))
}

render.categorical <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y,
                                                                        sprintf("%s (%s%%)", prettyNum(FREQ, big.mark=","), PCT))))
}

render.strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N=%s)</span></span>", 
          label, prettyNum(n, big.mark=","))
}

render.strat.subj_trials <- function (label, n, ...) {
  x <- print(unlist(str_split(names(label), '\\.')))
  n_subj <- rep(NA, length(n)) 
  for(i in 1:length(n)) {
    n_subj[i] <- n_distinct(df1$subject_id[df1$surgery == x[2*i-1] & df1$elig == x[2*i]])
  }
  
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N Subject Trials=%s<br>N Subjects=%s)</span></span>", 
          label, prettyNum(n, big.mark=","),  prettyNum(n_subj, big.mark=","))
}

render.strat.subj_trials_durable <- function (label, n, ...) {
  x <- print(unlist(str_split(names(label), '\\.')))
  n_subj <- rep(NA, length(n)) 
  for(i in 1:length(n)) {
    n_subj[i] <- n_distinct(df_final$subject_id[df_final$surgery == x[2*i-1] & df_final$elig_status == x[2*i]])
  }
  
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N Subject Trials=%s<br>N Subjects=%s)</span></span>", 
          label, prettyNum(n, big.mark=","),  prettyNum(n_subj, big.mark=","))
}

winsorize <- function(x, q) {
  ### Get upper and lower quantiles of the data
  q_data <- quantile(x[!is.na(x)], q)
  lower <- q_data[1]
  upper <- q_data[2]
  
  ### Replace extreme values by quantiles
  x[which(x < lower)] <- lower
  x[which(x > upper)] <- upper
  
  return(x)
  
}

### Work Around for gtsave on the cluster
cluster_gtsave <- function(table, file_path, keep_html = F, ...) {
  file <- gsub('\\.html', '', gsub('\\.png', '', file_path))
  gtsave(table, paste0(file, '.html'))
  webshot::webshot(paste0(file, '.html'), paste0(file, '.png'), ...)
  if(!keep_html) {
    file.remove(paste0(file, '.html'))
  }
}
