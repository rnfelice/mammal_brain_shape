
getStones <- function(stones, labels = NULL) {

  res <- matrix(ncol = 2, nrow = length(stones))
  if (!is.null(labels)) {
    if (length(labels) != length(stones)) {
      stop("The number of labels must equal the number of stones files.")
    }
  }
  colnames(res) <- c("logfile", "marginalLh")
  for (i in 1:length(stones)) {
    raw <- readLines(stones[[i]])
    res[i, 1] <- stones[[i]]
    res[i, 2] <- as.numeric(strsplit(raw[length(raw)], "\t")[[1]][2])
  }
  res <- data.frame(res, stringsAsFactors = FALSE)
  if (!is.null(labels)) {
    res$logfile <- labels
  }
  res$marginalLh <- as.numeric(as.character(res$marginalLh))
  res <- tibble::as_tibble(res)
  class(res) <- append("bt_stones", class(res))
  return(res)
}

plot.bt_stones <- function(x) {
  ggplot2::autoplot(x)
}

autoplot.bt_stones <- function(x) {
  d <- matrix(ncol = nrow(x), nrow = nrow(x))
  colnames(d) <- rownames(d) <- x$logfile
  
  for (i in seq_len(nrow(d))) {
    for (j in seq_len(ncol(d))) {
      d[i, j] <- round(2 * (x$marginalLh[i] - x$marginalLh[j]), 2)
    }
  }
  
  d <- reshape2::melt(d[nrow(d):1, ])
  #we use the 2ln(BF) scale  from Kass and Raferty 1995
  #because this puts the values on the same scale as likelihood ratios 
  
  d$bf_cat <- cut(d$value, 
                  breaks = c(min(d$value) - 1, 0, 2, 5, 10, max(d$value)),
                  labels = c("<0", "0-2", "2-5", "5-10", ">10"))
  d$bf_cat <- factor(d$bf_cat, levels = c("<0", "0-2", "2-5", "5-10", ">10"))
  d$value[d$value == 0] <- d$bf_cat[d$value == 0] <- NA
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Var2, y = Var1, fill = bf_cat)) +
    ggplot2::geom_tile(colour = "white", size = 0.25, show.legend = TRUE) +
    ggplot2::geom_text(data = d, ggplot2::aes(Var2, Var1, label = value), 
                       na.rm = TRUE) +
    ggplot2::scale_y_discrete(expand=c(0,0))+
    ggplot2::scale_x_discrete(expand=c(0,0), position = "top") +
    ggplot2::scale_fill_manual(
      values = c( viridis::viridis(8)[3:7]),
      na.value     = "grey90",
      drop         = FALSE,
      limits       = c("<0", "0-2", "2-5", "5-10", ">10")
    )+
    ggplot2::labs(x = "", y = "", fill = "2 × log(Bayes factor)") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_family = "Helvetica")
  p
}