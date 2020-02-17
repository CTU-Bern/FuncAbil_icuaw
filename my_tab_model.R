my_tab_model <- 
function (..., transform, show.intercept = TRUE, show.est = TRUE, 
          show.ci = 0.95, show.hdi50 = TRUE, show.se = NULL, show.std = NULL, 
          show.p = TRUE, show.stat = FALSE, show.df = FALSE, show.zeroinf = TRUE, 
          show.r2 = TRUE, show.icc = TRUE, show.adj.icc = FALSE, show.re.var = TRUE, 
          show.fstat = FALSE, show.aic = FALSE, show.aicc = FALSE, 
          show.dev = FALSE, show.obs = TRUE, terms = NULL, rm.terms = NULL, 
          group.terms = TRUE, order.terms = NULL, title = NULL, pred.labels = NULL, 
          dv.labels = NULL, wrap.labels = 25, string.pred = "Predictors", 
          string.std = "std. Beta", string.ci = "CI", string.se = "std. Error", 
          string.p = "p", string.df = "df", string.stat = "Statistic", 
          ci.hyphen = "&nbsp;to&nbsp;", minus.sign = "&#45;", # ci.hyphen = "&nbsp;&ndash;&nbsp;"
          collapse.ci = FALSE, collapse.se = FALSE, linebreak = TRUE, 
          col.order = c("est", "se", "std.est", "std.se", "ci", "std.ci", 
                        "hdi.inner", "hdi.outer", "stat", "p", "df", "response.level"), 
          digits = 2, digits.p = 3, emph.p = TRUE, p.val = c("wald", 
                                                             "kr"), p.style = c("numeric", "asterisk"), p.threshold = c(0.05, 
                                                                                                                        0.01, 0.001), case = "parsed", auto.label = TRUE, prefix.labels = c("none", 
                                                                                                                                                                                            "varname", "label"), bpe = "median", CSS = css_theme("regression"), 
          file = NULL, use.viewer = TRUE) 
{
  p.val <- match.arg(p.val)
  p.style <- match.arg(p.style)
  prefix.labels <- match.arg(prefix.labels)
  if (missing(case)) {
    if (prefix.labels == "none") 
      case <- "parsed"
    else case <- NULL
  }
  if (p.style == "asterisk") 
    show.p <- FALSE
  show.se <- sjPlot:::check_se_argument(se = show.se, type = NULL)
  models <- list(...)
  if (length(class(models[[1]])) == 1 && class(models[[1]]) == 
      "list") 
    models <- lapply(models[[1]], function(x) x)
  names(models) <- unlist(lapply(match.call(expand.dots = F)$..., 
                                 function(.x) deparse(.x, width.cutoff = 500L)))
  auto.transform <- missing(transform)
  ci.lvl <- ifelse(is.null(show.ci), 0.95, show.ci)
  copos <- which("est" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "estimate"
  copos <- which("se" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "std.error"
  copos <- which("ci" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "conf.int"
  copos <- which("std.est" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "std.estimate"
  copos <- which("std.ci" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "std.conf.int"
  copos <- which("p" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "p.value"
  copos <- which("stat" == col.order)
  if (!sjmisc::is_empty(copos)) 
    col.order[copos] <- "statistic"
  model.list <- purrr::map2(models, 1:length(models), function(model, 
                                                               i) {
    fam.info <- sjstats::model_family(model)
    if (auto.transform) {
      if (fam.info$is_linear) 
        transform <- NULL
      else transform <- "exp"
    }
    dat <- sjPlot:::tidy_model(model, ci.lvl = ci.lvl, transform, 
                      type = "est", bpe, se = show.se, facets = FALSE, 
                      show.zeroinf = show.zeroinf, p.val = p.val, robust = NULL)
    if (!sjPlot:::is.stan(model) && !is.null(transform)) {
      funtrans <- match.fun(transform)
      dat[["estimate"]] <- funtrans(dat[["estimate"]])
      dat[["conf.low"]] <- funtrans(dat[["conf.low"]])
      dat[["conf.high"]] <- funtrans(dat[["conf.high"]])
    }
    dat <- dat %>% dplyr::mutate(conf.int = sprintf("%.*f%s%.*f", 
                                                    digits, .data$conf.low, ci.hyphen, digits, .data$conf.high)) %>% 
      dplyr::select(-.data$conf.low, -.data$conf.high) %>% 
      dplyr::mutate(p.stars = sjPlot:::get_p_stars(.data$p.value, 
                                          p.threshold), p.sig = .data$p.value < 0.05, p.value = sprintf("%.*f", 
                                                                                                        digits.p, .data$p.value))
    if (emph.p && !all(dat$p.value == "NA")) 
      dat$p.value[dat$p.sig] <- sprintf("<strong>%s</strong>", 
                                        dat$p.value[dat$p.sig])
    dat <- dplyr::select(dat, -.data$p.sig)
    if (sjPlot:::is.stan(model)) {
      dat <- dat %>% sjmisc::var_rename(conf.int = "hdi.outer") %>% 
        dplyr::mutate(hdi.inner = sprintf("%.*f%s%.*f", 
                                          digits, .data$conf.low50, ci.hyphen, digits, 
                                          .data$conf.high50)) %>% dplyr::select(-.data$conf.low50, 
                                                                                -.data$conf.high50)
    }
    pv <- paste0("0.", paste(rep("0", digits.p), collapse = ""))
    dat$p.value[dat$p.value == pv] <- "&lt;0.001"
    pv <- paste0("<strong>0.", paste(rep("0", digits.p), 
                                     collapse = ""), "</strong>")
    dat$p.value[dat$p.value == pv] <- "<strong>&lt;0.001"
    if (!is.null(show.std) && fam.info$is_linear && !sjPlot:::is.stan(model)) {
      dat <- model %>% sjstats::std_beta(type = show.std, 
                                         ci.lvl = ci.lvl) %>% sjmisc::var_rename(std.error = "std.se", 
                                                                                 conf.low = "std.conf.low", conf.high = "std.conf.high") %>% 
        sjmisc::add_case(.after = -1) %>% dplyr::select(-1) %>% 
        sjmisc::add_columns(dat) %>% dplyr::mutate(std.conf.int = sprintf("%.*f%s%.*f", 
                                                                          digits, .data$std.conf.low, ci.hyphen, digits, 
                                                                          .data$std.conf.high)) %>% dplyr::select(-.data$std.conf.low, 
                                                                                                                  -.data$std.conf.high)
    }
    if (p.style == "asterisk") {
      if (obj_has_name(dat, "estimate")) 
        dat$estimate <- sprintf("%.*f <sup>%s</sup>", 
                                digits, dat$estimate, dat$p.stars)
      if (!show.est && obj_has_name(dat, "std.estimate")) 
        dat$std.estimate <- sprintf("%.*f <sup>%s</sup>", 
                                    digits, dat$std.estimate, dat$p.stars)
    }
    dat <- dplyr::select(dat, -.data$p.stars)
    dat <- dat[, sjPlot:::sort_columns(colnames(dat), sjPlot:::is.stan(model), 
                              col.order)]
    cn <- colnames(dat)[2:ncol(dat)]
    colnames(dat)[2:ncol(dat)] <- sprintf("%s_%i", cn, i)
    dat <- dat %>% purrr::map_if(is.numeric, ~sprintf("%.*f", 
                                                      digits, .x)) %>% as.data.frame(stringsAsFactors = FALSE)
    if (!show.hdi50) 
      dat <- dplyr::select(dat, -sjPlot:::string_starts_with("hdi.inner", 
                                                    colnames(dat)))
    if (collapse.ci) {
      if (linebreak) 
        lb <- "<br>"
      else lb <- " "
      est.cols <- string_starts_with("estimate", x = colnames(dat))
      dat[[est.cols]] <- sprintf("%s%s(%s)", dat[[est.cols]], 
                                 lb, dat[[est.cols + 2]])
      if (!sjmisc::is_empty(string_starts_with("hdi", x = colnames(dat)))) {
        dat <- dplyr::select(dat, -sjPlot:::string_starts_with("hdi.outer", 
                                                      x = colnames(dat)))
        dat[[est.cols]] <- sprintf("%s%s(%s)", dat[[est.cols]], 
                                   lb, dat[[est.cols + 2]])
        dat <- dplyr::select(dat, -sjPlot:::string_starts_with("hdi.inner", 
                                                      x = colnames(dat)))
      }
      else {
        dat <- dplyr::select(dat, -sjPlot:::string_starts_with("conf.int", 
                                                      x = colnames(dat)))
      }
      std.cols <- sjPlot:::string_starts_with("std.estimate", x = colnames(dat))
      if (!sjmisc::is_empty(std.cols)) {
        dat[[std.cols]] <- sprintf("%s%s(%s)", dat[[std.cols]], 
                                   lb, dat[[std.cols + 2]])
        dat <- dplyr::select(dat, -sjPlot:::string_starts_with("std.conf.int", 
                                                      x = colnames(dat)))
      }
    }
    if (collapse.se) {
      if (linebreak) 
        lb <- "<br>"
      else lb <- " "
      est.cols <- sjPlot:::string_starts_with("estimate", x = colnames(dat))
      dat[[est.cols]] <- sprintf("%s%s(%s)", dat[[est.cols]], 
                                 lb, dat[[est.cols + 1]])
      dat <- dplyr::select(dat, -sjPlot:::string_starts_with("std.error", 
                                                    x = colnames(dat)))
      std.cols <- sjPlot:::string_starts_with("std.estimate", x = colnames(dat))
      if (!sjmisc::is_empty(std.cols)) {
        dat[[std.cols]] <- sprintf("%s%s(%s)", dat[[std.cols]], 
                                   lb, dat[[std.cols + 1]])
        dat <- dplyr::select(dat, -sjPlot:::string_starts_with("std.se", 
                                                      x = colnames(dat)))
      }
    }
    zidat <- NULL
    wf <- sjPlot:::string_starts_with("wrap.facet", x = colnames(dat))
    if (!sjmisc::is_empty(wf)) {
      zi <- which(dat[[wf]] %in% c("Zero-Inflated Model", 
                                   "Zero Inflation Model"))
      if (show.zeroinf && !sjmisc::is_empty(zi)) {
        zidat <- dat %>% dplyr::slice(!(!zi)) %>% dplyr::select(!(!-wf))
      }
      if (!sjmisc::is_empty(zi)) 
        dat <- dplyr::slice(dat, !(!-zi))
      dat <- dplyr::select(dat, !(!-wf))
    }
    n_obs <- NULL
    if (show.obs) {
      n_obs <- tryCatch({
        stats::nobs(model)
      }, error = function(x) {
        NULL
      })
    }
    icc <- NULL
    if ((show.icc || show.re.var) && sjPlot:::is_mixed_model(model)) {
      icc <- tryCatch({
        suppressWarnings(sjstats::icc(model))
      }, error = function(x) {
        NULL
      })
    }
    icc.adjusted <- NULL
    if ((show.adj.icc) && sjPlot:::is_mixed_model(model)) {
      icc.adjusted <- tryCatch({
        sjstats::icc(model, adjusted = TRUE, type = "all")
      }, error = function(x) {
        NULL
      })
    }
    rsq <- NULL
    if (show.r2) {
      if (sjPlot:::is_mixed_model(model) && !is.null(icc.adjusted)) {
        rsq <- icc.adjusted[[1]]
      }
      else {
        rsq <- tryCatch({
          suppressWarnings(sjstats::r2(model))
        }, error = function(x) {
          NULL
        })
      }
    }
    if (!is.null(icc.adjusted)) 
      icc.adjusted <- icc.adjusted[[2]]
    dev <- NULL
    if (show.dev) 
      dev <- model_deviance(model)
    aic <- NULL
    if (show.aic) 
      aic <- model_aic(model)
    if (inherits(model, "brmsfit")) {
      dat$term <- gsub("^b_", "", dat$term)
      if (!is.null(zidat)) 
        zidat$term <- gsub("^b_", "", zidat$term)
    }
    list(dat = dat, transform = transform, zeroinf = zidat, 
         rsq = rsq, n_obs = n_obs, icc = icc, dev = dev, aic = aic, 
         icc.adj = icc.adjusted)
  })
  na.vals <- c("NA", sprintf("NA%sNA", ci.hyphen), sprintf("NA (NA%sNA)", 
                                                           ci.hyphen), sprintf("NA (NA%sNA) (NA)", ci.hyphen))
  model.data <- purrr::map(model.list, ~.x[[1]])
  transform.data <- purrr::map(model.list, ~.x[[2]])
  zeroinf.data <- purrr::map(model.list, ~.x[[3]])
  rsq.data <- purrr::map(model.list, ~.x[[4]])
  n_obs.data <- purrr::map(model.list, ~.x[[5]])
  icc.data <- purrr::map(model.list, ~.x[[6]])
  dev.data <- purrr::map(model.list, ~.x[[7]])
  aic.data <- purrr::map(model.list, ~.x[[8]])
  icc.adj.data <- purrr::map(model.list, ~.x[[9]])
  is.zeroinf <- purrr::map_lgl(model.list, ~!is.null(.x[[3]]))
  zeroinf.data <- purrr::compact(zeroinf.data)
  if (!show.zeroinf) 
    zeroinf.data <- NULL
  model.data <- purrr::map(model.data, function(.x) {
    resp.col <- sjPlot:::string_starts_with("response.level", x = colnames(.x))
    if (!sjmisc::is_empty(resp.col)) 
      .x[order(match(.x[[resp.col]], unique(.x[[resp.col]]))), 
         ]
    else .x
  })
  show.response <- TRUE
  if (length(model.data) == 1) {
    fi <- sjstats::model_family(models[[1]])
    if (fi$is_multivariate || fi$is_categorical) {
      show.response <- FALSE
      if (fi$is_categorical) {
        dv.labels <- sprintf("%s: %s", sjstats::resp_var(models[[1]]), 
                             unique(model.data[[1]][["response.level_1"]]))
        model.data <- split(model.data[[1]], model.data[[1]]["response.level_1"])
      }
      else {
        dv.labels <- sjmisc::word_wrap(sjlabelled::get_dv_labels(models, 
                                                                 mv = TRUE, case = case), wrap = wrap.labels, 
                                       linesep = "<br>")
        if (sjmisc::is_empty(dv.labels) || !isTRUE(auto.label)) 
          dv.labels <- sjstats::resp_var(models[[1]])
        model.data <- split(model.data[[1]], model.data[[1]]["response.level_1"])
        dv.labels <- dv.labels[match(names(dv.labels), 
                                     names(model.data))]
      }
      model.data <- purrr::map2(model.data, 1:length(model.data), 
                                function(x, y) {
                                  colnames(x) <- gsub(pattern = "_1", replacement = sprintf("_%i", 
                                                                                            y), x = colnames(x))
                                  x
                                })
    }
  }
  dat <- model.data %>% purrr::reduce(~dplyr::full_join(.x, 
                                                        .y, by = "term")) %>% purrr::map_df(~dplyr::if_else(.x %in% 
                                                                                                              na.vals | is.na(.x), "", .x))
  dat <- sjPlot:::remove_unwanted(dat, show.intercept, show.est, show.std, 
                         show.ci, show.se, show.stat, show.p, show.df, show.response, 
                         terms, rm.terms)
  zeroinf <- NULL
  if (!sjmisc::is_empty(zeroinf.data)) {
    zeroinf <- zeroinf.data %>% purrr::reduce(~dplyr::full_join(.x, 
                                                                .y, by = "term")) %>% purrr::map_df(~dplyr::if_else(.x %in% 
                                                                                                                      na.vals | is.na(.x), "", .x))
    zeroinf <- sjPlot:::remove_unwanted(zeroinf, show.intercept, show.est, 
                               show.std, show.ci, show.se, show.stat, show.p, show.df, 
                               show.response, terms, rm.terms)
  }
  if (isTRUE(auto.label) && sjmisc::is_empty(pred.labels)) {
    pred.labels <- sjlabelled::get_term_labels(models, case = case, 
                                               mark.cat = TRUE, prefix = prefix.labels)
    no.dupes <- !duplicated(names(pred.labels))
    pred.labels <- prepare.labels(pred.labels[no.dupes], 
                                  grp = group.terms)
  }
  else {
    group.terms <- FALSE
  }
  if (!sjmisc::is_empty(pred.labels)) {
    if (!is.null(names(pred.labels))) {
      labs <- sjmisc::word_wrap(pred.labels, wrap = wrap.labels, 
                                linesep = "<br>")
      tr <- 1:nrow(dat)
      find.matches <- match(dat$term, names(pred.labels))
      find.na <- which(is.na(find.matches))
      if (!sjmisc::is_empty(find.na)) 
        tr <- tr[-find.na]
      rp <- as.vector(stats::na.omit(find.matches))
      dat$term[tr] <- unname(labs[rp])
      if (!is.null(zeroinf)) {
        tr <- 1:nrow(zeroinf)
        find.matches <- match(zeroinf$term, names(pred.labels))
        find.na <- which(is.na(find.matches))
        if (!sjmisc::is_empty(find.na)) 
          tr <- tr[-find.na]
        rp <- as.vector(stats::na.omit(find.matches))
        zeroinf$term[tr] <- unname(labs[rp])
      }
    }
    else {
      if (length(pred.labels) == nrow(dat)) 
        dat$term <- pred.labels
      else message("Length of `pred.labels` does not equal number of predictors, no labelling applied.")
    }
  }
  if (isTRUE(auto.label) && sjmisc::is_empty(dv.labels)) {
    dv.labels <- sjmisc::word_wrap(sjlabelled::get_dv_labels(models, 
                                                             case = case), wrap = wrap.labels, linesep = "<br>")
  }
  else if (sjmisc::is_empty(dv.labels)) {
    dv.labels <- purrr::map(models, sjstats::resp_var) %>% 
      purrr::flatten_chr()
  }
  if (!is.null(order.terms)) {
    if (length(order.terms) == nrow(dat)) {
      dat <- dat[order.terms, ]
    }
    else {
      message("Number of values in `order.terms` does not match number of terms. Terms are not sorted.")
    }
  }
  col.header <- purrr::map_chr(colnames(dat), function(x) {
    pos <- grep("^estimate_", x)
    if (!sjmisc::is_empty(pos)) {
      i <- as.numeric(sub("estimate_", "", x = x, fixed = T))
      if (i <= length(models)) {
        x <- sjPlot:::estimate_axis_title(models[[i]], axis.title = NULL, 
                                 type = "est", transform = transform.data[[i]], 
                                 multi.resp = NULL, include.zeroinf = FALSE)
      }
      else if (length(models) == 1) {
        fi <- sjstats::model_family(models[[1]])
        if (fi$is_multivariate) 
          mr <- i
        else mr <- NULL
        x <- sjPlot:::estimate_axis_title(models[[1]], axis.title = NULL, 
                                 type = "est", transform = transform.data[[1]], 
                                 multi.resp = mr, include.zeroinf = FALSE)
      }
      else {
        x <- "Estimate"
      }
    }
    pos <- grep("^term", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.pred
    pos <- grep("^conf.int", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.ci
    pos <- grep("^std.error", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.se
    pos <- grep("^std.estimate", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.std
    pos <- grep("^std.se", x)
    if (!sjmisc::is_empty(pos)) 
      x <- paste("standardized", string.se)
    pos <- grep("^std.conf.int", x)
    if (!sjmisc::is_empty(pos)) 
      x <- paste("standardized", string.ci)
    pos <- grep("^p.value", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.p
    pos <- grep("^df", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.df
    pos <- grep("^statistic", x)
    if (!sjmisc::is_empty(pos)) 
      x <- string.stat
    pos <- grep("^response.level", x)
    if (!sjmisc::is_empty(pos)) 
      x <- "Response"
    pos <- grep("^hdi.inner", x)
    if (!sjmisc::is_empty(pos)) 
      x <- "HDI (50%)"
    pos <- grep("^hdi.outer", x)
    if (!sjmisc::is_empty(pos)) 
      x <- sprintf("HDI (%i%%)", round(100 * show.ci))
    x
  })
  if (p.style == "asterisk") 
    footnote <- sprintf("* p&lt;%s&nbsp;&nbsp;&nbsp;** p&lt;%s&nbsp;&nbsp;&nbsp;*** p&lt;%s", 
                        format(p.threshold[1]), format(p.threshold[2]), format(p.threshold[3]))
  else footnote <- NULL
  # tab_model_df(x = dat, zeroinf = zeroinf, is.zeroinf = is.zeroinf, 
  #              title = title, col.header = col.header, dv.labels = dv.labels, 
  #              rsq.list = rsq.data, n_obs.list = n_obs.data, icc.list = icc.data, 
  #              icc.adj.list = icc.adj.data, dev.list = dev.data, aic.list = aic.data, 
  #              n.models = length(model.list), show.re.var = show.re.var, 
  #              show.icc = show.icc, show.adj.icc = show.adj.icc, CSS = CSS, 
  #              file = file, use.viewer = use.viewer, footnote = footnote)
  return(as.data.frame(dat))
}


##### function for qrjoint plot ######


my.coef.qrjoint <- function (object, burn.perc = 0.5, nmc = 200, plot = FALSE, show.intercept = TRUE, 
          reduce = TRUE, ...) 
{
  nsamp <- object$dim[10]
  pars <- matrix(object$parsamp, ncol = nsamp)
  ss <- unique(round(nsamp * seq(burn.perc, 1, len = nmc + 
                                   1)[-1]))
  n <- object$dim[1]
  p <- object$dim[2]
  L <- object$dim[3]
  mid <- object$dim[4] + 1
  nknots <- object$dim[5]
  ngrid <- object$dim[6]
  a.sig <- object$hyper[1:2]
  a.kap <- matrix(object$hyper[-c(1:2)], nrow = 3)
  tau.g <- object$tau.g
  reg.ix <- object$reg.ix
  x.ce <- outer(rep(1, L), attr(object$x, "scaled:center"))
  x.sc <- outer(rep(1, L), attr(object$x, "scaled:scale"))
  base.bundle <- list()
  if (object$fbase.choice == 1) {
    base.bundle$q0 <- function(u, nu = Inf) return(1/(dt(qt(unitFn(u), 
                                                            df = nu), df = nu) * qt(0.9, df = nu)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qt(unitFn(u), 
                                                      df = nu)/qt(0.9, df = nu))
    base.bundle$F0 <- function(x, nu = Inf) return(pt(x * 
                                                        qt(0.9, df = nu), df = nu))
  }
  else if (object$fbase.choice == 2) {
    base.bundle$q0 <- function(u, nu = Inf) return(1/dlogis(qlogis(unitFn(u))))
    base.bundle$Q0 <- function(u, nu = Inf) return(qlogis(unitFn(u)))
    base.bundle$F0 <- function(x, nu = Inf) return(plogis(x))
  }
  else {
    base.bundle$q0 <- function(u, nu = Inf) return(1/(dunif(qunif(u, 
                                                                  -1, 1), -1, 1)))
    base.bundle$Q0 <- function(u, nu = Inf) return(qunif(u, 
                                                         -1, 1))
    base.bundle$F0 <- function(x, nu = Inf) return(punif(x, 
                                                         -1, 1))
  }
  beta.samp <- sapply(ss, function(p1) estFn(pars[, p1], object$x, 
                                             object$y, object$gridmats, L, mid, nknots, ngrid, a.kap, 
                                             a.sig, tau.g, reg.ix, reduce, x.ce, x.sc, base.bundle), 
                      simplify = "array")
  if (reduce) 
    tau.g <- tau.g[reg.ix]
  L <- length(tau.g)
  if (plot) {
    nr <- ceiling(sqrt(p + show.intercept))
    nc <- ceiling((p + show.intercept)/nr)
    cur.par <- par(no.readonly = TRUE)
    par(mfrow = c(nr, nc))
  }
  beta.hat <- array(0, c(L, p + 1, 3))
  plot.titles <- c("Intercept", object$xnames)
  beta.hat[, 1, ] <- getBands(beta.samp[, 1, ], plot = (plot & 
                                                          show.intercept), add = FALSE, x = tau.g, xlab = "tau", 
                              ylab = "Coefficient", bty = "n")
  if (plot & show.intercept) 
    title(main = plot.titles[1])
  for (j in 2:(p + 1)) {
    beta.hat[, j, ] <- getBands(beta.samp[, j, ], plot = plot, 
                                add = FALSE, x = tau.g, xlab = "tau", ylab = "Coefficient", 
                                bty = "n")
    if (plot) {
      title(main = plot.titles[j])
      abline(h = 0, lty = 2, col = 4)
    }
  }
  if (plot) 
    suppressWarnings(par(cur.par, no.readonly = TRUE))
  dimnames(beta.hat) <- list(tau = tau.g, beta = plot.titles, 
                             summary = c("b.lo", "b.med", "b.hi"))
  dimnames(beta.samp) <- list(tau = tau.g, beta = plot.titles, 
                              iter = 1:length(ss))
  invisible(list(beta.samp = beta.samp, beta.est = beta.hat))
  mid.red <- which(tau.g == object$tau.g[mid])
  parametric.list <- rbind(beta.samp[mid.red, , , drop = TRUE], 
                           sigma = sigFn(pars[nknots * (p + 1) + (p + 1) + 1, ss], 
                                         a.sig), nu = nuFn(pars[(nknots + 1) * (p + 1) + 2, 
                                                                ss]))
  dimnames(parametric.list)[[1]][1 + 0:p] <- c("Intercept", 
                                               object$xnames)
  gamsignu <- t(apply(parametric.list, 1, quantile, pr = c(0.5, 
                                                           0.025, 0.975)))
  dimnames(gamsignu)[[2]] <- c("Estimate", "Lo95%", "Up95%")
  invisible(list(beta.samp = beta.samp, beta.est = beta.hat, 
                 parametric = gamsignu))
}

getBands <- function (b, col = 2, lwd = 1, plot = TRUE, add = FALSE, 
                      x = seq(0, 1, len = nrow(b)), remove.edges = TRUE, ...) 
{
  colRGB <- col2rgb(col)/255
  colTrans <- rgb(colRGB[1], colRGB[2], colRGB[3], alpha = 0.2)
  b.med <- apply(b, 1, quantile, pr = 0.5)
  b.lo <- apply(b, 1, quantile, pr = 0.025)
  b.hi <- apply(b, 1, quantile, pr = 1 - 0.025)
  L <- nrow(b)
  ss <- 1:L
  ss.rev <- L:1
  if (remove.edges) {
    ss <- 2:(L - 1)
    ss.rev <- (L - 1):2
  }
  if (plot) {
    if (!add) 
      plot(x[ss], b.med[ss], ty = "n", ylim = range(c(b.lo[ss], 
                                                      b.hi[ss])), ...)
    polygon(x[c(ss, ss.rev)], c(b.lo[ss], b.hi[ss.rev]), 
            col = colTrans, border = colTrans)
    lines(x[ss], b.med[ss], col = col, lwd = lwd)
  }
  invisible(cbind(b.lo, b.med, b.hi))
}


#### my_termplot for plotting ####

termplot <- function (model, data = NULL, envir = environment(formula(model)), 
          partial.resid = FALSE, rug = FALSE, terms = NULL, se = FALSE, 
          xlabs = NULL, ylabs = NULL, main = NULL, col.term = 2, lwd.term = 1.5, 
          col.se = "orange", lty.se = 2, lwd.se = 1, col.res = "gray", 
          cex = 1, pch = par("pch"), col.smth = "darkred", lty.smth = 2, 
          span.smth = 2/3, ask = dev.interactive() && nb.fig < n.tms, 
          use.factor.levels = TRUE, smooth = NULL, ylim = "common", 
          plot = TRUE, transform.x = FALSE, ...) 
{
  which.terms <- terms
  terms <- if (is.null(terms)) 
    predict(model, type = "terms", se.fit = se)
  else predict(model, type = "terms", se.fit = se, terms = terms)
  n.tms <- ncol(tms <- as.matrix(if (se) 
    terms$fit
    else terms))
  transform.x <- rep_len(transform.x, n.tms)
  mf <- model.frame(model)
  if (is.null(data)) 
    data <- eval(model$call$data, envir)
  if (is.null(data)) 
    data <- mf
  use.rows <- if (NROW(tms) < NROW(data)) 
    match(rownames(tms), rownames(data))
  nmt <- colnames(tms)
  if (any(grepl(":", nmt, fixed = TRUE))) 
    warning("'model' appears to involve interactions: see the help page", 
            domain = NA, immediate. = TRUE)
  cn <- parse(text = nmt, keep.source = FALSE)
  if (!is.null(smooth)) 
    smooth <- match.fun(smooth)
  if (is.null(ylabs)) 
    ylabs <- paste("Partial for", nmt)
  if (is.null(main)) 
    main <- ""
  else if (is.logical(main)) 
    main <- if (main) 
      deparse(model$call, 500)
  else ""
  else if (!is.character(main)) 
    stop("'main' must be TRUE, FALSE, NULL or character (vector).")
  main <- rep_len(main, n.tms)
  pf <- envir
  carrier <- function(term, transform) {
    if (length(term) > 1L) {
      if (transform) 
        tms[, i]
      else carrier(term[[2L]], transform)
    }
    else eval(term, data, enclos = pf)
  }
  carrier.name <- function(term) {
    if (length(term) > 1L) 
      carrier.name(term[[2L]])
    else as.character(term)
  }
  in.mf <- nmt %in% names(mf)
  is.fac <- sapply(nmt, function(i) i %in% names(mf) && is.factor(mf[, 
                                                                     i]))
  if (!plot) {
    outlist <- vector("list", sum(in.mf))
    for (i in 1L:n.tms) {
      if (!in.mf[i]) 
        next
      if (is.fac[i]) {
        xx <- mf[, nmt[i]]
        if (!is.null(use.rows)) 
          xx <- xx[use.rows]
        ww <- match(levels(xx), xx, nomatch = 0L)
      }
      else {
        xx <- carrier(cn[[i]], transform.x[i])
        if (!is.null(use.rows)) 
          xx <- xx[use.rows]
        ww <- match(sort(unique(xx)), xx)
      }
      outlist[[i]] <- if (se) 
        data.frame(x = xx[ww], y = tms[ww, i], se = terms$se.fit[ww, 
                                                                 i], row.names = NULL)
      else data.frame(x = xx[ww], y = tms[ww, i], row.names = NULL)
    }
    attr(outlist, "constant") <- attr(terms, "constant")
    if (se && is.null(attr(outlist, "constant"))) 
      attr(outlist, "constant") <- attr(terms$fit, "constant")
    names(outlist) <- sapply(cn, carrier.name)[in.mf]
    return(outlist)
  }
  if (!is.null(smooth)) 
    smooth <- match.fun(smooth)
  if (is.null(ylabs)) 
    ylabs <- paste("Partial for", nmt)
  if (is.null(main)) 
    main <- ""
  else if (is.logical(main)) 
    main <- if (main) 
      deparse(model$call, 500)
  else ""
  else if (!is.character(main)) 
    stop("'main' must be TRUE, FALSE, NULL or character (vector).")
  main <- rep_len(main, n.tms)
  if (is.null(xlabs)) {
    xlabs <- unlist(lapply(cn, carrier.name))
    if (any(transform.x)) 
      xlabs <- ifelse(transform.x, lapply(cn, deparse), 
                      xlabs)
  }
  if (partial.resid || !is.null(smooth)) {
    pres <- residuals(model, "partial")
    if (!is.null(which.terms)) 
      pres <- pres[, which.terms, drop = FALSE]
  }
  se.lines <- function(x, iy, i, ff = 2) {
    tt <- ff * terms$se.fit[iy, i]
    lines(x, tms[iy, i] + tt, lty = lty.se, lwd = lwd.se, 
          col = col.se)
    lines(x, tms[iy, i] - tt, lty = lty.se, lwd = lwd.se, 
          col = col.se)
  }
  nb.fig <- prod(par("mfcol"))
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  ylims <- ylim
  if (identical(ylims, "common")) {
    ylims <- if (!se) 
      range(tms, na.rm = TRUE)
    else range(tms + 1.05 * 2 * terms$se.fit, tms - 1.05 * 
                 2 * terms$se.fit, na.rm = TRUE)
    if (partial.resid) 
      ylims <- range(ylims, pres, na.rm = TRUE)
    if (rug) 
      ylims[1L] <- ylims[1L] - 0.07 * diff(ylims)
  }
  for (i in 1L:n.tms) {
    if (identical(ylim, "free")) {
      ylims <- range(tms[, i], na.rm = TRUE)
      if (se) 
        ylims <- range(ylims, tms[, i] + 1.05 * 2 * terms$se.fit[, 
                                                                 i], tms[, i] - 1.05 * 2 * terms$se.fit[, i], 
                       na.rm = TRUE)
      if (partial.resid) 
        ylims <- range(ylims, pres[, i], na.rm = TRUE)
      if (rug) 
        ylims[1L] <- ylims[1L] - 0.07 * diff(ylims)
    }
    if (!in.mf[i]) 
      next
    if (is.fac[i]) {
      ff <- mf[, nmt[i]]
      if (!is.null(model$na.action)) 
        ff <- naresid(model$na.action, ff)
      ll <- levels(ff)
      xlims <- range(seq_along(ll)) + c(-0.5, 0.5)
      xx <- as.numeric(ff)
      if (rug) {
        xlims[1L] <- xlims[1L] - 0.07 * diff(xlims)
        xlims[2L] <- xlims[2L] + 0.03 * diff(xlims)
      }
      plot(1, 0, type = "n", xlab = xlabs[i], ylab = ylabs[i], 
           xlim = xlims, ylim = ylims, main = main[i], xaxt = "n", 
           ...)
      if (use.factor.levels) 
        axis(1, at = seq_along(ll), labels = ll, ...)
      else axis(1)
      for (j in seq_along(ll)) {
        ww <- which(ff == ll[j])[c(1, 1)]
        jf <- j + c(-0.4, 0.4)
        lines(jf, tms[ww, i], col = col.term, lwd = lwd.term, 
              ...)
        if (se) 
          se.lines(jf, iy = ww, i = i)
      }
    }
    else {
      xx <- carrier(cn[[i]], transform.x[i])
      if (!is.null(use.rows)) 
        xx <- xx[use.rows]
      xlims <- range(xx, na.rm = TRUE)
      if (rug) 
        xlims[1L] <- xlims[1L] - 0.07 * diff(xlims)
      oo <- order(xx)
      plot(xx[oo], tms[oo, i], type = "l", xlab = xlabs[i], 
           ylab = ylabs[i], xlim = xlims, ylim = ylims, 
           main = main[i], col = col.term, lwd = lwd.term, 
           ...)
      if (se) 
        se.lines(xx[oo], iy = oo, i = i)
    }
    if (partial.resid) {
      if (!is.fac[i] && !is.null(smooth)) {
        smooth(xx, pres[, i], lty = lty.smth, cex = cex, 
               pch = pch, col = col.res, col.smooth = col.smth, 
               span = span.smth)
      }
      else points(xx, pres[, i], cex = cex, pch = pch, 
                  col = col.res)
    }
    if (rug) {
      n <- length(xx)
      lines(rep.int(jitter(xx), rep.int(3, n)), rep.int(ylims[1L] + 
                                                          c(0, 0.05, NA) * diff(ylims), n))
      if (partial.resid) 
        lines(rep.int(xlims[1L] + c(0, 0.05, NA) * diff(xlims), 
                      n), rep.int(pres[, i], rep.int(3, n)))
    }
  }
  invisible(n.tms)
}

