cast2 <- function (x, ref) {
  if (is.factor(ref)) {
    x <- factor(x, levels = levels(ref), ordered = is.ordered(ref))
  }
  else {
    if (is.integer(ref) & !is.factor(ref)) {
      x <- as.integer(round(x, 0))
    }
  }
  x
}

weighted_median_impl2 <- function (x, wts)  {
  order_x <- order(x)
  x <- x[order_x]
  wts <- wts[order_x]
  wts_norm <- cumsum(wts)/sum(wts)
  ps <- min(which(wts_norm > 0.5))
  x[ps]
}

cov2pca2 <- function (cv_mat) {
  res <- eigen(cv_mat)
  list(sdev = sqrt(res$values), rotation = res$vectors)
}

#' @importFrom rlang arg_match
#' @importFrom stats complete.cases cov.wt
wt_calcs2 <- function (x, wts, statistic = "mean") {
  statistic <- rlang::arg_match(statistic, c("mean", "var",
                                             "cor", "cov", "pca", "median"))
  if (!is.data.frame(x)) {
    x <- data.frame(x)
  }
  if (is.null(wts)) {
    wts <- rep(1L, nrow(x))
  }
  complete <- stats::complete.cases(x) & !is.na(wts)
  wts <- wts[complete]
  x <- x[complete, , drop = FALSE]
  res <- stats::cov.wt(x, wt = wts, cor = statistic == "cor")
  if (statistic == "mean") {
    res <- unname(res[["center"]])
  }
  else if (statistic == "median") {
    res <- weighted_median_impl2(x$x, wts)
  }
  else if (statistic == "var") {
    res <- unname(diag(res[["cov"]]))
  }
  else if (statistic == "pca") {
    res <- cov2pca2(res$cov)
  }
  else if (statistic == "cov") {
    res <- res[["cov"]]
  }
  else {
    res <- res[["cor"]]
  }
  res
}

#' @importFrom purrr map_dbl
minimums <- function(x, wts = NULL) {
  if (NCOL(x) == 0) {
    return(vapply(x, min, c(min = 0), na.rm = TRUE))
  }
  if (is.null(wts)) {
    res <- apply(x, 2, min, na.rm = TRUE)
  } else {
    wts <- as.double(wts)
    res <- purrr::map_dbl(x, ~ wt_calcs2(.x, wts, statistic = "min"))
  }
  res
}


#' Impute numeric data using the minimum (or half-minimum)
#'
#' `step_impute_minimum()` creates a *specification* of a recipe step that will
#' substitute missing values of numeric variables by the training set minimum of
#' those variables.
#'
#' @param recipe A recipe object. The step will be added to the sequence
#' of operations for this recipe.
#' @param ... One or more selector functions to choose variables for this step.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing have
#' been estimated.
#' @param minimums A named numeric vector of minimums. This is `NULL` until
#'  computed by [prep()]. Note that, if the original data are integers,
#'  the minimum will be converted to an integer to maintain the same data type.
#' @param factor A numeric value greater than 0 and less than or equal to 1 to
#' scale the minimum values. To replace values with the half-minimum, set this
#' argument to 0.5. Default value of 1.
#' @param rand_value Logical value (`TRUE` or `FALSE` (default)) as to whether
#' a random uniform value greater than 0 and less than the minimum should be
#' used for NA replacement. Note that the random value is set for all `NA`
#' values for each numeric column being imputed, i.e., each `NA` does not get
#' its own random value. If set to FALSE, the minimum value is used.
#' @param skip A logical. Should the step be skipped when the recipe is baked by
#' bake()? While all operations are baked when prep() is run, some operations
#' may not be able to be conducted on new data (e.g. processing the outcome
#' variable(s)). Care should be taken when using skip = TRUE as it may affect
#' the computations for subsequent operations.
#' @param id A character string that is unique to this step to identify it.
#' @export
#' @importFrom recipes add_step rand_id prep bake recipes_pkg_check step
#' @importFrom recipes printer check_type is_trained sel2char ellipse_check
#' @importFrom generics tidy required_pkgs
#' @importFrom rlang enquos
#' @details `step_impute_minimum` estimates the variable minimums from the data
#'  used in the `training` argument of `prep.recipe`. `bake.recipe` then applies
#'  the new values to new data sets using these minimums.
#'
#'
#' # Tidying
#'
#' When you [`tidy()`][tidy.recipe()] this step, a tibble with
#' columns `terms` (the selectors or variables selected) and `model`
#' (the minimum value) is returned.
#'
#' @examplesIf rlang::is_installed("modeldata")
#' data("credit_data", package = "modeldata")
#'
#' ## missing data per column
#' vapply(credit_data, function(x) mean(is.na(x)), c(num = 0))
#'
#' set.seed(342)
#' in_training <- sample(1:nrow(credit_data), 2000)
#'
#' credit_tr <- credit_data[in_training, ]
#' credit_te <- credit_data[-in_training, ]
#' missing_examples <- c(14, 394, 565)
#' library(recipes)
#' rec <- recipe(Price ~ ., data = credit_tr)
#'
#' impute_rec <- rec %>%
#'   step_impute_minimum(Income, Assets, Debt)
#'
#' imp_models <- prep(impute_rec, training = credit_tr)
#'
#' imputed_te <- bake(imp_models, new_data = credit_te, everything())
#'
#' credit_te[missing_examples, ]
#' imputed_te[missing_examples, names(credit_te)]
#'
#' tidy(impute_rec, number = 1)
#' tidy(imp_models, number = 1)
step_impute_minimum <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           minimums = NULL,
           factor = 1,
           rand_value = FALSE,
           skip = FALSE,
           id = recipes::rand_id("impute_minimum")) {
    add_step(
      recipe,
      step_impute_minimum_new(
        terms = rlang::enquos(...),
        role = role,
        trained = trained,
        minimums = minimums,
        factor = factor,
        rand_value = rand_value,
        skip = skip,
        id = id,
        case_weights = NULL
      )
    )
  }

#' @importFrom recipes step
step_impute_minimum_new <-
  function(terms, role, trained, minimums, factor, rand_value, skip, id, case_weights) {
    step(
      subclass = "impute_minimum",
      terms = terms,
      role = role,
      trained = trained,
      minimums = minimums,
      factor = factor,
      rand_value = rand_value,
      skip = skip,
      id = id,
      case_weights = case_weights
    )
  }

#' @importFrom purrr map2
#' @importFrom recipes recipes_eval_select get_case_weights are_weights_used
#' @importFrom recipes variances check_type prep
#' @importFrom rlang warn
#' @importFrom stats runif
#' @export
prep.step_impute_minimum <- function(x, training, info = NULL, ...) {
  col_names <- recipes_eval_select(x$terms, training, info)
  check_type(training[, col_names], types = c("double", "integer"))

  wts <- get_case_weights(info, training)
  were_weights_used <- are_weights_used(wts, unsupervised = TRUE)
  if (isFALSE(were_weights_used)) {
    wts <- NULL
  }

  if (x$factor <= 0 | x$factor > 1) {
    rlang::warn(paste0("Scaling `factor` should take either a value greater ",
                       "than 0 and less than 1"))
  }

  minimums <- minimums(training[, col_names], wts = wts)
  minimums <- minimums * x$factor
  if (!is.logical(x$rand_value)) {
    rlang::warn("`rand_value` must be FALSE or TRUE. Setting to FALSE")
    x$rand_value <- FALSE
  }
  if (x$rand_value == TRUE) {
    nm <- names(minimums)
    minimums <- stats::runif(length(minimums),
                             min = rep(0, length(minimums)),
                             max = minimums)
    names(minimums) <- nm
  }

  minimums <- purrr::map2(minimums, training[, col_names], cast2)

  step_impute_minimum_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    minimums = minimums,
    factor = x$factor,
    rand_value = x$rand_value,
    skip = x$skip,
    id = x$id,
    case_weights = were_weights_used
  )
}

#' @importFrom recipes check_new_data bake
#' @importFrom vctrs vec_cast
#' @export
bake.step_impute_minimum <- function(object, new_data, ...) {
  check_new_data(names(object$minimums), object, new_data)

  for (i in names(object$minimums)) {
    if (any(is.na(new_data[[i]]))) {
      new_data[[i]] <- vctrs::vec_cast(new_data[[i]], object$minimums[[i]])
    }
    new_data[is.na(new_data[[i]]), i] <- object$minimums[[i]]
  }
  new_data
}

#' @importFrom recipes print_step
#' @export
print.step_impute_minimum <-
  function(x, width = max(20, options()$width - 30), ...) {
    title <- "minimum imputation for "
    print_step(names(x$minimums), x$terms, x$trained, title, width,
               case_weights = x$case_weights)
    invisible(x)
  }

#' @rdname step_impute_minimum
#' @param x A `step_impute_minimum` object.
#' @importFrom recipes sel2char is_trained
#' @importFrom tibble tibble
#' @importFrom rlang na_dbl
#' @importFrom generics tidy
#' @importFrom vctrs list_unchop
#' @export
tidy.step_impute_minimum <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble::tibble(
      terms = names(x$minimums),
      value = vctrs::list_unchop(unname(x$minimums), ptype = double())
    )
  } else {
    term_names <- sel2char(x$terms)
    res <- tibble::tibble(terms = term_names, value = rlang::na_dbl)
  }
  res$id <- x$id
  res
}

#' @export
#' @importFrom generics required_pkgs
required_pkgs.step_impute_minimum <- function(x, ...) {
  c("scalePareto")
}
