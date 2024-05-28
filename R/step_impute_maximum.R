#' @importFrom purrr map_dbl
maximums <- function(x, wts = NULL) {
  if (NCOL(x) == 0) {
    return(vapply(x, max, c(max = 0), na.rm = TRUE))
  }
  res <- apply(x, 2, max, na.rm = TRUE)
  res
}


#' Impute numeric data using the maximum (or factor * maximum)
#'
#' `step_impute_maximum()` creates a *specification* of a recipe step that will
#' substitute missing values of numeric variables by the training set maximum of
#' those variables.
#'
#' @param recipe A recipe object. The step will be added to the sequence
#' of operations for this recipe.
#' @param ... One or more selector functions to choose variables for this step.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing have
#' been estimated.
#' @param maximums A named numeric vector of maximums. This is `NULL` until
#'  computed by [prep()]. Note that, if the original data are integers,
#'  the maximum will be converted to an integer to maintain the same data type.
#' @param factor A numeric value greater than or equal to 1 to scale the maximum
#' values. To replace values with the double-maximum, set this argument to 2.
#' Default value of 1.
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
#' @details `step_impute_maximum` estimates the variable maximums from the data
#'  used in the `training` argument of `prep.recipe`. `bake.recipe` then applies
#'  the new values to new data sets using these maximums.
#'
#'
#' # Tidying
#'
#' When you [`tidy()`][tidy.recipe()] this step, a tibble with
#' columns `terms` (the selectors or variables selected) and `model`
#' (the maximum value) is returned.
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
#' missing_examples <- c(933, 1042, 1671, 1796)
#' library(recipes)
#' rec <- recipe(Price ~ ., data = credit_tr)
#'
#' impute_rec <- rec %>%
#'   step_impute_maximum(Income, Assets, Debt)
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
step_impute_maximum <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           maximums = NULL,
           factor = 1,
           skip = FALSE,
           id = recipes::rand_id("impute_maximum")) {
    add_step(
      recipe,
      step_impute_maximum_new(
        terms = rlang::enquos(...),
        role = role,
        trained = trained,
        maximums = maximums,
        factor = factor,
        skip = skip,
        id = id,
        case_weights = NULL
      )
    )
  }

#' @importFrom recipes step
step_impute_maximum_new <-
  function(terms, role, trained, maximums, factor, skip, id, case_weights) {
    step(
      subclass = "impute_maximum",
      terms = terms,
      role = role,
      trained = trained,
      maximums = maximums,
      factor = factor,
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
prep.step_impute_maximum <- function(x, training, info = NULL, ...) {
  col_names <- recipes_eval_select(x$terms, training, info)
  check_type(training[, col_names], types = c("double", "integer"))

  # wts <- get_case_weights(info, training)
  # were_weights_used <- are_weights_used(wts, unsupervised = TRUE)
  # if (isFALSE(were_weights_used)) {
  #   wts <- NULL
  # }
  were_weights_used <- FALSE
  wts <- NULL
  if (x$factor < 1) {
    rlang::warn(paste0("Scaling `factor` should take a value greater ",
                       "than or equal to 1. Setting to 1."))
    x$factor <- 1
  }

  maximums <- maximums(training[, col_names], wts = wts)
  maximums <- maximums * x$factor

  maximums <- purrr::map2(maximums, training[, col_names], cast2)

  step_impute_maximum_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    maximums = maximums,
    factor = x$factor,
    skip = x$skip,
    id = x$id,
    case_weights = were_weights_used
  )
}

#' @importFrom recipes check_new_data bake
#' @importFrom vctrs vec_cast
#' @export
bake.step_impute_maximum <- function(object, new_data, ...) {
  check_new_data(names(object$maximums), object, new_data)

  for (i in names(object$maximums)) {
    if (any(is.na(new_data[[i]]))) {
      new_data[[i]] <- vctrs::vec_cast(new_data[[i]], object$maximums[[i]])
    }
    new_data[is.na(new_data[[i]]), i] <- object$maximums[[i]]
  }
  new_data
}

#' @importFrom recipes print_step
#' @export
print.step_impute_maximum <-
  function(x, width = max(20, options()$width - 30), ...) {
    title <- "maximum imputation for "
    print_step(names(x$maximums), x$terms, x$trained, title, width,
               case_weights = x$case_weights)
    invisible(x)
  }

#' @rdname step_impute_maximum
#' @param x A `step_impute_maximum` object.
#' @importFrom recipes sel2char is_trained
#' @importFrom tibble tibble
#' @importFrom rlang na_dbl
#' @importFrom generics tidy
#' @importFrom vctrs list_unchop
#' @export
tidy.step_impute_maximum <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble::tibble(
      terms = names(x$maximums),
      value = vctrs::list_unchop(unname(x$maximums), ptype = double())
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
required_pkgs.step_impute_maximum <- function(x, ...) {
  c("scaleParetoImputeMin")
}
