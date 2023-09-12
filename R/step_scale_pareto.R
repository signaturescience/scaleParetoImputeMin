# Consider replacing with https://github.com/CVUA-RRW/tidySpectR


#' Scaling Numeric Data by the square root of the standard deviation
#'
#' `step_scale_pareto()` creates a *specification* of a recipe step that will
#' normalize numeric data to have a square root of the standard deviation of
#' one.
#'
#' @param recipe A recipe object. The step will be added to the sequence
#' of operations for this recipe.
#' @param ... One or more selector functions to choose variables for this step.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing have
#' been estimated.
#' @param sds A named numeric vector of standard deviations. This is `NULL`
#'  until computed by [prep()].
#' @param na_rm A logical value indicating whether `NA` values should be removed
#' when computing the standard deviation.
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
#' @importFrom stats sd
#' @details Scaling data means that the standard deviation of a
#'  variable is divided out of the data. `step_scale_pareto` estimates
#'  the variable square root standard deviations from the data used in the
#'  `training` argument of `prep.recipe`.
#'  `bake.recipe` then applies the scaling to new data sets
#'  using these square root standard deviations.
#'
#'  # Tidying
#'
#'  When you [`tidy()`][tidy.recipe()] this step, a tibble with columns
#'  `terms` (the selectors or variables selected) and `value` (the
#'  square root standard deviations) is returned.
#'
#'
#' @examplesIf rlang::is_installed("modeldata")
#' data(biomass, package = "modeldata")
#'
#' biomass_tr <- biomass[biomass$dataset == "Training", ]
#' biomass_te <- biomass[biomass$dataset == "Testing", ]
#'
#' rec <- recipes::recipe(
#'   HHV ~ carbon + hydrogen + oxygen + nitrogen + sulfur,
#'   data = biomass_tr
#' )
#'
#' scaled_trans <- rec |>
#'   step_scale_pareto(carbon, hydrogen)
#'
#' scaled_obj <- recipes::prep(scaled_trans, training = biomass_tr)
#'
#' transformed_te <- recipes::bake(scaled_obj, biomass_te)
#'
#' biomass_te[1:10, names(transformed_te)]
#' transformed_te
#' generics::tidy(scaled_trans, number = 1)
#' generics::tidy(scaled_obj, number = 1)
step_scale_pareto <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           sds = NULL,
           na_rm = TRUE,
           skip = FALSE,
           id = recipes::rand_id("pareto")) {

    recipes_pkg_check(required_pkgs.step_scale_pareto())

    # terms <- recipes::ellipse_check(...)

    add_step(
      recipe,
      step_scale_pareto_new(
        terms = rlang::enquos(...),
        role = role,
        trained = trained,
        sds = sds,
        na_rm = na_rm,
        skip = skip,
        id = id
      )
    )
  }

#' @importFrom recipes step
step_scale_pareto_new <-
  function(terms, role, trained, sds, na_rm, skip, id) {
    step(
      subclass = "scale_pareto",
      terms = terms,
      role = role,
      trained = trained,
      sds = sds,
      na_rm = na_rm,
      skip = skip,
      id = id
    )
  }

#' @importFrom glue glue glue_collapse
#' @importFrom rlang warn
sd_check2 <- function (x) {
  zero_sd <- which(x < .Machine$double.eps)
  if (length(zero_sd) > 0) {
    glue_cols <- glue::glue_collapse(glue::glue("`{names(zero_sd)}`"),
                                     sep = ", ", last = " and ")
    rlang::warn(glue::glue("Column(s) have zero variance so scaling cannot be used: {glue_cols}. ",
                           "Consider using `step_zv()` to remove those columns before normalizing"))
    x[zero_sd] <- 1
  }
  x
}

#' @export
#' @importFrom recipes recipes_eval_select
#' @importFrom recipes variances check_type prep
prep.step_scale_pareto <- function(x, training, info = NULL, ...) {
  col_names <- recipes_eval_select(x$terms, training, info)
  check_type(training[, col_names], types = c("double", "integer"))

  vars <- variances(training[, col_names], wts = NULL, na_rm = x$na_rm)
  sds <- sqrt(vars)
  sds <- sd_check2(sds)
  sds <- sqrt(sds)

  step_scale_pareto_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    sds = sds,
    na_rm = x$na_rm,
    skip = x$skip,
    id = x$id
  )
}

#' @export
#' @importFrom recipes check_new_data bake
bake.step_scale_pareto <- function(object, new_data, ...) {
  check_new_data(names(object$sds), object, new_data)

  for (column in names(object$sds)) {
    sd <- object$sds[column]
    new_data[[column]] <- new_data[[column]] / sd
  }
  new_data
}

#' @importFrom recipes print_step
#' @export
print.step_scale_pareto <- function(x, width = max(20, options()$width - 30),
                                    ...) {
  title <- "Pareto scaling for "
  print_step(
    tr_obj = names(x$sds),
    untr_obj = x$terms,
    x$trained,
    title = title,
    width = width)
  invisible(x)
}


#' @rdname step_scale_pareto
#' @param x A `step_scale_pareto` object.
#' @importFrom recipes sel2char is_trained
#' @importFrom tibble tibble as_tibble
#' @importFrom rlang na_dbl
#' @importFrom generics tidy
#' @export
tidy.step_scale_pareto <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble::tibble(
      terms = names(x$sds),
      value = unname(x$sds)
    )
  } else {
    term_names <- sel2char(x$terms)
    res <- tibble::tibble(
      terms = term_names,
      value = rlang::na_dbl
    )
  }
  res$id <- x$id
  res
}

#' @export
#' @importFrom generics required_pkgs
required_pkgs.step_scale_pareto <- function(x, ...) {
  c("scaleParetoImputeMin")
}
