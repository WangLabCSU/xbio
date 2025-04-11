prop_match <- function(choices, ..., default = NULL) {
    new_property(class_character,
        validator = function(value) {
            if (any(value == choices)) {
                return(NULL)
            }
            paste("must be one of", oxford_or(choices))
        },
        ...,
        default = default %||% choices[1L]
    )
}

prop_number_decimal <- function(getter = NULL, setter = NULL, ...,
                                min = NULL, max = NULL,
                                allow_infinite = TRUE, allow_na = FALSE) {
    force(min)
    force(max)
    force(allow_infinite)
    force(allow_na)
    new_property(
        class_double,
        getter = getter,
        setter = setter,
        validator = function(value) {
            if (is.call(min)) min <- eval(min)
            if (is.call(max)) max <- eval(max)
            if (0 == (exit_code <- .rlang_check_number(
                value,
                allow_decimal = TRUE,
                min = min,
                max = max,
                allow_infinite = allow_infinite,
                allow_na = allow_na,
                allow_null = FALSE
            ))) {
                return(NULL)
            }
            sprintf(
                "must be %s, not %s",
                .rlang_allow_number(
                    value,
                    exit_code,
                    allow_decimal = TRUE,
                    min = min, max = max,
                    allow_na = allow_na,
                    allow_null = FALSE
                ),
                obj_type_friendly(value)
            )
        },
        ...
    )
}

#' @importFrom S7 new_property class_double
prop_number_whole <- function(getter = NULL, setter = NULL, ...,
                              min = NULL, max = NULL, allow_na = FALSE) {
    force(min)
    force(max)
    force(allow_na)
    new_property(
        class_integer,
        getter = getter,
        setter = setter,
        validator = function(value) {
            if (is.call(min)) min <- eval(min)
            if (is.call(max)) max <- eval(max)
            if (0 == (exit_code <- .rlang_check_number(
                value,
                allow_decimal = FALSE,
                min = min,
                max = max,
                allow_infinite = FALSE,
                allow_na = allow_na,
                allow_null = FALSE
            ))) {
                return(NULL)
            }
            sprintf(
                "must be %s, not %s",
                .rlang_allow_number(
                    value,
                    exit_code,
                    allow_decimal = FALSE,
                    min = min, max = max,
                    allow_na = allow_na,
                    allow_null = FALSE
                ),
                obj_type_friendly(value)
            )
        },
        ...
    )
}
