#' Get AnnotationHub object
#'
#' [`AnnotationHub()`][AnnotationHub::AnnotationHub()] restrict the
#' AnnotationHub snapshotDate before the Bioconductor version, this function
#' just bypass the restriction to let you query all AnnotationHub data.
#'
#' @param ... See [`AnnotationHub`][AnnotationHub::AnnotationHub] for additional
#' arguments.
#' @param version A valid Bioconductor version, see:
#' <https://bioconductor.org/config.yaml> for all version.
#'
#' @param hub The URL for the online AnnotationHub.
#' @param cache The file system location of the local AnnotationHub cache.
#' @param proxy Set the proxy.
#' @param local A bool, whether to use the cache.
#' @return A [`AnnotationHub`][AnnotationHub::AnnotationHub] instance.
#' @export
annohub <- function(..., version = NULL, hub = NULL, cache = NULL,
                    proxy = NULL, local = NULL) {
    check_bioc_installed("AnnotationHub", "to use `annohub()`")
    if (is.null(version)) {
        if (is_installed("BiocManager")) {
            version <- getExportedValue("BiocManager", "version")()
        } else {
            cli::cli_abort(paste(
                "Please provide {.arg version} or you should install",
                "{.pkg BiocManager} for us to get a version"
            ))
        }
    }
    url <- hub %||% AnnotationHub::getAnnotationHubOption("URL")
    cache <- cache %||% AnnotationHub::getAnnotationHubOption("CACHE")
    proxy <- proxy %||% AnnotationHub::getAnnotationHubOption("PROXY")
    localHub <- localHub %||% AnnotationHub::getAnnotationHubOption("LOCAL")
    db_path <- annohub_namespace(".create_cache")(
        "AnnotationHub",
        url = url, cache = cache, proxy = proxy, local = local,
        ask = AnnotationHub::getAnnotationHubOption("ASK")
    )
    if (!local) {
        rlang::try_fetch(
            {
                dates <- as.POSIXlt(
                    annohub_namespace(".possibleDates")(db_path),
                    format = "%Y-%m-%d"
                )
                restrict <- as.POSIXlt(
                    annohub_namespace(".biocVersionDate")(version),
                    format = "%Y-%m-%d"
                )
                if (length(restrict)) {
                    db_date <- as.character(max(dates[dates <= restrict]))
                } else {
                    db_date <- as.character(max(dates))
                }
            },
            error = function(err) {
                stop(
                    "failed to connect", "\n  reason: ", conditionMessage(err),
                    "\n  Consider rerunning with 'localHub=TRUE'"
                )
            }
        )
    } else {
        dates <- annohub_namespace(".possibleDates")(db_path)
        db_date <- dates[length(dates)]
    }
    db_uid <- annohub_namespace(".db_uid0")(db_path, db_date, local)
    hub <- methods::new("AnnotationHub",
        cache = cache, hub = url, date = db_date,
        .db_path = db_path, .db_uid = db_uid, isLocalHub = local,
        ...
    )
    cli::cli_inform("snapshotDate(): ", AnnotationHub::snapshotDate(hub))
    if (!local) {
        index <- annohub_namespace(".db_create_index")(hub)
        hub <- annohub_namespace(".db_index<-")(hub, value = index)
    } else {
        index <- annohub_namespace(".db_index_file")(hub)
        hub <- annohub_namespace(".db_index<-")(hub, value = index)
        hub <- annohub_namespace(".subsethub")(hub)
    }
    hub
}

annohub_namespace <- function(fn) from_namespace("AnnotationHub", fn)
