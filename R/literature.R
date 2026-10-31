# =============================================================================
# R/literature.R - optional literature mapping of prioritized genes
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
# Purpose : Count PubMed/Europe PMC records per prioritized gene and context
#           (needs internet). Optional; not used for any manuscript result.
# Used by : run_all.R when RUN_LITERATURE = TRUE
# =============================================================================

# Gene set for the literature mapping. 00_run_analysis.R sets LITERATURE_GENES to
# the network's prioritized nodes; otherwise they are read from
# AMR_reanalysis_outputs_v5*/network/prioritized_central_args.csv in getwd().
# Search terms: the FIRST term is the one queried (aliases are kept for audit only,
# as in v4). Add a row here for any newly prioritized gene.
LITERATURE_SEARCH_TERMS <- list(
  "ant(6)-Ia"   = c("ANT(6)-Ia", "aadE"),
  "ant(3'')-Ia" = c("ANT(3'')-Ia", "aadA1"),
  "erm(F)"      = c("ermF", "erm(F)"),
  "tet(W)"      = c("tet(W)", "tetW"),
  "tet(Q)"      = c("tet(Q)", "tetQ"),
  "tet(M)"      = c("tet(M)", "tetM"),
  "tet(O)"      = c("tet(O)", "tetO"),
  "tet(32)"     = c("tet(32)", "tet32"),
  "tet(40)"     = c("tet(40)", "tet40"),
  "msr(D)"      = c("msr(D)", "msrD"),
  "mef(A)"      = c("mef(A)", "mefA"),
  "aph(3')-III" = c("APH(3')-III", "aphA-3"),
  "aph(3'')-Ib" = c("APH(3'')-Ib", "strA"),
  "aph(6)-Id"   = c("APH(6)-Id", "strB"),
  # v4 labels (kept so the earlier cache/search table can still be read)
  "ANT(6)-Ia"   = c("ANT(6)-Ia", "aadE"),
  "ANT(3'')-Ia" = c("ANT(3'')-Ia", "aadA1"),
  "ErmF"        = c("ermF", "erm(F)"),
  "APH(3')-Ib"  = c("APH(3')-Ib")
)

literature_genes <- function() {
  if (exists("LITERATURE_GENES", envir = globalenv()) &&
      length(get("LITERATURE_GENES", envir = globalenv()))) {
    return(get("LITERATURE_GENES", envir = globalenv()))
  }
  candidates <- c(file.path("network", "prioritized_central_args.csv"),
                  Sys.glob(file.path(c("output", "AMR_reanalysis_outputs_v5*"), "network", "prioritized_central_args.csv")))
  candidates <- candidates[file.exists(candidates)]
  if (!length(candidates)) stop("Set LITERATURE_GENES or run the network step first.")
  utils::read.csv(candidates[1], check.names = FALSE, stringsAsFactors = FALSE)$node
}

literature_search_terms <- function(genes) {
  missing <- setdiff(genes, names(LITERATURE_SEARCH_TERMS))
  if (length(missing)) {
    stop("Add search terms to LITERATURE_SEARCH_TERMS for: ", paste(missing, collapse = ", "))
  }
  LITERATURE_SEARCH_TERMS[genes]
}

literature_palette <- function(genes) {
  okabe_ito <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
                 "#56B4E9", "#F0E442", "#000000")
  if (length(genes) > length(okabe_ito)) stop("Palette supports up to 8 genes.")
  stats::setNames(okabe_ito[seq_along(genes)], genes)
}

# run_literature_search keeps the original search or plotting rules.
run_literature_search <- function() {
  # Research in Veterinary Science | QUERY STRATUM x ARG | retrieved records
  # This figure summarizes unique bibliographic search results, NOT verified
  # Query strata are keywords only. No manual review is implied.
  # The five gene labels are provisional v4 network names pending annotation audit.
  # Search result cache is stored separately for this exact 10-year window.
  # Existing search tables are never overwritten unless RUN_SEARCH is TRUE.
  # R >= 4.2 recommended. All outputs contain actual retrieved counts only.

  CFG <- list(
    OUTDIR = "ARG_literature_search_2016_2026_exact",
    # Fresh output folder prevents mixing earlier (unrestricted-date) searches.
    RUN_SEARCH = NA,   # NA = search when the hit table is missing or lacks any requested gene
    SEARCH_PUBMED = TRUE,
    SEARCH_EUROPEPMC = TRUE,
    IMPORT_PREVIOUS_CANDIDATES = FALSE,    # do NOT merge unharmonized earlier 4-db searches
    PREVIOUS_CANDIDATES = "ARG_OneHealth_literature_RVS/literature_tables/screening_candidates.csv",
    USE_ALIASES_FOR_SCREENING = FALSE,    # alias matches NEVER count as gene identity
    MAX_RESULTS_PER_QUERY = 1500L,        # truncated retrieval is flagged, not interpreted as complete
    PAGE_SIZE = 100L,
    HTTP_PAUSE_SECONDS = 0.38,
    EMAIL = Sys.getenv("NCBI_EMAIL", unset = ""),
    SEED = 2026L,
    WIDTH_IN = 12.8,
    HEIGHT_IN = 6.7,
    DPI = 300L,
    MIN_YEAR = 2016L,
    MAX_YEAR = 2026L,
    DATE_START = "2016-10-15",  # inclusive first publication date
    DATE_END   = "2026-10-15",  # inclusive last publication date
    FIGURE_DIR = "ARG_literature_search_2016_2026_exact"
  )

  GENES <- literature_genes()
  SECTORS <- c("Human", "Avian", "Other animals", "Food", "Wastewater",
               "Other environment")
  hits_file <- file.path(CFG$OUTDIR, "tables", "bibliographic_search_hits.csv")
  if (is.na(CFG$RUN_SEARCH)) {
    cached_genes <- if (file.exists(hits_file)) unique(utils::read.csv(hits_file)$gene) else character()
    CFG$RUN_SEARCH <- !all(GENES %in% cached_genes)
    if (CFG$RUN_SEARCH) message("Searching because these genes are not in the hit table: ",
                                paste(setdiff(GENES, cached_genes), collapse = ", "))
  }

  # Search aliases are retrieval aids only. They are not automatic exact-allele
  # equivalents. APH(3')-Ib and APH(3'')-Ib MUST NOT be silently combined.
  SEARCH_TERMS <- literature_search_terms(GENES)

  # Search-stratum keywords are for RETRIEVAL only; not evidence of sample origin.
  QUERIES <- list(
    Human = c("human", "clinical", "patient", "hospital", "community"),
    Avian = c("poultry", "chicken", "avian", "bird"),
    `Other animals` = c("veterinary", "livestock", "pig", "swine", "cattle",
                        "companion animal", "fish", "wildlife"),
    Food = c("food", "meat", "milk", "seafood", "vegetable"),
    Wastewater = c("wastewater", "sewage", "wastewater treatment"),
    `Other environment` = c("compost", "manure", "soil", "river", "sediment"),
    Unrestricted = character(0)
  )

  # NOTE: the source label "feces" alone CANNOT be assigned to Human or Animal.
  # Classify samples by documented provenance, not the word in an abstract.

  needed <- c("httr2", "jsonlite", "xml2", "digest", "dplyr", "tidyr",
              "readr", "ggplot2", "ragg", "svglite")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))

  D <- list(
    root = CFG$OUTDIR,
    tables = file.path(CFG$OUTDIR, "tables"),
    cache = file.path(CFG$OUTDIR, "api_cache"),
    fig = file.path(CFG$OUTDIR, "figures")
  )
  invisible(lapply(D, function(path) dir.create(path, recursive = TRUE,
                                                showWarnings = FALSE)))
  FILES <- list(
    raw = file.path(D$tables, "bibliographic_search_hits.csv"),
    queries = file.path(D$tables, "search_query_log.csv"),
    candidates = file.path(D$tables, "candidate_papers_to_screen.csv"),
    review = file.path(D$tables, "REVIEW_source_ARG_studies_EDIT_ME.csv"),
    screening = file.path(D$tables, "SCREENING_completion_EDIT_ME.csv"),
    identity = file.path(D$tables, "GENE_identity_audit_EDIT_ME.csv"),
    checklist = file.path(D$tables, "PUBLICATION_checklist_EDIT_ME.csv"),
    excluded = file.path(D$tables, "publication_gate_issues.csv")
  )

  nonblank <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
  yes <- function(x) tolower(trimws(as.character(x))) %in% c("yes", "y", "true", "1")
  clean <- function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    trimws(gsub("[[:space:]]+", " ", gsub("<[^>]*>", " ", x)))
  }
  ch <- function(x) {
    if (is.null(x) || length(x) == 0) return("")
    z <- as.character(unlist(x, recursive = TRUE, use.names = FALSE))
    z <- z[!is.na(z) & nzchar(z)]
    if (length(z)) z[1] else ""
  }
  doi_clean <- function(x) {
    x <- tolower(clean(x))
    x <- sub("^https?://(dx\\.)?doi\\.org/", "", x)
    x <- sub("^doi:[[:space:]]*", "", x)
    x
  }
  paper_id <- function(doi, pmid, title, year) {
    if (nonblank(doi)) return(paste0("doi:", doi_clean(doi)))
    if (nonblank(pmid)) return(paste0("pmid:", trimws(pmid)))
    paste0("title:", digest::digest(paste(tolower(clean(title)), year)))
  }
  save_if_missing <- function(x, path) {
    if (!file.exists(path)) {
      readr::write_csv(x, path, na = "")
      message("NEW review file (will never overwrite): ", path)
    } else message("Preserving edited file: ", path)
  }

  fetch_string <- function(url, params, ext = "json") {
    key <- digest::digest(list(url, params), algo = "sha1")
    filename <- file.path(D$cache, paste0(key, ".", ext))
    if (file.exists(filename))
      return(paste(readLines(filename, warn = FALSE, encoding = "UTF-8"),
                   collapse = "\n"))
    req <- httr2::request(url)
    req <- do.call(httr2::req_url_query, c(list(req), params))
    req <- httr2::req_user_agent(req,
       "ARG-evidence-matrix-RVS/1.0 (academic research; contact via NCBI_EMAIL)")
    req <- httr2::req_timeout(req, seconds = 65)
    req <- httr2::req_retry(req, max_tries = 3)
    Sys.sleep(CFG$HTTP_PAUSE_SECONDS)
    response <- httr2::req_perform(req)
    body <- httr2::resp_body_string(response)
    writeLines(body, filename, useBytes = TRUE)
    body
  }
  fetch_json <- function(url, params) {
    jsonlite::fromJSON(fetch_string(url, params, "json"),
                       simplifyVector = FALSE)
  }

  terms_boolean <- function(gene) {
    x <- SEARCH_TERMS[[gene]]
    if (!CFG$USE_ALIASES_FOR_SCREENING) x <- x[1]
    paste0("(", paste(sprintf('"%s"', x), collapse = " OR "), ")")
  }
  build_query <- function(gene, stratum) {
    gn <- terms_boolean(gene)
    context <- QUERIES[[stratum]]
    if (!length(context)) return(gn)
    paste0(gn, " AND (", paste(sprintf('"%s"', context), collapse = " OR "), ")")
  }

  # Apply the exact time window in the database query (not only after retrieval).
  build_pubmed_query <- function(gene, stratum) {
    start <- gsub("-", "/", CFG$DATE_START, fixed = TRUE)
    end <- gsub("-", "/", CFG$DATE_END, fixed = TRUE)
    paste0("(", build_query(gene, stratum), ") AND (\"", start,
           "\"[Date - Publication] : \"", end, "\"[Date - Publication])")
  }
  build_europepmc_query <- function(gene, stratum) {
    paste0("(", build_query(gene, stratum), ") AND FIRST_PDATE:[",
           CFG$DATE_START, " TO ", CFG$DATE_END, "]")
  }
  stopifnot(length(CFG$MIN_YEAR) == 1L, length(CFG$MAX_YEAR) == 1L,
            length(CFG$DATE_START) == 1L, length(CFG$DATE_END) == 1L,
            !is.na(as.Date(CFG$DATE_START)), !is.na(as.Date(CFG$DATE_END)),
            as.Date(CFG$DATE_START) <= as.Date(CFG$DATE_END))
  empty_hit <- function() data.frame(
    gene = character(), query_stratum = character(), database = character(),
    title = character(), abstract = character(), year = character(),
    doi = character(), pmid = character(), journal = character(), url = character(),
    stringsAsFactors = FALSE
  )
  new_hit <- function(gene, stratum, db, title = "", abstract = "", year = "",
                      doi = "", pmid = "", journal = "", url = "") {
    data.frame(gene = gene, query_stratum = stratum, database = db,
               title = clean(ch(title)), abstract = clean(ch(abstract)),
               year = ch(year), doi = doi_clean(ch(doi)), pmid = ch(pmid),
               journal = clean(ch(journal)), url = ch(url),
               stringsAsFactors = FALSE)
  }

  search_pubmed <- function(gene, stratum) {
    query <- build_pubmed_query(gene, stratum)
    url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    fetch_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params <- list(db = "pubmed", term = query, retmode = "json", retmax = 0,
                   tool = "ARG_RVS_review")
    if (nzchar(CFG$EMAIL)) params$email <- CFG$EMAIL
    js <- fetch_json(url, params)
    total <- suppressWarnings(as.integer(ch(js$esearchresult$count)))
    if (!is.finite(total)) stop("PubMed did not provide a valid search count")
    limit <- min(total, as.integer(CFG$MAX_RESULTS_PER_QUERY))
    # PubMed ESearch limits a single query to the first 10,000 identifiers.
    limit <- min(limit, 10000L)
    if (limit == 0L)
      return(list(hits = empty_hit(), total = total, truncated = FALSE))
    ids <- character()
    for (start in seq.int(0L, limit - 1L, by = CFG$PAGE_SIZE)) {
      par2 <- params
      par2$retmax <- min(CFG$PAGE_SIZE, limit - start)
      par2$retstart <- start
      page <- fetch_json(url, par2)
      ids <- c(ids, unlist(page$esearchresult$idlist, use.names = FALSE))
    }
    ids <- unique(as.character(ids))
    if (!length(ids)) stop("PubMed returned no IDs despite positive count")
    hits <- list()
    ix <- 0L
    for (start in seq.int(1L, length(ids), by = CFG$PAGE_SIZE)) {
      selected <- ids[start:min(start + CFG$PAGE_SIZE - 1L, length(ids))]
      pars <- list(db = "pubmed", id = paste(selected, collapse = ","),
                   retmode = "xml", tool = "ARG_RVS_review")
      if (nzchar(CFG$EMAIL)) pars$email <- CFG$EMAIL
      xml <- xml2::read_xml(fetch_string(fetch_url, pars, "xml"))
      articles <- xml2::xml_find_all(xml, ".//PubmedArticle")
      for (a in articles) {
        xt <- function(xpath) {
          obj <- xml2::xml_find_first(a, xpath)
          if (inherits(obj, "xml_missing")) "" else xml2::xml_text(obj)
        }
        abstract <- paste(xml2::xml_text(xml2::xml_find_all(
          a, ".//Abstract/AbstractText")), collapse = " ")
        doi <- xt(".//ArticleId[@IdType='doi']")
        if (!nzchar(doi)) doi <- xt(".//ELocationID[@EIdType='doi']")
        pmid <- xt(".//MedlineCitation/PMID")
        yr <- xt(".//PubDate/Year")
        if (!nzchar(yr)) yr <- xt(".//PubDate/MedlineDate")
        ix <- ix + 1L
        hits[[ix]] <- new_hit(gene, stratum, "PubMed",
          title = xt(".//Article/ArticleTitle"), abstract = abstract,
          year = yr, doi = doi, pmid = pmid,
          journal = xt(".//Journal/Title"),
          url = paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/"))
      }
    }
    list(hits = dplyr::bind_rows(hits), total = total,
         truncated = total > length(ids))
  }

  search_europepmc <- function(gene, stratum) {
    url <- "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    query <- build_europepmc_query(gene, stratum)
    cursor <- "*"
    total <- NA_integer_
    seen <- 0L
    hits <- list()
    ix <- 0L
    repeat {
      remaining <- CFG$MAX_RESULTS_PER_QUERY - seen
      if (remaining <= 0L) break
      page_size <- min(CFG$PAGE_SIZE, remaining)
      params <- list(query = query, format = "json", resultType = "core",
                     pageSize = page_size, cursorMark = cursor)
      js <- fetch_json(url, params)
      if (is.na(total)) total <- suppressWarnings(as.integer(ch(js$hitCount)))
      rows <- js$resultList$result
      if (is.null(rows) || !length(rows)) break
      for (x in rows) {
        ix <- ix + 1L
        pmid <- ch(x$pmid)
        link <- if (nzchar(pmid))
          paste0("https://pubmed.ncbi.nlm.nih.gov/", pmid, "/") else
          paste0("https://europepmc.org/article/", ch(x$source), "/", ch(x$id))
        hits[[ix]] <- new_hit(gene, stratum, "Europe PMC",
          title = x$title, abstract = x$abstractText, year = x$pubYear,
          doi = x$doi, pmid = pmid, journal = x$journalTitle, url = link)
      }
      seen <- seen + length(rows)
      next_cursor <- ch(js$nextCursorMark)
      if (!nzchar(next_cursor) || identical(next_cursor, cursor)) break
      cursor <- next_cursor
      if (is.finite(total) && seen >= total) break
    }
    list(hits = if (length(hits)) dplyr::bind_rows(hits) else empty_hit(),
         total = total, truncated = is.finite(total) && total > seen)
  }

  if (CFG$RUN_SEARCH) {
    message("Retrieving bibliographic records; keyword matches are NOT verified biological occurrences.")
    engines <- list()
    if (CFG$SEARCH_PUBMED) engines$PubMed <- search_pubmed
    if (CFG$SEARCH_EUROPEPMC) engines$EuropePMC <- search_europepmc
    if (!length(engines)) stop("Enable at least one search engine.")
    all_hits <- list()
    all_logs <- list()
    index <- 0L
    for (gene in GENES) for (stratum in names(QUERIES)) for (db in names(engines)) {
      message("Searching: ", gene, " | ", stratum, " | ", db)
      result <- tryCatch(engines[[db]](gene, stratum), error = function(e) {
        warning("Search failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })
      index <- index + 1L
      if (!is.null(result)) all_hits[[index]] <- result$hits
      all_logs[[index]] <- data.frame(
        gene = gene, query_stratum = stratum, database = db,
        exact_query = if (db == "PubMed") build_pubmed_query(gene, stratum)
                      else build_europepmc_query(gene, stratum),
        date_start = CFG$DATE_START,
        date_end = CFG$DATE_END,
        total_hits_reported = if (is.null(result)) NA_integer_ else result$total,
        retrieved = if (is.null(result)) NA_integer_ else nrow(result$hits),
        truncated = if (is.null(result)) NA else isTRUE(result$truncated),
        status = if (is.null(result)) "ERROR" else "OK",
        search_date = as.character(Sys.Date()), stringsAsFactors = FALSE)
    }
    raw <- dplyr::bind_rows(all_hits)
    logs <- dplyr::bind_rows(all_logs)
    readr::write_csv(logs, FILES$queries, na = "")
    if (nrow(raw)) readr::write_csv(raw, FILES$raw, na = "")
    else warning("No literature records retrieved. Check query log.")
  } else if (file.exists(FILES$raw)) {
    raw <- readr::read_csv(FILES$raw, show_col_types = FALSE)
  } else {
    raw <- empty_hit()
  }

  # This helps find papers but does not upgrade their evidence status.
  if (CFG$IMPORT_PREVIOUS_CANDIDATES &&
      file.exists(CFG$PREVIOUS_CANDIDATES)) {
    prev <- readr::read_csv(CFG$PREVIOUS_CANDIDATES, show_col_types = FALSE)
    required_prev <- c("gene", "title", "doi", "pmid")
    if (all(required_prev %in% names(prev))) {
      getcol <- function(x, name, default = "")
        if (name %in% names(x)) as.character(x[[name]]) else rep(default, nrow(x))
      prev_hits <- data.frame(
        gene = getcol(prev, "gene"), query_stratum = getcol(prev, "query_context"),
        database = "Prior search: PubMed/EuropePMC/OpenAlex/Crossref",
        title = getcol(prev, "title"), abstract = getcol(prev, "abstract"),
        year = getcol(prev, "year"), doi = getcol(prev, "doi"),
        pmid = getcol(prev, "pmid"), journal = getcol(prev, "journal"),
        url = getcol(prev, "url"), stringsAsFactors = FALSE)
      raw <- dplyr::bind_rows(raw, prev_hits)
      message("Included previous bibliographic CANDIDATES (not verified findings).")
    }
  }

  if (nrow(raw)) {
    raw <- raw |>
      dplyr::mutate(
        gene = as.character(gene), doi = doi_clean(doi), pmid = clean(pmid),
        title = clean(title), abstract = clean(abstract),
        year = suppressWarnings(as.integer(substr(as.character(year), 1, 4)))
      ) |>
      dplyr::filter(gene %in% GENES, nonblank(title))
    raw$paper_key <- vapply(seq_len(nrow(raw)), function(i)
      paper_id(raw$doi[i], raw$pmid[i], raw$title[i], raw$year[i]), character(1))
    # Retain each gene-paper association but remove repeat database/query hits.
    candidates <- raw |>
      dplyr::group_by(gene, paper_key) |>
      dplyr::mutate(databases = paste(sort(unique(database)), collapse = " | "),
                    query_strata = paste(sort(unique(query_stratum)), collapse = " | ")) |>
      dplyr::arrange(dplyr::desc(nchar(abstract)), .by_group = TRUE) |>
      dplyr::slice(1) |>
      dplyr::ungroup() |>
      dplyr::arrange(match(gene, GENES), dplyr::desc(year))
    readr::write_csv(candidates, FILES$candidates, na = "")
  } else {
    candidates <- data.frame(
      gene = character(), paper_key = character(), title = character(),
      doi = character(), pmid = character(), year = integer(),
      abstract = character(), url = character(), stringsAsFactors = FALSE)
  }

  # All columns correspond to query STRATA, NOT documented sampling origins.
  # "Unrestricted" is excluded from the plot because it is not a source stratum.
  # No manual screening is performed; the only automated exclusions are missing

  plot_strata <- c("Human", "Avian", "Other animals", "Food",
                   "Wastewater", "Other environment")

  figure_dir <- file.path(CFG$FIGURE_DIR, "figures")
  audit_dir <- file.path(CFG$FIGURE_DIR, "tables")
  dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)

  # The source code above saves all query hits as FILES$raw and query logs as
  # FILES$queries. An exploratory plot DOES NOT use REVIEW_source_ARG_studies.
  if (!file.exists(FILES$raw)) {
    stop("No bibliographic CSV at: ", FILES$raw,
         "\nSet CFG$RUN_SEARCH = TRUE and run again with internet access.")
  }
  if (!file.exists(FILES$queries)) {
    warning("No query log was found: query completeness cannot be established. ",
            "All cells will be flagged as potentially incomplete.")
  }

  hits <- readr::read_csv(FILES$raw, show_col_types = FALSE,
                          col_types = readr::cols(.default = "c"))
  needed_hits <- c("gene", "query_stratum", "database", "title",
                   "year", "doi", "pmid")
  if (!all(needed_hits %in% names(hits))) {
    stop("Missing columns in search hits: ",
         paste(setdiff(needed_hits, names(hits)), collapse = ", "))
  }
  for (column in c("abstract", "journal", "url")) {
    if (!(column %in% names(hits))) hits[[column]] <- ""
  }

  # Keep a copy BEFORE title/year filtering for reproducibility.
  readr::write_csv(hits, file.path(audit_dir, "01_raw_retrieved_records.csv"),
                   na = "")

  hits$gene <- clean(hits$gene)
  hits$query_stratum <- clean(hits$query_stratum)
  hits$title <- clean(hits$title)
  hits$year_int <- suppressWarnings(as.integer(substr(hits$year, 1, 4)))
  hits$doi <- doi_clean(hits$doi)
  hits$pmid <- clean(hits$pmid)
  hits$database <- clean(hits$database)
  # These are year-level checks only; the exact dates are applied in API queries.
  # All comparisons are length nrow(hits), even if the input table has 0 rows.
  if (length(CFG$MIN_YEAR) != 1L || length(CFG$MAX_YEAR) != 1L ||
      is.na(CFG$MIN_YEAR) || is.na(CFG$MAX_YEAR)) {
    stop("CFG$MIN_YEAR and CFG$MAX_YEAR must be single, nonmissing years.")
  }
  hits$included_metadata <-
    (hits$gene %in% GENES) &
    (hits$query_stratum %in% names(QUERIES)) &
    nonblank(hits$title) &
    !is.na(hits$year_int) &
    (hits$year_int >= CFG$MIN_YEAR) &
    (hits$year_int <= CFG$MAX_YEAR)
  hits$included_metadata[is.na(hits$included_metadata)] <- FALSE
  readr::write_csv(hits, file.path(audit_dir, "02_metadata_screening_audit.csv"),
                   na = "")
  hits <- hits[hits$included_metadata, , drop = FALSE]

  # can guarantee that similarly titled items are the exact same publication.
  canonical_title <- function(x) {
    x <- tolower(clean(x))
    x <- gsub("[^[:alnum:]]+", " ", x)
    trimws(gsub("[[:space:]]+", " ", x))
  }
  if (nrow(hits) == 0L) stop("No retained records after metadata filtering.")
  hits$publication_key <- paste0(
    "titleyear:", vapply(seq_len(nrow(hits)), function(i)
      digest::digest(paste(canonical_title(hits$title[i]),
                           ifelse(is.na(hits$year_int[i]), "unknown",
                                  hits$year_int[i]), sep = "|"), algo = "sha1"),
      character(1)))

  # Deduplicate WITHIN each gene-query-stratum combination. Do not deduplicate
  # across source-query strata, because the very same record can legitimately
  # match more than one set of search terms. Count publications, not API hits.
  first_populated <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & nzchar(trimws(x))]
    if (length(x)) x[1] else ""
  }

  unique_records <- hits |>
    dplyr::group_by(gene, query_stratum, publication_key) |>
    dplyr::summarise(
      title = first_populated(title),
      publication_year = first_populated(as.character(year_int)),
      doi = first_populated(doi),
      pmid = first_populated(pmid),
      databases = paste(sort(unique(database)), collapse = " | "),
      database_record_count = dplyr::n(),
      url = first_populated(url),
      .groups = "drop"
    ) |>
    dplyr::arrange(match(gene, GENES), match(query_stratum, names(QUERIES)),
                   dplyr::desc(publication_year))

  readr::write_csv(unique_records,
                   file.path(audit_dir, "03_unique_retrieved_publications_by_query.csv"),
                   na = "")

  # Query completeness: a missing/error/truncated search must never be silently
  # em dash, and their raw counts are lower bounds of RETRIEVED records only.
  query_grid <- expand.grid(gene = GENES, query_stratum = plot_strata,
                            stringsAsFactors = FALSE)

  if (file.exists(FILES$queries)) {
    logs <- readr::read_csv(FILES$queries, show_col_types = FALSE,
                            col_types = readr::cols(.default = "c"))
  } else {
    logs <- data.frame(gene = character(), query_stratum = character(),
                       database = character(), status = character(),
                       truncated = character(), search_date = character())
  }
  required_log <- c("gene", "query_stratum", "database", "status", "truncated")
  if (!all(required_log %in% names(logs))) {
    warning("Incomplete query log schema; flagging cells as incomplete.")
    logs <- data.frame(gene = character(), query_stratum = character(),
                       database = character(), status = character(),
                       truncated = character(), search_date = character())
  }
  if (!("search_date" %in% names(logs))) logs$search_date <- ""
  if (!("total_hits_reported" %in% names(logs))) logs$total_hits_reported <- ""
  if (!("retrieved" %in% names(logs))) logs$retrieved <- ""
  if (!("exact_query" %in% names(logs))) logs$exact_query <- ""

  # Infer only databases whose query log is available; if none are recorded,
  # the results remain incomplete (not zero). Check every expected database.
  expected_databases <- c(if (CFG$SEARCH_PUBMED) "PubMed",
                          if (CFG$SEARCH_EUROPEPMC) "EuropePMC")
  if (!length(expected_databases)) expected_databases <- unique(logs$database)

  query_grid$log_complete <- vapply(seq_len(nrow(query_grid)), function(i) {
    subset_logs <- logs[logs$gene == query_grid$gene[i] &
      logs$query_stratum == query_grid$query_stratum[i] &
      logs$database %in% expected_databases, , drop = FALSE]
    if (!length(expected_databases) || nrow(subset_logs) != length(expected_databases))
      return(FALSE)
    if (!("date_start" %in% names(subset_logs)) ||
        !("date_end" %in% names(subset_logs))) return(FALSE)
    all(!is.na(subset_logs$status) & subset_logs$status == "OK" &
          !is.na(subset_logs$truncated) &
          tolower(subset_logs$truncated) %in% c("false", "f", "0") &
          !is.na(subset_logs$date_start) &
          !is.na(subset_logs$date_end) &
          subset_logs$date_start == CFG$DATE_START &
          subset_logs$date_end == CFG$DATE_END)
  }, logical(1))
  query_grid$incomplete <- !query_grid$log_complete

  count_df <- unique_records |>
    dplyr::filter(query_stratum %in% plot_strata) |>
    dplyr::count(gene, query_stratum, name = "n_records")
  plot_data <- dplyr::left_join(query_grid, count_df,
                                by = c("gene", "query_stratum"))
  plot_data$n_records[is.na(plot_data$n_records)] <- 0L

  plot_data$cell_label <- ifelse(
    plot_data$incomplete & plot_data$n_records == 0,
    "—",
    as.character(plot_data$n_records)
  )
  plot_data$plot_value <- ifelse(plot_data$incomplete &
                                  plot_data$n_records == 0,
                                NA_real_, plot_data$n_records)

  plot_data$query_stratum <- factor(plot_data$query_stratum,
                                    levels = plot_strata)
  plot_data$gene <- factor(plot_data$gene, levels = rev(GENES))
  plot_data$source_display <- factor(as.character(plot_data$query_stratum),
    levels = plot_strata,
    labels = c("Human", "Avian", "Other\nanimals", "Food",
               "Wastewater", "Other\nenvironment"))

  upper <- max(plot_data$n_records, 1L, na.rm = TRUE)
  plot_data$text_tone <- ifelse(!is.na(plot_data$plot_value) &
                                 plot_data$plot_value >= upper * 0.62,
                               "white", "dark")
  readr::write_csv(plot_data,
                   file.path(audit_dir, "04_MAIN_FIGURE_heatmap_counts_and_flags.csv"),
                   na = "")
  readr::write_csv(logs,
                   file.path(audit_dir, "05_query_completeness_log.csv"),
                   na = "")

  # This heatmap shows keyword-based RETRIEVAL, not verified sample sources.
  n_cells <- nrow(plot_data)
  n_partial <- sum(plot_data$incomplete, na.rm = TRUE)
  plot_subtitle <- paste0(
    "Deduplicated bibliographic records from ", CFG$DATE_START,
    " to ", CFG$DATE_END, " (source-context keyword queries)")
  plot_note <- paste0(
    "Cells show deduplicated retrieved records, not confirmed source-gene associations. ",
    "Orange outline = truncated, failed, or unauditable query; ",
    "a single publication can appear in multiple columns. ",
    "Gene labels are provisional pending annotation verification.")

  p <- ggplot2::ggplot(plot_data,
                       ggplot2::aes(x = source_display, y = gene,
                                    fill = plot_value)) +
    ggplot2::geom_tile(
      data = plot_data[!plot_data$incomplete, , drop = FALSE],
      colour = "#FFFFFF", linewidth = 1.8, width = 0.96, height = 0.91
    ) +
    ggplot2::geom_tile(
      data = plot_data[plot_data$incomplete, , drop = FALSE],
      colour = "#E58C44", linewidth = 1.4, width = 0.96, height = 0.91
    ) +
    ggplot2::geom_text(ggplot2::aes(label = cell_label, colour = text_tone),
                       size = 5.3, fontface = "bold", show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(dark = "#18364B",
                                              white = "#FFFFFF"),
                                 guide = "none") +
    ggplot2::scale_fill_gradientn(
      colours = c("#F3F8FA", "#BFDDE7", "#6CB8CE", "#2288AB", "#064968"),
      na.value = "#EBEEF0", name = "Retrieved\nrecords",
      limits = c(0, upper), breaks = pretty(c(0, upper), n = 5)
    ) +
    ggplot2::scale_x_discrete(position = "top", expand = c(0.02, 0)) +
    ggplot2::scale_y_discrete(expand = c(0.02, 0)) +
    ggplot2::labs(
      title = paste0("Literature retrieval landscape of ", length(GENES), " prioritized ARG indicators"),
      subtitle = plot_subtitle,
      x = NULL, y = NULL, caption = plot_note
    ) +
    ggplot2::coord_fixed(ratio = 0.65, clip = "off") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 12, colour = "#1D4055",
                                          face = "bold", lineheight = 0.96,
                                          margin = ggplot2::margin(b = 10)),
      axis.text.y = ggplot2::element_text(size = 12, colour = "#1D4055",
                                          face = "bold", margin = ggplot2::margin(r = 10)),
      plot.title = ggplot2::element_text(size = 18, colour = "#123247",
                                         face = "bold", margin = ggplot2::margin(b = 8)),
      plot.subtitle = ggplot2::element_text(size = 10.5, colour = "#4C6475",
                                            margin = ggplot2::margin(b = 18)),
      plot.caption = ggplot2::element_text(size = 8.5, colour = "#536D7D",
                                           lineheight = 1.05,
                                           margin = ggplot2::margin(t = 15)),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 11),
      legend.text = ggplot2::element_text(size = 10),
      plot.margin = ggplot2::margin(18, 25, 20, 18)
    )

  figure_base <- file.path(figure_dir, "FIGURE_MAIN_ARG_literature_retrieval_heatmap")
  ggplot2::ggsave(paste0(figure_base, ".png"), plot = p,
    width = CFG$WIDTH_IN, height = CFG$HEIGHT_IN, units = "in",
    dpi = CFG$DPI, bg = "white", device = ragg::agg_png)
  ggplot2::ggsave(paste0(figure_base, ".tiff"), plot = p,
    width = CFG$WIDTH_IN, height = CFG$HEIGHT_IN, units = "in",
    dpi = CFG$DPI, compression = "lzw", bg = "white", device = "tiff")
  ggplot2::ggsave(paste0(figure_base, ".svg"), plot = p,
    width = CFG$WIDTH_IN, height = CFG$HEIGHT_IN, units = "in",
    bg = "white", device = svglite::svglite)
  ggplot2::ggsave(paste0(figure_base, ".pdf"), plot = p,
    width = CFG$WIDTH_IN, height = CFG$HEIGHT_IN, units = "in", bg = "white")

  cat("Flagged incomplete search cells: ", n_partial, " / ", n_cells, "\n")
  if (n_partial > 0L) warning("Some search cells are incomplete/capped; see orange outlines and query log.")
  print(p)

  invisible(TRUE)
}

# plot_literature_bars keeps the original search or plotting rules.
plot_literature_bars <- function() {
  ROOT <- "ARG_literature_search_2016_2026_exact"
  INPUT <- file.path(ROOT, "tables", "04_MAIN_FIGURE_heatmap_counts_and_flags.csv")
  OUTDIR <- file.path(ROOT, "figures")
  DPI <- 300

  packages <- c("ggplot2", "readr", "dplyr", "scales")
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages)) {
    stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
  }

  if (!requireNamespace("ragg", quietly = TRUE)) {
    message("Optional package 'ragg' not installed; using base PNG device.")
  }
  if (!requireNamespace("svglite", quietly = TRUE)) {
    message("Optional package 'svglite' not installed; SVG output will be skipped.")
  }

  # ------------------------- 3. READ INPUT ------------------------------------
  if (!file.exists(INPUT)) {
    stop("Could not find input file:\n", INPUT,
         "\nRun the literature heatmap script first.")
  }

  dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

  dat <- readr::read_csv(INPUT, show_col_types = FALSE)

  need <- c("gene", "query_stratum", "n_records", "incomplete")
  if (!all(need %in% names(dat))) {
    stop("Input CSV must contain: ", paste(need, collapse = ", "),
         "\nMissing: ", paste(setdiff(need, names(dat)), collapse = ", "))
  }

  GENES <- literature_genes()
  STRATA <- c("Human", "Avian", "Other animals", "Food", "Wastewater", "Other environment")

  PAL_GENE <- literature_palette(GENES)

  to_flag <- function(x) {
    if (is.logical(x)) return(replace(x, is.na(x), TRUE))
    value <- trimws(tolower(as.character(x)))
    value[is.na(value)] <- ""
    !(value %in% c("false", "f", "0", "no"))
  }

  dat <- dat |>
    dplyr::mutate(
      gene = as.character(gene),
      query_stratum = as.character(query_stratum),
      n_records = suppressWarnings(as.numeric(n_records)),
      incomplete = to_flag(incomplete)
    ) |>
    dplyr::filter(gene %in% GENES, query_stratum %in% STRATA)

  if (anyDuplicated(paste(dat$gene, dat$query_stratum))) {
    stop("Duplicate gene/query-stratum rows found. Check input table.")
  }

  expected <- expand.grid(
    gene = GENES,
    query_stratum = STRATA,
    stringsAsFactors = FALSE
  )

  dat <- dplyr::left_join(expected, dat,
                          by = c("gene", "query_stratum")) |>
    dplyr::mutate(
      incomplete = ifelse(is.na(incomplete), TRUE, incomplete),
      n_records = ifelse(is.na(n_records) & !incomplete, 0, n_records),
      available = !is.na(n_records),
      value = ifelse(available, n_records, NA_real_),
      label = ifelse(available, scales::comma(n_records), "n/a"),
      gene = factor(gene, levels = rev(GENES)),
      query_stratum = factor(query_stratum, levels = STRATA)
    )

  if (!any(dat$available)) {
    stop("No numeric counts available to plot.")
  }

  max_count <- max(dat$value, na.rm = TRUE)
  label_pad <- max(35, max_count * 0.035)
  x_limit <- max_count + label_pad * 3.5

  plot_data <- dat[dat$available, , drop = FALSE]
  missing_data <- dat[!dat$available, , drop = FALSE]

  p <- ggplot2::ggplot() +
    ggplot2::geom_col(
      data = plot_data,
      ggplot2::aes(x = value, y = gene, fill = gene),
      width = 0.68,
      colour = NA
    ) +
    ggplot2::geom_text(
      data = plot_data,
      ggplot2::aes(x = value + label_pad, y = gene, label = label),
      hjust = 0,
      size = 3.4,
      fontface = "bold",
      colour = "#18384D"
    ) +
    ggplot2::geom_text(
      data = missing_data,
      ggplot2::aes(x = label_pad, y = gene, label = "n/a"),
      hjust = 0,
      size = 3.35,
      colour = "#7A7A7A"
    ) +
    ggplot2::scale_fill_manual(values = PAL_GENE, drop = FALSE, guide = "none") +
    ggplot2::facet_wrap(~query_stratum, ncol = 3) +
    ggplot2::scale_x_continuous(
      limits = c(0, x_limit),
      breaks = pretty(c(0, max_count), n = 4),
      labels = scales::comma,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      x = "Retrieved records",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        face = "bold",
        colour = "#173B51",
        size = 12
      ),
      strip.background = ggplot2::element_rect(
        fill = "#ECF3F6",
        colour = NA
      ),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(
        colour = "#E4ECF0",
        linewidth = 0.30
      ),
      axis.text.y = ggplot2::element_text(
        colour = "#173B51",
        size = 10.5
      ),
      axis.text.x = ggplot2::element_text(
        colour = "#516776",
        size = 9
      ),
      axis.title.x = ggplot2::element_text(
        colour = "#314A5A",
        size = 11,
        margin = ggplot2::margin(t = 10)
      ),
      panel.spacing = grid::unit(1.2, "lines"),
      plot.margin = ggplot2::margin(15, 18, 15, 18)
    )

  print(p)

  # ------------------------- 7. EXPORT ----------------------------------------
  base <- file.path(OUTDIR, "FIGURE_MAIN_ARG_literature_faceted_barplot_by_gene")

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      paste0(base, ".png"),
      p,
      width = 14,
      height = 9,
      units = "in",
      dpi = DPI,
      bg = "white",
      device = ragg::agg_png
    )
  } else {
    ggplot2::ggsave(
      paste0(base, ".png"),
      p,
      width = 14,
      height = 9,
      units = "in",
      dpi = DPI,
      bg = "white"
    )
  }

  ggplot2::ggsave(
    paste0(base, ".tiff"),
    p,
    width = 14,
    height = 9,
    units = "in",
    dpi = DPI,
    compression = "lzw",
    bg = "white",
    device = "tiff"
  )

  ggplot2::ggsave(
    paste0(base, ".pdf"),
    p,
    width = 14,
    height = 9,
    units = "in",
    bg = "white"
  )

  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(
      paste0(base, ".svg"),
      p,
      width = 14,
      height = 9,
      units = "in",
      bg = "white",
      device = svglite::svglite
    )
  }

  readr::write_csv(dat, paste0(base, "_counts_and_flags.csv"), na = "")

  cat("\nSaved figure:\n", base, "\n")

  invisible(TRUE)
}

# plot_literature_two_panel keeps the original search or plotting rules.
plot_literature_two_panel <- function() {
  library(patchwork)
  # FIGURE: Computational literature mapping of five prioritized ARG labels
  # Research in Veterinary Science | Two-panel manuscript figure
  # B: Five ARGs x six source-context keyword queries (record counts)
  # INPUT: FIGURE_MAIN_ARG_literature_faceted_barplot_by_gene_counts_and_flags.csv
  # SVG is only a rendered plot; do NOT infer data from SVG dimensions.
  # Data are keyword-query retrieval yields, NOT confirmed biological occurrences.
  # Publication year per record is absent: DO NOT plot annual trends from this CSV.

  DATA_FILE <- "FIGURE_MAIN_ARG_literature_faceted_barplot_by_gene_counts_and_flags.csv"
  OUTPUT_DIR <- "RVS_literature_mapping_figure"
  FIGURE_NAME <- "Figure_literature_mapping_two_panel"
  DPI <- 300L
  WIDTH_IN <- 13.8
  HEIGHT_IN <- 8.5
  DATE_START <- "2016-10-15"
  DATE_END <- "2026-10-15"

  # Use exact IDs from the original table; final ARG identity verification remains
  # necessary before claiming the allele annotations are biologically confirmed.
  GENES <- literature_genes()
  STRATA <- c("Human", "Avian", "Other animals", "Food",
              "Wastewater", "Other environment")

  PAL_GENE <- literature_palette(GENES)

  needed <- c("ggplot2", "dplyr", "readr", "scales", "patchwork")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))

  if (!file.exists(DATA_FILE)) {
    backup <- file.path("ARG_literature_search_2016_2026_exact", "tables",
                        "04_MAIN_FIGURE_heatmap_counts_and_flags.csv")
    if (file.exists(backup)) DATA_FILE <- backup
  }
  if (!file.exists(DATA_FILE)) {
    stop("Cannot locate the gene x context count CSV. Place it in getwd(): ", getwd(),
         "\nExpected filename: ", basename(DATA_FILE),
         "\nThe SVG image alone does not contain a dependable analysis table.")
  }
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

  # 3. LOAD + VALIDATE THE EXISTING SEARCH YIELDS -------------------------------
  dat <- readr::read_csv(DATA_FILE, show_col_types = FALSE)
  required <- c("gene", "query_stratum", "n_records", "incomplete")
  if (!all(required %in% colnames(dat))) {
    stop("Missing CSV fields: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }

  parse_incomplete <- function(x) {
    if (is.logical(x)) return(replace(x, is.na(x), TRUE))
    v <- tolower(trimws(as.character(x)))
    is.na(v) | !(v %in% c("false", "f", "0", "no"))
  }

  dat <- dat |>
    dplyr::transmute(
      gene = as.character(gene),
      query_stratum = as.character(query_stratum),
      n_records = suppressWarnings(as.numeric(n_records)),
      incomplete = parse_incomplete(incomplete)
    ) |>
    dplyr::filter(gene %in% GENES, query_stratum %in% STRATA)

  if (anyDuplicated(paste(dat$gene, dat$query_stratum, sep = "||"))) {
    stop("Duplicate gene x query combinations detected; inspect the input CSV.")
  }

  expected <- expand.grid(gene = GENES, query_stratum = STRATA,
                          stringsAsFactors = FALSE)
  dat <- dplyr::left_join(expected, dat, by = c("gene", "query_stratum"))
  # If a combination is absent from the CSV, do NOT convert it to zero.
  dat$incomplete[is.na(dat$incomplete)] <- TRUE
  n_cells <- length(GENES) * length(STRATA)
  if (nrow(dat) != n_cells) stop("Expected ", n_cells, " gene x query cells.")
  if (!any(is.finite(dat$n_records))) stop("No numeric counts were found.")
  if (any(is.finite(dat$n_records) & dat$n_records < 0)) {
    stop("Negative bibliographic counts are invalid.")
  }

  n_incomplete <- sum(dat$incomplete)
  if (n_incomplete) {
    message("QUERY AUDIT: ", n_incomplete, " of ", n_cells, " cells are incomplete.")
    print(dat[dat$incomplete,
              c("gene", "query_stratum", "n_records", "incomplete"),
              drop = FALSE])
  }

  readr::write_csv(dat, file.path(OUTPUT_DIR, "gene_context_counts_and_flags.csv"), na = "")

  steps <- data.frame(
    x = seq_len(5),
    label = c(
      paste0(length(GENES), " network-\nprioritized ARGs"),
      "PubMed and\nEurope PMC",
      "2016–2026\ndate-restricted search",
      "Title–year\ndeduplication + audit",
      paste0(n_cells, " gene–context\nretrieval summaries")
    ),
    stringsAsFactors = FALSE
  )

  arrows <- data.frame(
    x = 1:4 + 0.46,
    xend = 2:5 - 0.46,
    y = 0,
    yend = 0
  )

  panel_a <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = steps,
      ggplot2::aes(xmin = x - 0.42, xmax = x + 0.42,
                   ymin = -0.52, ymax = 0.52),
      fill = "#F3F7FA", colour = "#527D98", linewidth = 0.65
    ) +
    ggplot2::geom_segment(
      data = arrows,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
      colour = "#397CA1", linewidth = 0.8,
      arrow = grid::arrow(type = "closed", length = grid::unit(0.13, "cm"))
    ) +
    ggplot2::geom_text(
      data = steps,
      ggplot2::aes(x = x, y = 0, label = label),
      size = 3.7, colour = "#1B4057", fontface = "bold", lineheight = 1.06
    ) +
    ggplot2::coord_cartesian(xlim = c(0.48, 5.52), ylim = c(-0.65, 0.65),
                             expand = FALSE, clip = "off") +
    ggplot2::labs(title = "A   Computational literature-mapping workflow") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#12374D",
                                         size = 14, margin = ggplot2::margin(b = 8)),
      plot.margin = ggplot2::margin(8, 18, 0, 18)
    )

  # The six columns represent SEARCH KEYWORDS, NOT verified sampling origins.
  # Bubble area/diameter is a visual encoding of retrieval volume, not ARG risk.
  dat$gene <- factor(dat$gene, levels = rev(GENES))
  dat$query_stratum <- factor(dat$query_stratum, levels = STRATA)
  available <- dat[is.finite(dat$n_records), , drop = FALSE]
  unavailable <- dat[!is.finite(dat$n_records), , drop = FALSE]

  panel_b <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = dat,
      ggplot2::aes(x = query_stratum, y = gene),
      fill = "#F6F9FB", colour = "#DFE9EF", linewidth = 0.55,
      width = 0.93, height = 0.88
    ) +
    ggplot2::geom_point(
      data = available,
      ggplot2::aes(x = query_stratum, y = gene, size = n_records,
                   fill = gene, colour = gene),
      shape = 21, stroke = 0.95, alpha = 0.60
    ) +
    ggplot2::geom_text(
      data = available,
      ggplot2::aes(x = query_stratum, y = gene,
                   label = scales::comma(n_records, accuracy = 1)),
      size = 4.1, fontface = "bold", colour = "#16374B"
    ) +
    ggplot2::geom_text(
      data = unavailable,
      ggplot2::aes(x = query_stratum, y = gene),
      label = "NA", colour = "#5C7080", size = 3.7
    ) +
    ggplot2::scale_size_continuous(
      range = c(12, 26), guide = "none"
    ) +
    ggplot2::scale_fill_manual(values = PAL_GENE, guide = "none") +
    ggplot2::scale_colour_manual(values = PAL_GENE, guide = "none") +
    ggplot2::scale_x_discrete(
      position = "top",
      labels = c("Human", "Avian", "Other\nanimals", "Food",
                 "Wastewater", "Other\nenvironment"),
      expand = ggplot2::expansion(add = 0.55)
    ) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.42)) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = "B   Bibliographic retrieval by ARG and keyword context",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 11.5, face = "bold",
                                         colour = "#173B51", lineheight = 0.98,
                                         margin = ggplot2::margin(b = 9)),
      axis.text.y = ggplot2::element_text(size = 12.5, face = "bold",
                                         colour = "#173B51",
                                         margin = ggplot2::margin(r = 9)),
      plot.title = ggplot2::element_text(size = 14, face = "bold",
                                         colour = "#12374D",
                                         margin = ggplot2::margin(b = 13)),
      plot.margin = ggplot2::margin(0, 23, 5, 14)
    )

  # No overarching title or in-figure warning note; place explanation in caption.
  p <- panel_a / panel_b +
    patchwork::plot_layout(heights = c(1.1, 3.0)) &
    ggplot2::theme(plot.background = ggplot2::element_rect(
      fill = "white", colour = NA))

  print(p)

  # 7. EXPORT FIGURE AT 300 DPI + VECTOR FORMATS --------------------------------
  base <- file.path(OUTPUT_DIR, FIGURE_NAME)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(paste0(base, ".png"), p,
                    width = WIDTH_IN, height = HEIGHT_IN,
                    units = "in", dpi = DPI, bg = "white", device = ragg::agg_png)
  } else {
    ggplot2::ggsave(paste0(base, ".png"), p,
                    width = WIDTH_IN, height = HEIGHT_IN,
                    units = "in", dpi = DPI, bg = "white")
  }

  ggplot2::ggsave(paste0(base, ".tiff"), p,
                  width = WIDTH_IN, height = HEIGHT_IN,
                  units = "in", dpi = DPI, compression = "lzw",
                  device = "tiff", bg = "white")
  ggplot2::ggsave(paste0(base, ".pdf"), p,
                  width = WIDTH_IN, height = HEIGHT_IN,
                  units = "in", bg = "white")
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(paste0(base, ".svg"), p,
                    width = WIDTH_IN, height = HEIGHT_IN,
                    units = "in", bg = "white", device = svglite::svglite)
  }

  # The following manuscript caption is accurate for this input and figure.
  caption <- paste0(
    "Figure X. Computational literature mapping of ", length(GENES), " network-prioritized ",
    "antimicrobial resistance gene (ARG) labels. (A) Bibliographic retrieval ",
    "and analysis workflow implemented in R, including date-restricted ",
    "queries of PubMed and Europe PMC, provisional title-year deduplication ",
    "and query auditing. (B) Numbers of retrieved bibliographic records for ",
    "each gene and source-context keyword stratum, 15 October 2016–15 October ",
    "2026. Colors identify gene labels; bubble size indicates retrieval volume. ",
    "Publications may be included in multiple strata. ",
    if (n_incomplete) paste0("Retrieval reached the per-query limit or was incomplete for ",
      paste(unique(paste0(dat$gene[dat$incomplete], " (", dat$query_stratum[dat$incomplete], ")")), collapse = ", "), ". ")
    else "", "Counts represent ",
    "bibliographic search yields rather than verified ARG occurrence, host ",
    "identity, prevalence, transmission or clinical risk. Gene-to-reference ",
    "sequence annotation remains subject to verification."
  )
  writeLines(caption, file.path(OUTPUT_DIR, "Figure_caption.txt"), useBytes = TRUE)
  cat("\nSaved figure files and audit table in: ", OUTPUT_DIR,
      "\nIncomplete query-context cells: ", n_incomplete,
      " / ", nrow(dat),
      "\nThese counts DO NOT establish annual publication growth or verified ARG occurrence.\n",
      sep = "")

  invisible(TRUE)
}

# Run search first, then the two figure types.
run_literature_mapping <- function() {
  run_literature_search()
  plot_literature_bars()
  plot_literature_two_panel()
  invisible(TRUE)
}
