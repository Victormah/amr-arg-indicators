# =============================================================================
# scripts/build_card_mapping.R - rebuild the ResFinder-to-CARD mapping table
# Project : Candidate antimicrobial resistance gene indicators for veterinary
#           public health surveillance (Mahawan et al., Research in Veterinary
#           Science, submitted)
# Purpose : Rebuild data/AMRgenelist_mapped_CARD_ALL.csv from the CARD ontology
#           (data/aro.tsv) and the online CARD database. Optional: run only if
#           you want to regenerate the mapping instead of downloading it.
# Inputs  : data/aro.tsv (CARD download)
# Outputs : CARD_full_mapping/ (copy the final CSV to data/)
# Run     : source("scripts/build_card_mapping.R") from the project folder
# =============================================================================


# ============================================================
# Complete CARD annotation mapping from the online database
# Inputs: data/aro.tsv
# Outputs: CARD_full_mapping/
# ============================================================

# Install required packages once:
# install.packages(c("jsonlite", "writexl"))

library(jsonlite)

# ---------- Helper functions ----------

card_text <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)

  y <- tryCatch(
    as.character(unlist(x, use.names = FALSE)),
    error = function(e) character()
  )

  y <- trimws(y)
  y <- y[!is.na(y) & nzchar(y)]

  if (length(y)) y[1L] else NA_character_
}

card_join <- function(x) {
  if (is.null(x) || !length(x)) return(NA_character_)

  x <- trimws(as.character(x))
  x <- sort(unique(x[!is.na(x) & nzchar(x)]))

  if (length(x)) paste(x, collapse = "; ") else NA_character_
}

card_norm_acc <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  sub("\\.[0-9]+$", "", x)
}

card_first <- function(primary, backup) {
  primary <- as.character(primary)
  backup <- as.character(backup)

  replace <- is.na(primary) | !nzchar(trimws(primary))
  primary[replace] <- backup[replace]

  primary
}

card_write <- function(x, path) {
  write.csv(
    x,
    path,
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

# ---------- Download CARD archives ----------

card_download <- function(url, destination, refresh = TRUE) {
  dir.create(dirname(destination),
             recursive = TRUE, showWarnings = FALSE)

  if (refresh || !file.exists(destination)) {
    message("Downloading: ", url)

    old_timeout <- getOption("timeout")
    options(timeout = max(1200, old_timeout))
    on.exit(options(timeout = old_timeout), add = TRUE)

    download.file(
      url,
      destfile = destination,
      mode = "wb",
      method = "libcurl"
    )
  } else {
    message("Using cached file: ", destination)
  }

  if (!file.exists(destination) ||
      file.size(destination) < 10000) {
    stop("CARD archive is missing or too small.")
  }

  destination
}

card_extract <- function(archive, filename, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  members <- untar(archive, list = TRUE)
  selected <- members[basename(members) == filename]

  if (length(selected) != 1L) {
    stop("Could not uniquely locate ", filename,
         " in ", archive)
  }

  untar(archive, files = selected, exdir = out_dir)

  extracted <- file.path(out_dir, selected)
  destination <- file.path(out_dir, filename)

  if (!file.exists(extracted)) {
    stop("Extraction failed for ", filename)
  }

  if (normalizePath(extracted, mustWork = TRUE) !=
      normalizePath(destination, mustWork = FALSE)) {
    file.copy(extracted, destination, overwrite = TRUE)
  }

  destination
}

# ---------- Read original aro.tsv ----------

card_read_aro <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)

  x <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )

  names(x) <- sub("^\ufeff", "", names(x))

  required <- c(
    "Accession",
    "Name",
    "Description",
    "ID",
    "CARD Short Name"
  )

  if (!all(required %in% names(x))) {
    stop("Required columns missing from aro.tsv")
  }

  x <- x[, required]

  names(x) <- c(
    "aro_accession",
    "aro_name_input",
    "description_input",
    "ontology_numeric_id",
    "card_short_name"
  )

  x$aro_accession <- trimws(x$aro_accession)

  if (anyNA(x$aro_accession) ||
      anyDuplicated(x$aro_accession)) {
    stop("ARO accessions are missing or duplicated.")
  }

  x
}

# ---------- Parse CARD ontology ----------

card_parse_obo <- function(path) {
  message("Parsing CARD ontology...")

  lines <- readLines(
    path, warn = FALSE, encoding = "UTF-8"
  )

  starts <- which(trimws(lines) == "[Term]")

  boundaries <- which(
    grepl("^\\[[^]]+\\]$", trimws(lines))
  )

  terms <- vector("list", length(starts))
  counter <- 0L

  for (st in starts) {
    next_boundary <- boundaries[boundaries > st]

    en <- if (length(next_boundary)) {
      next_boundary[1L] - 1L
    } else {
      length(lines)
    }

    block <- if (st < en) {
      lines[(st + 1L):en]
    } else {
      character()
    }

    extract <- function(prefix) {
      substring(
        block[startsWith(block, prefix)],
        nchar(prefix) + 1L
      )
    }

    id <- extract("id: ")

    if (!length(id) ||
        !grepl("^ARO:[0-9]+$", id[1L])) {
      next
    }

    nm <- extract("name: ")
    parents <- sub(" .*$", "", extract("is_a: "))
    relations <- extract("relationship: ")
    obsolete <- extract("is_obsolete: ")

    counter <- counter + 1L

    terms[[counter]] <- list(
      id = id[1L],
      name = if (length(nm)) nm[1L] else NA_character_,
      parents = as.character(parents),
      rels = as.character(relations),
      obsolete = length(obsolete) > 0L &&
        identical(obsolete[1L], "true")
    )
  }

  terms <- terms[seq_len(counter)]

  names(terms) <- vapply(
    terms, `[[`, character(1), "id"
  )

  terms
}

# ---------- Derive drug classes and mechanisms ----------

card_ontology_categories <- function(obo, ids) {
  labels <- vapply(
    obo, `[[`, character(1), "name"
  )

  parents <- lapply(obo, `[[`, "parents")
  rels <- lapply(obo, `[[`, "rels")

  names(parents) <- names(rels) <- names(obo)

  ancestor_cache <- new.env(parent = emptyenv())

  ancestors <- function(id) {
    if (is.na(id) || !nzchar(id)) {
      return(character())
    }

    if (exists(
      id, envir = ancestor_cache, inherits = FALSE
    )) {
      return(get(
        id, envir = ancestor_cache,
        inherits = FALSE
      ))
    }

    seen <- character()
    todo <- id

    while (length(todo)) {
      current <- todo[1L]
      todo <- todo[-1L]

      if (is.na(current) ||
          current %in% seen) next

      seen <- c(seen, current)

      next_ids <- parents[[current]]

      if (!is.null(next_ids) && length(next_ids)) {
        next_ids <- as.character(next_ids)
        next_ids <- next_ids[!is.na(next_ids)]
        todo <- c(todo, next_ids)
      }
    }

    assign(id, seen, envir = ancestor_cache)
    seen
  }

  # FIX: safely handle missing ontology relationships.
  relationship_targets <- function(ids, relation) {
    if (!length(ids)) {
      return(character())
    }

    r <- unlist(
      rels[ids],
      use.names = FALSE
    )

    if (is.null(r) || !length(r)) {
      return(character())
    }

    r <- as.character(r)
    r <- r[!is.na(r)]

    r <- r[
      startsWith(r, paste0(relation, " "))
    ]

    if (!length(r)) {
      return(character())
    }

    target <- sub(
      "^[^ ]+ (ARO:[0-9]+).*$",
      "\\1",
      r
    )

    unique(
      target[grepl("^ARO:[0-9]+$", target)]
    )
  }

  drug_class_one <- function(id) {
    a <- ancestors(id)

    direct <- relationship_targets(
      a, "confers_resistance_to_drug_class"
    )

    if (!length(direct)) {
      drugs <- relationship_targets(
        a, "confers_resistance_to_antibiotic"
      )

      direct <- unique(unlist(
        lapply(drugs, function(d) {
          aa <- ancestors(d)

          aa[vapply(aa, function(t) {
            "ARO:1000003" %in% parents[[t]]
          }, logical(1))]
        }),
        use.names = FALSE
      ))
    }

    card_join(unname(labels[direct]))
  }

  mechanism_one <- function(id) {
    a <- ancestors(id)

    reactions <- relationship_targets(
      a, "participates_in"
    )

    target <- unique(unlist(
      lapply(reactions, function(d) {
        aa <- ancestors(d)

        aa[vapply(aa, function(t) {
          "ARO:1000002" %in% parents[[t]]
        }, logical(1))]
      }),
      use.names = FALSE
    ))

    card_join(unname(labels[target]))
  }

  message(
    "Deriving CARD annotations for ",
    length(ids), " ARO entries..."
  )

  data.frame(
    aro_accession = ids,

    ontology_drug_class = vapply(
      ids, drug_class_one, character(1)
    ),

    ontology_mechanism = vapply(
      ids, mechanism_one, character(1)
    ),

    obo_name = unname(labels[ids]),

    ontology_obsolete = vapply(ids, function(z) {
      if (is.null(obo[[z]])) return(NA)
      obo[[z]]$obsolete
    }, logical(1)),

    stringsAsFactors = FALSE
  )
}

# ---------- Parse CARD reference models ----------

card_parse_models <- function(path) {
  message("Reading CARD card.json...")
  message("This step may require several GB of RAM.")

  raw <- jsonlite::fromJSON(
    path, simplifyVector = FALSE
  )

  valid <- vapply(raw, function(z) {
    is.list(z) &&
      !is.null(z$ARO_accession) &&
      !is.null(z$ARO_name)
  }, logical(1))

  models <- raw[valid]

  model_rows <- vector(
    "list", length(models)
  )

  sequence_rows <- vector(
    "list", length(models)
  )

  for (i in seq_along(models)) {
    z <- models[[i]]

    aro_id <- card_text(z$ARO_accession)

    if (is.na(aro_id)) next

    if (!grepl("^ARO:", aro_id)) {
      aro_id <- paste0("ARO:", aro_id)
    }

    mid <- card_text(z$model_id)
    mtype <- card_text(z$model_type)

    if (is.na(mtype)) {
      mtype <- card_text(z$model_type_id)
    }

    categories <- z$ARO_category

    if (is.null(categories)) {
      categories <- list()
    }

    category_values <- function(class_name) {
      terms <- vapply(categories, function(c) {
        if (!is.list(c)) return(NA_character_)

        kind <- card_text(
          c$category_aro_class_name
        )

        if (!is.na(kind) &&
            identical(kind, class_name)) {
          card_text(c$category_aro_name)
        } else {
          NA_character_
        }
      }, character(1))

      card_join(terms)
    }

    model_rows[[i]] <- data.frame(
      aro_accession = aro_id,
      model_id = mid,
      model_type = mtype,

      aro_name_json = card_text(
        z$ARO_name
      ),

      description_json = card_text(
        z$ARO_description
      ),

      drug_class_json = category_values(
        "Drug Class"
      ),

      mechanism_json = category_values(
        "Resistance Mechanism"
      ),

      family_json = category_values(
        "AMR Gene Family"
      ),

      stringsAsFactors = FALSE
    )

    seqs <- z$model_sequences$sequence

    if (is.null(seqs) || !length(seqs)) next

    rows <- vector(
      "list", length(seqs)
    )

    for (k in seq_along(seqs)) {
      s <- seqs[[k]]

      if (!is.list(s)) next

      dna <- card_text(
        s$dna_sequence$accession
      )

      protein <- card_text(
        s$protein_sequence$accession
      )

      if (is.na(dna) &&
          is.na(protein)) next

      rows[[k]] <- data.frame(
        aro_accession = aro_id,
        model_id = mid,
        model_type = mtype,
        dna_accession = dna,
        protein_accession = protein,

        taxon = card_text(
          s$NCBI_taxonomy$NCBI_taxonomy_name
        ),

        stringsAsFactors = FALSE
      )
    }

    rows <- Filter(Negate(is.null), rows)

    if (length(rows)) {
      sequence_rows[[i]] <- do.call(
        rbind, rows
      )
    }
  }

  model_rows <- Filter(
    Negate(is.null), model_rows
  )

  sequence_rows <- Filter(
    Negate(is.null), sequence_rows
  )

  if (!length(model_rows)) {
    stop("No CARD models were found.")
  }

  model_table <- unique(
    do.call(rbind, model_rows)
  )

  sequence_table <- if (length(sequence_rows)) {
    unique(do.call(rbind, sequence_rows))
  } else {
    data.frame(
      aro_accession = character(),
      model_id = character(),
      model_type = character(),
      dna_accession = character(),
      protein_accession = character(),
      taxon = character()
    )
  }

  rm(raw, models)
  gc()

  list(
    models = model_table,
    sequences = sequence_table
  )
}

# ---------- Collapse model annotations by ARO ----------

card_group_models <- function(models) {
  groups <- split(
    seq_len(nrow(models)),
    models$aro_accession
  )

  ans <- lapply(groups, function(idx) {
    x <- models[idx, , drop = FALSE]

    data.frame(
      aro_accession = x$aro_accession[1],

      aro_name_json = card_join(
        x$aro_name_json
      ),

      description_json = card_join(
        x$description_json
      ),

      json_drug_class = card_join(
        x$drug_class_json
      ),

      json_mechanism = card_join(
        x$mechanism_json
      ),

      json_family = card_join(
        x$family_json
      ),

      model_ids = card_join(
        x$model_id
      ),

      model_types = card_join(
        x$model_type
      ),

      number_of_models = length(unique(
        x$model_id[!is.na(x$model_id)]
      )),

      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, ans)
  rownames(out) <- NULL

  out
}

# ============================================================
# MAIN FUNCTION
# ============================================================

build_full_card_mapping <- function(
    aro_tsv = "data/aro.tsv",
    out_dir = "CARD_full_mapping",
    refresh = TRUE,
    excel = TRUE) {

  dir.create(
    out_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  cache_dir <- file.path(
    out_dir, "CARD_download_cache"
  )

  # 1. Read all input ARO terms

  original <- card_read_aro(aro_tsv)

  message(
    "ARO terms in input: ",
    nrow(original)
  )

  # 2. Download current CARD files

  data_archive <- card_download(
    "https://card.mcmaster.ca/latest/data",
    file.path(
      cache_dir, "CARD_latest_data.tar.bz2"
    ),
    refresh
  )

  ontology_archive <- card_download(
    "https://card.mcmaster.ca/latest/ontology",
    file.path(
      cache_dir, "CARD_latest_ontology.tar.bz2"
    ),
    refresh
  )

  # 3. Extract JSON and ontology

  json_file <- card_extract(
    data_archive, "card.json", cache_dir
  )

  obo_file <- card_extract(
    ontology_archive, "aro.obo", cache_dir
  )

  # 4. Parse CARD models

  parsed <- card_parse_models(json_file)

  model_catalog <- card_group_models(
    parsed$models
  )

  # 5. Parse ontology and derive annotations

  ontology <- card_parse_obo(obo_file)

  inferred <- card_ontology_categories(
    ontology,
    original$aro_accession
  )

  # 6. Join using ARO accession, not DNA accession

  x <- merge(
    original,
    model_catalog,
    by = "aro_accession",
    all.x = TRUE,
    sort = FALSE
  )

  x <- merge(
    x,
    inferred,
    by = "aro_accession",
    all.x = TRUE,
    sort = FALSE
  )

  x <- x[
    match(original$aro_accession,
          x$aro_accession),
    ,
    drop = FALSE
  ]

  if (nrow(x) != nrow(original) ||
      anyDuplicated(x$aro_accession)) {
    stop("ARO mapping changed the number of entries.")
  }

  # 7. Prepare comprehensive annotation columns

  x$`Source` <- NA_character_

  x$`ARO Accession` <- x$aro_accession

  x$`ARO Name` <- card_first(
    x$aro_name_json,
    x$aro_name_input
  )

  x$`Gene` <- x$`ARO Name`

  x$`Drug Class` <- card_first(
    x$json_drug_class,
    x$ontology_drug_class
  )

  x$`Resistance Mechanism` <- card_first(
    x$json_mechanism,
    x$ontology_mechanism
  )

  x$`AMR Gene Family` <- x$json_family

  x$`CARD Description` <- card_first(
    x$description_json,
    x$description_input
  )

  x$`CARD Short Name` <- x$card_short_name

  x$`Model IDs` <- x$model_ids

  x$`Model Types` <- x$model_types

  x$`Number of Models` <- ifelse(
    is.na(x$number_of_models),
    0L,
    x$number_of_models
  )

  x$`Has CARD Model` <-
    x$`Number of Models` > 0L

  x$`ARO Name Conflict` <-
    !is.na(x$aro_name_json) &
    !is.na(x$aro_name_input) &
    tolower(trimws(x$aro_name_json)) !=
    tolower(trimws(x$aro_name_input))

  x$`Drug Class Source` <- ifelse(
    !is.na(x$json_drug_class),
    "CARD card.json category",
    ifelse(
      !is.na(x$ontology_drug_class),
      "CARD ontology relation",
      NA_character_
    )
  )

  x$`Mechanism Source` <- ifelse(
    !is.na(x$json_mechanism),
    "CARD card.json category",
    ifelse(
      !is.na(x$ontology_mechanism),
      "CARD ontology relation",
      NA_character_
    )
  )

  x$`CARD Term URL` <- paste0(
    "https://card.mcmaster.ca/aro/",
    sub("^ARO:", "", x$aro_accession)
  )

  x$`CARD Release Downloaded` <-
    as.character(Sys.Date())

  # 8. Prepare reference sequence crosswalk

  seqs <- parsed$sequences

  seqs <- seqs[
    seqs$aro_accession %in% x$aro_accession,
    ,
    drop = FALSE
  ]

  seqs$dna_key <- card_norm_acc(
    seqs$dna_accession
  )

  seqs$protein_key <- card_norm_acc(
    seqs$protein_accession
  )

  dna_by_aro <- split(
    seqs$dna_accession,
    seqs$aro_accession
  )

  protein_by_aro <- split(
    seqs$protein_accession,
    seqs$aro_accession
  )

  x$`DNA Reference Accessions` <- vapply(
    x$aro_accession,
    function(id) card_join(dna_by_aro[[id]]),
    character(1)
  )

  x$`Protein Reference Accessions` <- vapply(
    x$aro_accession,
    function(id) card_join(protein_by_aro[[id]]),
    character(1)
  )

  x$`Number of DNA References` <- vapply(
    x$aro_accession,
    function(id) {
      v <- dna_by_aro[[id]]
      if (is.null(v)) return(0L)
      as.integer(length(unique(
        v[!is.na(v) & nzchar(v)]
      )))
    },
    integer(1)
  )
  # Rename the original ontology ID column
  x$`Ontology Numeric ID` <- x$ontology_numeric_id

  # 9. Full ARO catalog: one row per ARO accession

  # Ensure the ontology numeric ID has the correct name
  if (!"Ontology Numeric ID" %in% names(x)) {
    if ("ontology_numeric_id" %in% names(x)) {
      x[["Ontology Numeric ID"]] <- x[["ontology_numeric_id"]]
    }
  }

  # Columns required in the final catalog
  required_cols <- c(
    "Source",
    "Gene",
    "ARO Accession",
    "ARO Name",
    "Drug Class",
    "Resistance Mechanism",
    "AMR Gene Family",
    "CARD Description",
    "CARD Short Name",
    "Ontology Numeric ID",
    "Has CARD Model",
    "Number of Models",
    "Model IDs",
    "Model Types",
    "DNA Reference Accessions",
    "Protein Reference Accessions",
    "Number of DNA References",
    "Drug Class Source",
    "Mechanism Source",
    "ARO Name Conflict",
    "CARD Term URL",
    "CARD Release Downloaded"
  )

  # Identify missing columns before selecting
  missing_cols <- setdiff(required_cols, names(x))

  if (length(missing_cols) > 0) {
    message(
      "Missing columns: ",
      paste(missing_cols, collapse = ", ")
    )

    stop(
      "Cannot build CARD catalog: required columns missing."
    )
  }

  # Build the final catalog
  all_aro <- x[, required_cols, drop = FALSE]

  # Verify completeness
  stopifnot(nrow(all_aro) == nrow(original))

  message(
    "CARD catalog created successfully: ",
    nrow(all_aro),
    " ARO entries"
  )

  # 10. Sequence-level reference mapping

  crosswalk <- merge(
    seqs,
    all_aro[, c(
      "ARO Accession",
      "ARO Name",
      "Drug Class",
      "Resistance Mechanism",
      "AMR Gene Family"
    )],
    by.x = "aro_accession",
    by.y = "ARO Accession",
    all.x = TRUE,
    sort = FALSE
  )

  names(crosswalk)[
    names(crosswalk) == "dna_accession"
  ] <- "ID"

  crosswalk$`ARO Accession` <-
    crosswalk$aro_accession

  crosswalk$`Gene` <- crosswalk$`ARO Name`

  crosswalk$`Source` <- NA_character_

  crosswalk <- crosswalk[, c(
    "Source",
    "Gene",
    "ID",
    "ARO Accession",
    "ARO Name",
    "Drug Class",
    "Resistance Mechanism",
    "AMR Gene Family",
    "model_id",
    "model_type",
    "protein_accession",
    "taxon",
    "dna_key",
    "protein_key"
  ), drop = FALSE]

  # 11. Mapping table in original Excel format.
  # AROs without DNA accessions are retained.

  mapping <- all_aro[, c(
    "Source",
    "Gene",
    "ARO Accession",
    "ARO Name",
    "Drug Class",
    "Resistance Mechanism",
    "AMR Gene Family",
    "CARD Description",
    "CARD Short Name",
    "Has CARD Model"
  ), drop = FALSE]

  id_index <- split(
    crosswalk$ID,
    crosswalk$`ARO Accession`
  )

  repeated <- lapply(
    seq_len(nrow(mapping)),
    function(i) {
      ids <- unique(
        id_index[[mapping$`ARO Accession`[i]]]
      )

      ids <- sort(
        ids[!is.na(ids) & nzchar(ids)]
      )

      if (!length(ids)) {
        ids <- NA_character_
      }

      data.frame(
        index = i,
        ID = ids
      )
    }
  )

  repeated <- do.call(rbind, repeated)

  mapping <- cbind(
    mapping[repeated$index, , drop = FALSE],
    ID = repeated$ID
  )

  mapping <- mapping[, c(
    "Source",
    "Gene",
    "ID",
    "ARO Accession",
    "ARO Name",
    "Drug Class",
    "Resistance Mechanism",
    "AMR Gene Family",
    "CARD Description",
    "CARD Short Name",
    "Has CARD Model"
  ), drop = FALSE]

  rownames(mapping) <- NULL

  # 12. Identify DNA accessions linked to multiple AROs

  audit_seqs <- unique(
    seqs[
      !is.na(seqs$dna_key) &
        nzchar(seqs$dna_key),
      c("dna_key", "aro_accession")
    ]
  )

  grouped <- split(
    audit_seqs$aro_accession,
    audit_seqs$dna_key
  )

  ambiguous_list <- grouped[
    vapply(
      grouped,
      function(z) length(unique(z)) > 1L,
      logical(1)
    )
  ]

  ambiguous <- if (length(ambiguous_list)) {
    data.frame(
      accession_key = names(ambiguous_list),

      number_of_aro_terms = vapply(
        ambiguous_list,
        function(z) length(unique(z)),
        integer(1)
      ),

      aro_candidates = vapply(
        ambiguous_list,
        card_join,
        character(1)
      )
    )
  } else {
    data.frame(
      accession_key = character(),
      number_of_aro_terms = integer(),
      aro_candidates = character()
    )
  }

  # 13. Flag incomplete annotations

  review <- all_aro[
    all_aro$`ARO Name Conflict` |
      !all_aro$`Has CARD Model` |
      is.na(all_aro$`Drug Class`) |
      is.na(all_aro$`Resistance Mechanism`),
    ,
    drop = FALSE
  ]

  review$`Review Reason` <- ifelse(
    review$`ARO Name Conflict`,
    "Name differs between input TSV and CARD model",
    ifelse(
      !review$`Has CARD Model`,
      "No CARD detection model",
      "Missing drug class or mechanism"
    )
  )

  # 14. Save CSV files

  card_write(
    all_aro,
    file.path(out_dir, "CARD_all_ARO_terms.csv")
  )

  card_write(
    mapping,
    file.path(
      out_dir,
      "AMRgenelist_mapped_CARD_ALL.csv"
    )
  )

  card_write(
    crosswalk,
    file.path(
      out_dir,
      "CARD_reference_sequence_crosswalk.csv"
    )
  )

  card_write(
    parsed$models,
    file.path(
      out_dir,
      "CARD_model_details.csv"
    )
  )

  card_write(
    ambiguous,
    file.path(
      out_dir,
      "CARD_ambiguous_DNA_accessions.csv"
    )
  )

  card_write(
    review,
    file.path(
      out_dir,
      "CARD_annotation_review.csv"
    )
  )

  # 15. Save annotation coverage statistics

  summary <- c(
    paste0("CARD download date: ", Sys.Date()),
    "CARD data: https://card.mcmaster.ca/latest/data",
    "CARD ontology: https://card.mcmaster.ca/latest/ontology",
    paste0(
      "Data MD5: ",
      unname(tools::md5sum(data_archive))
    ),
    paste0(
      "Ontology MD5: ",
      unname(tools::md5sum(ontology_archive))
    ),
    paste0("Input ARO terms: ", nrow(original)),
    paste0("Output ARO terms: ", nrow(all_aro)),
    paste0(
      "Terms with CARD models: ",
      sum(all_aro$`Has CARD Model`)
    ),
    paste0(
      "Terms with drug class: ",
      sum(!is.na(all_aro$`Drug Class`))
    ),
    paste0(
      "Terms with resistance mechanism: ",
      sum(!is.na(all_aro$`Resistance Mechanism`))
    ),
    paste0(
      "Terms with DNA references: ",
      sum(all_aro$`Number of DNA References` > 0)
    ),
    paste0("Mapping rows: ", nrow(mapping)),
    paste0("Reference rows: ", nrow(crosswalk)),
    paste0("Ambiguous DNA accessions: ", nrow(ambiguous))
  )

  writeLines(
    summary,
    file.path(out_dir, "CARD_mapping_summary.txt")
  )

  # 16. Create Excel workbook

  if (excel) {
    if (requireNamespace("writexl", quietly = TRUE)) {

      workbook <- list(
        All_ARO_Terms = all_aro,
        Mapping_Style = mapping,
        Reference_Sequences = crosswalk,
        Ambiguous_Accessions = ambiguous,
        Annotation_Review = review
      )

      writexl::write_xlsx(
        workbook,
        file.path(
          out_dir,
          "AMRgenelist_mapped_CARD_ALL.xlsx"
        )
      )

      message("Excel workbook saved.")

    } else {
      message(
        "Install writexl to generate Excel output."
      )
    }
  }

  message(
    "Completed: ",
    nrow(all_aro),
    " ARO entries; ",
    nrow(mapping),
    " mapping rows."
  )

  invisible(list(
    catalog = all_aro,
    mapping = mapping,
    crosswalk = crosswalk,
    ambiguous = ambiguous,
    review = review,
    summary = summary
  ))
}

# ============================================================
# RUN
# ============================================================

result <- build_full_card_mapping(
  aro_tsv = "data/aro.tsv",
  out_dir = "CARD_full_mapping",
  refresh = FALSE,
  excel = TRUE
)

# Check output
cat("\nNumber of ARO terms:",
    nrow(result$catalog), "\n")

cat("Number of mapping rows:",
    nrow(result$mapping), "\n")

cat("Number of reference sequence records:",
    nrow(result$crosswalk), "\n")

cat("Number of accessions with multiple ARO matches:",
    nrow(result$ambiguous), "\n")

print(result$summary)


source("scripts/build_full_CARD_mapping_online.R")

result <- build_full_card_mapping(
  aro_tsv = "data/aro.tsv",
  out_dir = "CARD_full_mapping",
  refresh = FALSE,
  excel = TRUE
)

save.image("build_full_CARD_mapping_online.rdata")
