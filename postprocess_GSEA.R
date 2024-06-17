# Script to postprocess GSEA results to keep only those gene sets below a given
# FDR q cutoff and to make the output compatible with Brainarray Entrez Gene
# probeset mapping (instead of Affymetrix probeset mapping)
#
# Adam Gower

arglist <- commandArgs(trailing=TRUE)
if (length(arglist) != 3) {
  stop(
    "Usage: postprocess_GSEA.R --args ",
    "[run path] [description of comparison] [output filename]"
  )
} else {
  run.path <- arglist[1]
  comparison <- arglist[2]
  output.filename <- arglist[3]
}

q.cutoff <- 0.25

cat("Working on run folder:", run.path)
cat("\n")
setwd(run.path)

# Set phenotype labels based on comparison
phenotypes <- list(
  pos = paste("upregulated with respect to", comparison),
  neg = paste("downregulated with respect to", comparison)
)

# Rename the "na_pos"/"na_neg" fields in the index.html file
cat("Editing index.html.\n")
report <- readLines("index.html")
for (direction in c("pos", "neg")) {
  report <- sub(
    paste0(
      "na</b></h4><ul><li>",
      "([0-9]+ / [0-9]+|None of the) gene sets are ",
      "(upregulated|enriched) in phenotype <b>",
      "na_", direction
    ),
    paste0(
      phenotypes[[direction]],
      "</b></h4><ul><li>\\1 gene sets are <b>",
      phenotypes[[direction]]
    ),
    report
  )
}
report <- sub(
  "na_pos</b><i> versus </i><b>na_neg",
  paste0(phenotypes$pos, "</b><i> versus </i><b>", phenotypes$neg),
  report
)
write(report, "index.html")

for (direction in c("pos", "neg")) {
  # Get name of '.xls' (TSV) file containing output for given direction
  report_filename <- list.files(
    pattern=paste("gsea_report_for_na", direction, "[0-9]+.xls", sep="_")
  )
  # Read TSV file into a data frame
  report <- read.delim(
    report_filename, stringsAsFactors=FALSE, check.names=FALSE, row.names=1
  )

  # Identify the rows (gene sets) that pass the FDR q cutoff
  gene_sets_pass <- report[["FDR q-val"]] < q.cutoff

  # Iterate over each gene set passing FDR q cutoff
  for (gene.set in rownames(report)[gene_sets_pass]) {
    # Edit Details .html file
    filename <- paste0(gene.set, ".html")
    cat("Editing file:", filename)
    cat("\n")
    html <- readLines(filename)
    # Remove Phenotype row and rename "na_pos"/"na_neg" phenotype labels
    html <- sub(
      "<tr><td>Phenotype</td><td>NoPhenotypeAvailable</td></tr>", "", html
    )
    html <- sub(paste0("na_", direction), phenotypes[[direction]], html)
    # Replace incorrect Affymetrix link with Entrez Gene link
    # and remove unneeded 'Entrez' and 'Source' links
    html <- gsub(
      "https\\://www.affymetrix.com/LinkServlet\\?probeset=([0-9]+)",
      "http://www.ncbi.nlm.nih.gov/gene?term=\\1[uid]",
      html
    )
    html <- gsub(
      "<br><a href='[^']+'>Entrez</a>, &nbsp<a href='[^']+'>Source</a>",
      "",
      html
    )
    # Remove " [Source:...]" suffix from end of GENE_TITLE column
    html <- gsub(" \\[Source:[^]]+\\]", "", html)
    # Reconstitute HTML and write back to file
    write(paste(html, collapse="\n"), filename)

    # Edit .xls (TSV) file containing GSEA details table
    filename <- paste0(gene.set, ".xls")
    cat("Editing file:", filename)
    cat("\n")
    xls <- readLines(filename)
    # Remove " [Source:...]" suffix from end of GENE_TITLE column
    xls <- gsub(" \\[Source: [^\\]]+\\]", "", xls)
    # Write table back to file
    write(paste(xls, collapse="\n"), filename)
  }

  # Remove the .html, .xls, and .png files that correspond to the gene sets
  # that fail the FDR q cutoff
  for (gene.set in rownames(report)[!gene_sets_pass]) {
    cat("Removing Details files for gene set:", gene.set)
    cat("\n")
    unlink(paste0(gene.set, ".html"), force=TRUE)
    unlink(paste0(gene.set, ".xls"), force=TRUE)
    filename <- list.files(
      pattern=paste("enplot", gene.set, "[0-9]+.png", sep="_")
    )
    if (length(filename) == 1) {
      id <- as.integer(sub("^enplot_.+_([0-9]+).png$", "\\1", filename))
      unlink(paste0("enplot_", gene.set, "_", id, ".png"), force=TRUE)
      unlink(paste0("gset_rnd_es_dist_", id+1, ".png"), force=TRUE)
    }
  }

  # Lightly parse the "snapshot" .html file to extract header, table elements,
  # and footer, removing row delimiters
  filename <- paste(direction, "snapshot.html", sep="_")
  cat("Editing:", filename)
  cat("\n")
  snapshot <- readLines(filename)
  snapshot[[2]] <- gsub("</?tr>", "", snapshot[[2]])
  snapshot[2] <- strsplit(snapshot[2], "</?td>(<td>)?")
  # Keep only table elements corresponding to gene sets that passed FDR q cutoff
  i <- which(gene_sets_pass)
  n <- sum(gene_sets_pass)
  snapshot[[2]] <- snapshot[[2]][c(1, i+1, length(snapshot[[2]]))]
  # Fix the title of the HTML file
  snapshot[[2]][1] <-
    sub(
      "Snapshot of [0-9]+ enrichment plots",
      paste("Snapshot of", n, "enrichment plots"),
      snapshot[[2]][1]
    )
  # Reconstitute HTML and write back to file
  snapshot[[2]][1+(1:n)] <- paste0("<td>", snapshot[[2]][1+(1:n)], "</td>")
  snapshot[[2]] <- paste(
    c(
      snapshot[[2]][1],
      na.omit(
        c(
          rbind(
            "<tr>",
            matrix(
              c(snapshot[[2]][(1:n)+1], rep(NA, 3*ceiling(n/3)-n)),
              nrow=3, ncol=ceiling(n/3)
            ),
            "</tr>"
          )
        )
      ),
      snapshot[[2]][n+2]
    ),
    collapse=""
  )
  snapshot <- paste(snapshot, collapse="\n")
  write(snapshot, filename)

  # Lightly parse report .html file to extract header, table rows, and footer
  filename <- list.files(
    pattern=paste("gsea_report_for_na", direction, "[0-9]+.html", sep="_")
  )
  cat("Editing:", filename)
  cat("\n")
  html <- readLines(filename)
  html[2] <- strsplit(html[2], "</?tr>(<tr>)?")
  # Fix the title and caption of the HTML file
  html[[2]][1] <- sub(
    "Report for na_pos [0-9]+",
    paste("Report for", phenotypes[[direction]]),
    html[[2]][1]
  )
  html[[2]][length(html[[2]])] <- sub(
    "phenotype <b>na<b>",
    paste0("phenotype <b>", phenotypes[[direction]], "<b>"),
    html[[2]][length(html[[2]])]
  )
  # Remove hyperlinks from table rows that failed FDR q cutoff
  html[[2]][which(report[["FDR q-val"]] >= q.cutoff)+1] <- sub(
    pattern = paste0(
      "<td><a href='[^']+'>([^<]+)</a></td><td><a href='[^']+'>",
      "Details ...",
      "</a></td>"
    ),
    replacement = "<td>\\1</td><td></td>",
    x = html[[2]][which(report[["FDR q-val"]] >= q.cutoff)+1]
  )
  # Reconstitute HTML and write back to file
  html[[2]][2:(length(html[[2]])-1)] <- paste0(
    "<tr>", html[[2]][2:(length(html[[2]])-1)], "</tr>"
  )
  html[[2]] <- paste(html[[2]], collapse="")
  html <- paste(html, collapse="\n")
  write(html, filename)
}

# Combine '.xls' (TSV) output files into single TSV file for import to Excel
output <- NULL
for (direction in c("neg", "pos")) {
  filename <- list.files(
    pattern=paste("gsea_report_for_na", direction, "[0-9]+.xls", sep="_")
  )
  output <- rbind(
    output, read.delim(filename, check.names=FALSE, stringsAsFactors=FALSE)
  )
}
# Rename/remove columns
colnames(output)[colnames(output) == "NAME"] <- "Gene Set Name"
output[["GS<br> follow link to MSigDB"]] <- NULL
output[["GS DETAILS"]] <- NULL
colnames(output)[colnames(output) == "SIZE"] <- "Gene Set Size"
colnames(output)[colnames(output) == "ES"] <- "Enrichment Score (ES)"
colnames(output)[colnames(output) == "NES"] <-
  "Normalized Enrichment Score (NES)"
colnames(output)[colnames(output) == "NOM p-val"] <- "Nominal p value"
colnames(output)[colnames(output) == "FDR q-val"] <- "FDR q value"
output[["FWER p-val"]] <- NULL
output[["RANK AT MAX"]] <- NULL
output[["LEADING EDGE"]] <- NULL
# Remove column with empty column name
output[colnames(output) == ""] <- list(NULL)

# Extract the list of gmt filenames from the .rpt file
rpt.filename <- list.files(pattern=".+.GseaPreranked.[0-9]+.rpt")
rpt <- readLines(rpt.filename)
params <- strsplit(rpt[grep("^param", rpt)], "\t")
params <- structure(sapply(params, "[", 3), names=sapply(params, "[", 2))
gmt.filenames <- strsplit(params["gmx"], ",")[[1]]
# Define human-readable labels for the various MSigDB collections
group.labels <- c(
  "h.all"              = "H Hallmark",
  "c1.all"             = "C1 Cytogenetic bands",
  "c2.all"             = "C2 Curated gene sets",
  "c2.cgp"             = "C2 Chemical/genetic perturbations",
  "c2.cp"              = "C2 Canonical pathways",
  "c2.cp.biocarta"     = "C2 BioCarta pathway",
  "c2.cp.kegg"         = "C2 KEGG pathway",
  "c2.cp.kegg_legacy"  = "C2 KEGG legacy pathway",
  "c2.cp.kegg_medicus" = "C2 KEGG MEDICUS pathway",
  "c2.cp.pid"          = "C2 PID pathway",
  "c2.cp.reactome"     = "C2 Reactome pathway",
  "c2.cp.wikipathways" = "C2 WikiPathways pathway",
  "c3.all"             = "C3 Regulatory targets",
  "c3.mir"             = "C3 microRNA motifs",
  "c3.mir.mir_legacy"  = "C3 Legacy microRNA motif",
  "c3.mir.mirdb"       = "C3 miRDB microRNA motif",
  "c3.tft"             = "C3 TF motifs",
  "c3.tft.tft_legacy"  = "C3 Legacy TF motif",
  "c3.tft.gtrd"        = "C3 GTRD TF motif",
  "c4.all"             = "C4 Computational gene sets",
  "c4.3ca"             = "C4 Curated Cancer Cell Atlas",
  "c4.cgn"             = "C4 Cancer gene neighborhoods",
  "c4.cm"              = "C4 Cancer modules",
  "c5.all"             = "C5 Ontology gene sets",
  "c5.go"              = "C5 Gene Ontology",
  "c5.bp"              = "C5 GO Biological Process",
  "c5.go.bp"           = "C5 GO Biological Process",
  "c5.cc"              = "C5 GO Cellular Component",
  "c5.go.cc"           = "C5 GO Cellular Component",
  "c5.mf"              = "C5 GO Molecular Function",
  "c5.go.mf"           = "C5 GO Molecular Function",
  "c5.hpo"             = "C5 Human Phenotype Ontology",
  "c6.all"             = "C6 Oncogenic signatures",
  "c7.all"             = "C7 Immunologic signatures",
  "c7.immunesigdb"     = "C7 ImmuneSigDB",
  "c7.vax"             = "C7 HIPC vaccine response",
  "c8.all"             = "C8 Cell type signatures",
  "mh.all"             = "MH Hallmark",
  "m1.all"             = "M1 Cytogenetic bands",
  "m2.all"             = "M2 Curated gene sets",
  "m2.cgp"             = "M2 Chemical/genetic perturbations",
  "m2.cp"              = "M2 Canonical pathways",
  "m2.cp.biocarta"     = "M2 BioCarta pathway",
  "m2.cp.reactome"     = "M2 Reactome pathway",
  "m2.cp.wikipathways" = "M2 WikiPathways pathway",
  "m3.all"             = "M3 Regulatory targets",
  "m3.mirdb"           = "M3 miRDB microRNA motif",
  "m3.gtrd"            = "M3 GTRD TF motif",
  "m5.all"             = "M5 Ontology gene sets",
  "m5.go"              = "M5 Gene Ontology",
  "m5.go.bp"           = "M5 GO Biological Process",
  "m5.go.cc"           = "M5 GO Cellular Component",
  "m5.go.mf"           = "M5 GO Molecular Function",
  "m5.mpt"             = "M5 Mouse Phenotype Ontology MP Tumor",
  "m8.all"             = "M8 Cell type signatures"
)
# Add 'Group' column (MSigDB collection name) to beginning of output data frame
group.table <- NULL
for (filename in gmt.filenames) {
  gene.set.names <- sapply(strsplit(readLines(filename), "\t"), "[", 1)
  group.label <- group.labels[
    sub("\\.v.+\\.entrez\\.gmt$", "", basename(filename))
  ]
  group.table <- rbind(
    group.table,
    data.frame(
      Group=group.label, Name=gene.set.names,
      row.names=NULL, stringsAsFactors=FALSE
    )
  )
}
output <- cbind(
  Group=group.table$Group[match(output[["Gene Set Name"]], group.table$Name)],
  output
)

# Convert the 'Gene Set Name' column to an Excel HYPERLINK formula
output[["Gene Set Name"]] <- paste0(
  "=HYPERLINK(",
  dQuote(
    file.path(
      "http://www.broadinstitute.org/gsea/msigdb/cards",
      paste0(output[["Gene Set Name"]], ".html")
    )
  ),
  ",",
  dQuote(output[["Gene Set Name"]]),
  ")"
)
# Reorder by NES (in descending order by default)
output <- output[
  order(output[["Normalized Enrichment Score (NES)"]], decreasing=TRUE),
]
# Write output to a tab-delimited file
write.table(
  output, output.filename,
  quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE
)
