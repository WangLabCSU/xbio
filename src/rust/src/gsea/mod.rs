use std::collections::HashSet;

use extendr_api::prelude::*;

// Module used to do the actual work, calcualte the enrichment scores
mod algorithm;
mod input;
mod permutate_gene;
mod permutate_sample;

#[extendr]
// GSEA gene permutation analysis.
//
// # Arguments
// - `identifiers`: Character vector of gene identifiers.
// - `metrics`: Numeric vector of ranking metrics.
// - `genesets`: List of character vectors, each representing a geneset.
// - `exponent`: Weighting exponent for metrics.
// - `nperm`: Number of permutations.
//
// # Returns
// A list containing enrichment results.
fn gsea_gene_permutate(
    identifiers: Robj,
    metrics: Robj,
    genesets: Robj,
    exponent: f64,
    nperm: usize,
    threads: usize,
    seed: usize,
) -> std::result::Result<algorithm::GSEAOutput, String> {
    //  Check and parse `identifiers`
    let identifiers: Vec<&str> = identifiers
        .as_str_vector()
        .ok_or("`identifiers` must be a character")?;

    //  Check and parse `metrics`
    let metrics = metrics
        .as_real_slice()
        .ok_or("`metrics` must be a numeric")?;

    // Check and parse `genesets`
    let geneset_list = genesets.as_list().ok_or("`genesets` must be a list")?;
    let mut input_gs: Vec<HashSet<&str>> =
        Vec::with_capacity(geneset_list.len());
    for geneset in geneset_list.as_slice() {
        let gs = geneset
            .as_str_vector()
            .ok_or("`genesets` must be a list of character")?
            .into_iter()
            .collect::<HashSet<&str>>();
        input_gs.push(gs);
    }

    // Set rayon threads
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| e.to_string())?;

    // Run GSEA and Permutation
    let out = pool.install(|| {
        permutate_gene::gsea_gene(
            &identifiers,
            metrics,
            &input_gs,
            exponent,
            nperm,
            seed,
        )
    });

    Ok(out)
}

extendr_module! {
    mod gsea;
    fn gsea_gene_permutate;
}
