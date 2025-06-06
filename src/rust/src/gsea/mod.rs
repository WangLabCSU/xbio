use std::collections::HashSet;

use algorithm::GSEAOutput;
use extendr_api::prelude::*;

// Module used to do the actual work, calcualte the enrichment scores
mod algorithm;
mod input;
mod permutate_gene;
mod permutate_sample;
mod ssgsea;

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
) -> std::result::Result<GSEAOutput, String> {
    if nperm.le(&0usize) {
        return Err("`nperm` must be a positive integer".to_string());
    }

    //  Check and parse `identifiers`
    let identifiers: Vec<&str> = identifiers
        .as_str_vector()
        .ok_or("`identifiers` must be a character")?;

    //  Check and parse `metrics`
    let metrics = metrics
        .as_real_slice()
        .ok_or("`metrics` must be a numeric")?;

    // Check and parse `geneset_list`
    let input_gs = genesets.as_list().ok_or("`genesets` must be a list")?;
    let mut geneset_list: Vec<HashSet<&str>> =
        Vec::with_capacity(input_gs.len());
    for geneset in input_gs.as_slice() {
        let gs = geneset
            .as_str_vector()
            .ok_or("`genesets` must be a list of character")?
            .into_iter()
            .collect::<HashSet<&str>>();
        geneset_list.push(gs);
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
            &geneset_list,
            exponent,
            nperm,
            seed as u64,
        )
    });

    Ok(out)
}

// fn gsea_ssgsea(
//     array: Robj,
//     genesets: Robj,
//     exponent: f64,
//     nperm: usize,
//     threads: usize,
//     seed: usize,
// ) {
// }

extendr_module! {
    mod gsea;
    fn gsea_gene_permutate;
}
