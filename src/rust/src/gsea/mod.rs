use extendr_api::prelude::*;
use std::collections::HashSet;

// Module used to do the actual work, calcualte the enrichment scores
mod algorithm;

use algorithm::GSEAInput;

#[extendr]
/// GSEA gene permutation analysis.
///
/// # Arguments
/// - `identifiers`: Character vector of gene identifiers.
/// - `metrics`: Numeric vector of ranking metrics.
/// - `genesets`: List of character vectors, each representing a geneset.
/// - `exponent`: Weighting exponent for metrics.
/// - `nperm`: Number of permutations.
///
/// # Returns
/// A list containing enrichment results.
fn gsea_gene_permutate(
    identifiers: Robj,
    metrics: Vec<f64>,
    genesets: Robj,
    exponent: f64,
    nperm: usize,
    // threads: usize,
    // seed: usize,
) -> Result<List> {
    //  Check and parse `identifiers`
    let identifiers: Vec<&str> = identifiers
        .as_str_vector()
        .ok_or(Error::ExpectedRstr(r!("`identifiers` must be a character")))?;

    // Check and parse `genesets`
    let geneset_list = genesets
        .as_list()
        .ok_or(Error::ExpectedList(r!("`genesets` must be a list")))?;
    let mut input_gs: Vec<HashSet<&str>> = Vec::with_capacity(geneset_list.len());
    for geneset in geneset_list.as_slice() {
        let gs = geneset
            .as_str_vector()
            .ok_or(Error::ExpectedRstr(r!(
                "`genesets` must be a list of character"
            )))?
            .into_iter()
            .collect::<HashSet<&str>>();
        input_gs.push(gs);
    }

    // Run permutation and return result
    let input = GSEAInput {
        identifiers,
        metrics,
        genesets: input_gs,
        exponent,
        nperm,
    };
    Ok(input.gene_permutate())
}

extendr_module! {
    mod gsea;
    fn gsea_gene_permutate;
}
