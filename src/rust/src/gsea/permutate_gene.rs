use std::collections::HashSet;

use rand::prelude::IteratorRandom;
use rand_chacha::{rand_core::SeedableRng, ChaCha8Rng};
use rayon::prelude::*;

use super::algorithm::*;

pub(super) fn gsea_gene(
    identifiers: &[&str],
    metrics: &[f64],
    genesets: &Vec<HashSet<&str>>,
    exponent: f64,
    nperm: usize,
    seed: u64,
) -> GSEAOutput {
    let (_, ids, weights) = gsea_prerank(identifiers, metrics, exponent);
    let permutations: GSEAPermutate = genesets
        .par_iter()
        .map(|geneset| {
            let hits = gsea_hits(geneset, &ids);
            let input = GSEAInput::new(&weights, &hits);
            if input.hits.len() > 0 {
                let norm_neg = input.norm_neg(None);
                let score = input.score(None, None, Some(norm_neg));
                let null = permutate_hits(&input, norm_neg, nperm, seed);
                Some(GSEARunning { score, null })
            } else {
                None
            }
        })
        .collect();
    permutations.normalize_and_test()
}

fn permutate_hits(
    gsea_input: &GSEAInput,
    norm_neg: f64, // won't be changed when shuffling hits
    nperm: usize,
    seed: u64,
) -> Vec<f64> {
    let range = 0 .. gsea_input.weights.len();
    let nhits = gsea_input.hits.len();
    (0 .. nperm)
        .into_par_iter()
        .map(|i| {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            rng.set_stream(i as u64);
            let sample = range.clone().choose_multiple(&mut rng, nhits);
            gsea_input.es(Some(&sample), None, Some(norm_neg))
        })
        .collect()
}
