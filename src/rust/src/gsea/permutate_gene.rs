use std::collections::HashSet;

use rand::seq::SliceRandom;
use rand_chacha::{rand_core::SeedableRng, ChaCha8Rng};
use rayon::prelude::*;

use super::algorithm::*;

pub(super) fn gsea_gene(
    identifiers: &[&str],
    metrics: &[f64],
    genesets: &[HashSet<&str>],
    exponent: f64,
    nperm: usize,
    seed: usize,
) -> GSEAOutput {
    let (_, ids, weights) = gsea_par_prerank(identifiers, metrics, exponent);
    let mut output = GSEAOutput::new(genesets.len());
    for geneset in genesets {
        let hits = gsea_par_hits(geneset, &ids);
        if hits.iter().any(|hit| *hit) {
            let norm_neg = gsea_norm_neg(&hits);
            let norm_pos = gsea_norm_pos(&weights, &hits);
            let (running, es_pos, es) =
                gsea_running_es(&weights, &hits, norm_pos, norm_neg);
            let null = permutate_hits(
                &weights,
                hits,
                es >= 0.0,
                norm_neg,
                nperm,
                seed as u64,
            );
            output.add_output(running, es_pos, es, &null);
        } else {
            // If no gene found
            output.add_empty()
        }
    }
    return output;
}

fn permutate_hits(
    weights: &[f64],
    hits: Vec<bool>,
    direction: bool,
    norm_neg: f64, // won't be changed when shuffling hits
    nperm: usize,
    seed: u64,
) -> Vec<f64> {
    (0 .. nperm)
        .into_par_iter()
        .filter_map(|i| {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            rng.set_stream(i as u64);
            let mut shuffled = hits.clone();
            shuffled.shuffle(&mut rng);
            let norm_pos = gsea_norm_pos(weights, &shuffled);
            let es = gsea_es(weights, &shuffled, norm_pos, norm_neg);
            // Only keep the Enrichment scores with the same sign
            if direction == (es >= 0.0) {
                Some(es)
            } else {
                None
            }
        })
        .collect()
}
