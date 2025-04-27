use rand::seq::SliceRandom;
use rand_chacha::{rand_core::SeedableRng, ChaCha8Rng};
use rayon::prelude::*;
use std::collections::HashSet;

use extendr_api::prelude::{list, List};

// we should provide a method to convert this object into Robj
pub(super) struct GSEAInput<'a> {
    pub(super) identifiers: Vec<&'a str>,
    pub(super) metrics: Vec<f64>,
    pub(super) genesets: Vec<HashSet<&'a str>>,
    pub(super) exponent: f64,
    pub(super) nperm: usize,
    // pub(super) threads: usize,
    pub(super) seed: usize,
}

impl<'a> GSEAInput<'a> {
    fn prerank(&self) -> (Vec<usize>, Vec<f64>, Vec<Vec<bool>>) {
        let metrics = &self.metrics;
        // Create a list of indices
        let mut indices: Vec<usize> = (0..metrics.len()).collect();

        // Sort indices by descending metrics
        indices.sort_unstable_by(|&a, &b| f64::total_cmp(&metrics[b], &self.metrics[a]));

        let mut ranking: Vec<f64> = metrics
            .par_iter()
            .map(|x| x.abs().powf(self.exponent))
            .collect();

        let ids: Vec<_> = indices.par_iter().map(|&i| &self.identifiers[i]).collect();
        ranking = indices.par_iter().map(|&i| ranking[i]).collect();

        let hit_list = self
            .genesets
            .par_iter()
            .map(|gs| ids.iter().map(|id| gs.contains(*id)).collect::<Vec<bool>>())
            .collect();
        (indices, ranking, hit_list)
    }

    pub fn gene_permutate(&self) -> List {
        let (indices, ranking, hit_list) = self.prerank();
        let mut es_list: Vec<f64> = Vec::with_capacity(hit_list.len());
        let mut esindex_list: Vec<usize> = Vec::with_capacity(hit_list.len());
        let mut nes_list: Vec<f64> = Vec::with_capacity(hit_list.len());
        let mut pvalue_list: Vec<f64> = Vec::with_capacity(hit_list.len());
        let mut running_es_list: Vec<Vec<f64>> = Vec::with_capacity(hit_list.len());
        for hits in hit_list {
            let norm_neg = 1.0 / ((ranking.len() - hits.iter().filter(|hit| **hit).count()) as f64);
            let running_escores = running_es(&ranking, &hits, &norm_neg);
            let (esindex, escore) = es_and_index(&running_escores);
            esindex_list.push(esindex);
            es_list.push(escore);
            running_es_list.push(running_escores);
            // Do the permutations
            let es_perm_list = permutate_hits(
                &escore,
                &self.nperm,
                &ranking,
                &norm_neg,
                &hits,
                self.seed as u64,
            );
            // Compute the P-value
            pvalue_list.push(compute_pval(&escore, &es_perm_list));
            // rescale by the permutation mean
            let mean = es_perm_list.iter().sum::<f64>() / (es_perm_list.len() as f64);
            // Compute NES (Normalized ES)
            nes_list.push(escore / mean);
        }
        list!(
            indices = indices,
            es = es_list,
            es_index = esindex_list,
            nes = nes_list,
            pvalue = pvalue_list,
            running_es = List::from_values(running_es_list)
        )
    }
}

fn permutate_hits(
    es_obs: &f64,
    nperm: &usize,
    ranking: &[f64],
    norm_neg: &f64,
    hits: &[bool],
    seed: u64,
) -> Vec<f64> {
    let hits_vec = hits.to_vec();
    let direction = es_obs >= &0.0;
    (0..*nperm)
        .into_par_iter()
        .filter_map(|i| {
            let mut rng = ChaCha8Rng::seed_from_u64(seed);
            rng.set_stream(i as u64);

            let mut shuffled = hits_vec.clone();
            shuffled.shuffle(&mut rng);

            let score = es(ranking, &shuffled, norm_neg);
            if direction == (score >= 0.0) {
                Some(score)
            } else {
                None
            }
        })
        .collect()
}

fn running_es(ranking: &[f64], hits: &[bool], norm_neg: &f64) -> Vec<f64> {
    let sum_pos: f64 = hits
        .iter()
        .zip(ranking)
        .filter_map(|(hit, metric)| hit.then_some(metric))
        .sum();
    let norm_pos = 1.0 / sum_pos;
    hits.iter()
        .zip(ranking)
        .map(|(hit, metric)| if *hit { metric * norm_pos } else { -norm_neg })
        .scan(0.0, |acc, score| {
            *acc += score;
            Some(*acc)
        })
        .collect()
}

fn es_and_index(running: &[f64]) -> (usize, f64) {
    let out = running.iter().enumerate().fold((0usize, &0.0f64), |x, y| {
        if x.1.abs() > y.1.abs() {
            x
        } else {
            y
        }
    });
    (out.0, out.1.to_owned())
}

fn es(ranking: &[f64], hits: &[bool], norm_neg: &f64) -> f64 {
    let sum_pos: f64 = hits
        .iter()
        .zip(ranking)
        .filter_map(|(hit, metric)| hit.then_some(metric))
        .sum();
    let norm_pos = 1.0 / sum_pos;
    hits.iter()
        .zip(ranking)
        .map(|(hit, metric)| if *hit { metric * norm_pos } else { -norm_neg })
        .fold(0.0f64, |old, score| -> f64 {
            let new = old + score;
            if old.abs() > new.abs() {
                old
            } else {
                new
            }
        })
}

fn compute_pval(es_obs: &f64, null: &[f64]) -> f64 {
    let total: f64 = null.len() as f64;
    let n: f64;
    if *es_obs >= 0.0 {
        n = null.iter().filter(|score| *score >= es_obs).count() as f64;
    } else {
        n = null.iter().filter(|score| *score <= es_obs).count() as f64;
    }
    n / total
}
