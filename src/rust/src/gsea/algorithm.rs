use core::f64;
use std::{collections::HashSet, ops::Div};

use extendr_api::{list, List, Robj};
use rayon::prelude::*;

// For each geneset
pub struct GSEAInput<'a> {
    pub weights: &'a [f64],
    pub hits: &'a [usize],
}

impl<'a, 'b> GSEAInput<'a> {
    pub fn new(weights: &'a [f64], hits: &'a [usize]) -> Self {
        Self { weights, hits }
    }

    #[inline]
    pub fn switch_hits(&'b self, hits: Option<&'b [usize]>) -> &'b [usize] {
        hits.map_or_else(|| self.hits, |o| o)
    }

    #[inline]
    pub fn norm_neg(&self, hits: Option<&[usize]>) -> f64 {
        1.0 / ((self.weights.len() - self.switch_hits(hits).len()) as f64)
    }

    #[inline]
    pub fn norm_pos(&self, hits: Option<&[usize]>) -> f64 {
        let sum_pos: f64 = self
            .switch_hits(hits)
            .iter()
            .map(|hit| self.weights[*hit])
            .sum();
        1.0 / sum_pos
    }

    #[inline]
    pub fn score(
        &self,
        hits: Option<&[usize]>,
        norm_pos: Option<f64>,
        norm_neg: Option<f64>,
    ) -> GSEAEnrichmentScore {
        let norm_pos = norm_pos.map_or_else(|| self.norm_pos(hits), |s| s);
        let norm_neg = norm_neg.map_or_else(|| self.norm_neg(hits), |s| s);
        let hits = self.switch_hits(hits);
        let mut steps = vec![-norm_neg; self.weights.len()];
        for &i in hits {
            debug_assert!(i < self.weights.len());
            steps[i] = unsafe { self.weights.get_unchecked(i) * norm_pos }
        }
        let mut running = Vec::with_capacity(self.weights.len());
        let mut score: f64 = 0.0;
        let mut es: f64 = 0.0;
        let mut es_pos: usize = 0;
        for (i, delta) in steps.iter().enumerate() {
            score += delta;
            running.push(score);
            if score.abs() > es.abs() {
                es = score;
                es_pos = i
            }
        }
        GSEAEnrichmentScore {
            running,
            es_pos,
            es,
        }
    }

    // Used to create null distribution
    #[inline]
    pub fn es(
        &self,
        hits: Option<&[usize]>,
        norm_pos: Option<f64>,
        norm_neg: Option<f64>,
    ) -> f64 {
        let norm_pos = norm_pos.map_or_else(|| self.norm_pos(hits), |s| s);
        let norm_neg = norm_neg.map_or_else(|| self.norm_neg(hits), |s| s);
        let hits = self.switch_hits(hits);
        let mut nhits = hits.len();
        let mut hit_flags = vec![false; self.weights.len()];
        for &i in hits {
            debug_assert!(i < self.weights.len());
            hit_flags[i] = true;
        }
        let mut es = 0.0f64;
        let mut score = 0.0f64;
        for (i, hit_flag) in hit_flags.iter().enumerate() {
            if *hit_flag {
                nhits -= 1;
                // SAFETY: i is guaranteed in-bounds by `hit_flags`.
                score += unsafe { self.weights.get_unchecked(i) * norm_pos };
                if score.abs() > es.abs() {
                    es = score;
                }
                // early exit when no hit remains
                // which means all follows will move `es` to zero
                if nhits == 0 {
                    break;
                }
            } else {
                score += -norm_neg;
                if score.abs() > es.abs() {
                    es = score;
                }
            };
        }
        es
    }
}

pub struct GSEAEnrichmentScore {
    pub running: Vec<f64>,
    pub es_pos: usize,
    pub es: f64,
}

// Intermediate result combine the score and null distributions
// Used by GSEAPermutate to do normalization and statistical test
pub struct GSEARunning {
    pub score: GSEAEnrichmentScore,
    pub null: Vec<f64>,
}

impl GSEARunning {
    fn pos_null(&self) -> Vec<&f64> {
        self.null.iter().filter(|es| **es >= 0.0f64).collect()
    }
    fn neg_null(&self) -> Vec<&f64> {
        self.null.iter().filter(|es| **es <= 0.0f64).collect()
    }
    fn split_null(&self) -> (Vec<&f64>, Vec<&f64>) {
        let pos_null: Vec<&f64> = self.pos_null();
        let neg_null: Vec<&f64> = self.neg_null();
        (pos_null, neg_null)
    }

    fn normalize_null(&self, pos_mean: f64, neg_mean: f64) -> Vec<f64> {
        self.null
            .par_iter()
            .map(|es| {
                if es.ge(&0.0) {
                    es.div(pos_mean.abs())
                } else {
                    es.div(neg_mean.abs())
                }
            })
            .collect()
    }

    fn mean_null(null: &[&f64]) -> f64 {
        null.iter().copied().sum::<f64>() / (null.len() as f64)
    }

    fn test(&self, pos_null: Vec<&f64>, neg_null: Vec<&f64>) -> f64 {
        let total: f64;
        let n: f64;
        if self.score.es.ge(&0.0) {
            total = pos_null.len() as f64;
            n = pos_null
                .iter()
                .filter(|score| ***score >= self.score.es)
                .count() as f64;
        } else {
            total = neg_null.len() as f64;
            n = neg_null
                .iter()
                .filter(|score| ***score <= self.score.es)
                .count() as f64;
        }
        // https://stats.stackexchange.com/questions/299449/is-empirical-fdr-denominated-by-measurements-or-by-experiments
        // to be more precise, the permutation p-value
        // see Phipson & Smyth, 2010
        (n + 1.0f64) / (total + 1.0f64)
    }

    // Get nes and pvalue
    #[allow(unused)]
    pub fn normalize(&self) -> f64 {
        let reference;
        if self.score.es.ge(&0.0) {
            reference = self.pos_null();
        } else {
            reference = self.neg_null();
        }
        self.score.es.div(Self::mean_null(&reference).abs())
    }

    pub fn normalize_and_test(&self) -> (f64, (f64, Vec<f64>)) {
        let (pos_null, neg_null) = self.split_null();
        let pos_mean = Self::mean_null(&pos_null);
        let neg_mean = Self::mean_null(&neg_null);
        // normalize enrichment score
        let nes;
        if self.score.es.ge(&0.0) {
            nes = self.score.es.div(pos_mean.abs());
        } else {
            nes = self.score.es.div(neg_mean.abs());
        }
        let pvalue = self.test(pos_null, neg_null);
        let null_nes = self.normalize_null(pos_mean, neg_mean);
        (nes, (pvalue, null_nes))
    }
}

pub fn gsea_prerank<'a>(
    identifiers: &[&'a str],
    metrics: &[f64],
    exponent: f64,
) -> (Vec<usize>, Vec<&'a str>, Vec<f64>) {
    // Create a list of indices
    let mut indices: Vec<usize> = (0 .. metrics.len()).collect();

    // Sort indices by descending metrics
    indices.par_sort_unstable_by(|a, b| {
        f64::total_cmp(&metrics[*b], &metrics[*a])
    });

    // Define the weights
    let weights: Vec<f64> =
        metrics.par_iter().map(|x| x.abs().powf(exponent)).collect();

    // re-order the identifiers and weights
    let identifiers = indices
        .par_iter()
        .map(|&i| identifiers[i])
        .collect::<Vec<&str>>();
    let weights = indices.par_iter().map(|&i| weights[i]).collect();
    (indices, identifiers, weights)
}

pub fn gsea_hits(geneset: &HashSet<&str>, identifiers: &[&str]) -> Vec<usize> {
    identifiers
        .par_iter()
        .enumerate()
        .filter_map(|(i, id)| if geneset.contains(id) { Some(i) } else { None })
        .collect::<Vec<usize>>()
}

pub struct GSEAOutput {
    pub running_es: Vec<Vec<f64>>,
    pub es: Vec<f64>,
    pub es_pos: Vec<usize>,
    pub nes: Vec<f64>,
    pub pvalue: Vec<f64>,
    // pub null: Vec<Vec<f64>>,
    // pub mean: Vec<f64>,
    pub fdr: Vec<f64>,
}

impl From<GSEAOutput> for Robj {
    fn from(value: GSEAOutput) -> Self {
        let running_es_list: List = value.running_es.into_iter().collect();
        // let null_list: List = value.null.into_iter().collect();
        list!(
            running_es = running_es_list,
            es = value.es,
            es_pos = value.es_pos,
            nes = value.nes,
            pvalue = value.pvalue,
            // null_list = null_list,
            // mean = value.mean,
            fdr = value.fdr,
        )
        .into()
    }
}

pub struct GSEAPermutate(Vec<Option<GSEARunning>>);

impl FromIterator<Option<GSEARunning>> for GSEAPermutate {
    fn from_iter<I: IntoIterator<Item = Option<GSEARunning>>>(iter: I) -> Self {
        Self::new(iter.into_iter().collect::<Vec<Option<GSEARunning>>>())
    }
}

impl FromParallelIterator<Option<GSEARunning>> for GSEAPermutate {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = Option<GSEARunning>>,
    {
        Self::new(
            par_iter
                .into_par_iter()
                .collect::<Vec<Option<GSEARunning>>>(),
        )
    }
}

impl GSEAPermutate {
    fn new(value: Vec<Option<GSEARunning>>) -> Self {
        Self(value)
    }

    pub fn normalize_and_test(self) -> GSEAOutput {
        let (nes, (pvalue, null_nes_list)): (
            Vec<f64>,
            (Vec<f64>, Vec<Vec<f64>>),
        ) = self
            .0
            .par_iter()
            .map(|x| {
                x.as_ref().map_or_else(
                    || (f64::NAN, (f64::NAN, Vec::new())),
                    |running| running.normalize_and_test(),
                )
            })
            .unzip();
        let fdr = self.fdr(&nes, &null_nes_list);

        let n = self.0.len();
        let mut running_es = Vec::with_capacity(n);
        let mut es: Vec<f64> = Vec::with_capacity(n);
        let mut es_pos = Vec::with_capacity(n);
        for o in self.0 {
            match o {
                Some(running) => {
                    running_es.push(running.score.running);
                    es.push(running.score.es);
                    es_pos.push(running.score.es_pos);
                }
                None => {
                    running_es.push(Vec::new());
                    es.push(f64::NAN);
                    es_pos.push(0);
                }
            }
        }
        GSEAOutput {
            running_es,
            es,
            es_pos,
            nes,
            pvalue,
            fdr,
        }
    }

    fn fdr(&self, nes_list: &[f64], null_nes_list: &[Vec<f64>]) -> Vec<f64> {
        // we remove empty null (which means no gene are found in the ranking metrics),
        // and normalize the null enrichment scores
        let null_list: Vec<&Vec<f64>> = null_nes_list
            .par_iter()
            .filter(|null| null.len().gt(&0usize))
            .collect();

        // early exit when all genesets are failed
        if null_list.len().le(&0usize) {
            return vec![f64::NAN; self.0.len()];
        }

        let nperm = unsafe { null_list.get_unchecked(0).len() };
        // early exit when no permutations, but this shouldn't occur
        // Since we have check argument in the R-side
        if nperm.le(&0usize) {
            return vec![f64::NAN; self.0.len()];
        }

        // Note we won't use fraction, since the speed of rust is fast
        // From line 580-596?
        // The original method uses the last sample of sample permutations
        // as the observed enrichment score????
        // https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
        // From line 823-876
        // here we just regard `obs.phi.norm` is from the actual data
        // instead of using the last sample of sample permutations
        // prepare the observe enrichment score for both direction
        let obs_pos: Vec<&f64> =
            nes_list.iter().filter(|nes| **nes >= 0.0).collect(); // obs.count.col.norm
        let n_obs_pos = obs_pos.len() as f64;
        let obs_neg: Vec<&f64> =
            nes_list.iter().filter(|nes| **nes <= 0.0).collect(); // obs.count.col.norm
        let n_obs_neg = obs_neg.len() as f64;

        // in each permutation, we take `es` from all genesets as the null distribution
        // we split the null distribution into two direction
        let (reference_pos_list, reference_neg_list): (
            Vec<Vec<&f64>>,
            Vec<Vec<&f64>>,
        ) = (0 .. nperm)
            .into_par_iter()
            .map(|i| -> (Vec<&f64>, Vec<&f64>) {
                let null: Vec<&f64> = null_list
                    .iter()
                    .map(|null| unsafe { null.get_unchecked(i) }) // for performance
                    .collect::<Vec<&f64>>();
                let reference_pos: Vec<&f64> = null
                    .iter()
                    .filter(|score| ***score >= 0.0)
                    .copied()
                    .collect();
                let reference_neg: Vec<&f64> = null
                    .iter()
                    .filter(|score| ***score <= 0.0)
                    .copied()
                    .collect();
                (reference_pos, reference_neg)
            })
            .unzip();

        // see line 823-876: https://github.com/GSEA-MSigDB/GSEA_R/blob/master/R/GSEA.R
        self.0
            .par_iter()
            .zip(nes_list)
            .map(|(x, nes)| -> f64 {
                x.as_ref().map_or(f64::NAN, |_| {
                    let reference_props: Vec<f64>; // count.col
                    let obs_prop: f64; // obs.count.col
                    if *nes >= 0.0 {
                        let obs_count = obs_pos
                            .iter()
                            .filter(|score| **score >= nes)
                            .count();
                        if n_obs_pos.gt(&0.0) {
                            obs_prop = (obs_count as f64) / n_obs_pos;
                        } else {
                            obs_prop = 0.0;
                        }
                        reference_props = reference_pos_list
                            .iter()
                            .map(|reference_pos| {
                                let reference_count = reference_pos
                                    .iter()
                                    .filter(|score| **score >= nes)
                                    .count();
                                let n_reference = reference_pos.len() as f64;
                                if n_reference.gt(&0.0) {
                                    (reference_count as f64).div(n_reference)
                                } else {
                                    0.0f64
                                }
                            })
                            .collect();
                    } else {
                        let obs_count = obs_neg
                            .iter()
                            .filter(|score| **score <= nes)
                            .count();
                        if n_obs_neg.gt(&0.0) {
                            obs_prop = (obs_count as f64) / n_obs_neg;
                        } else {
                            obs_prop = 0.0;
                        }
                        reference_props = reference_neg_list
                            .iter()
                            .map(|reference_neg| {
                                let reference_count = reference_neg
                                    .iter()
                                    .filter(|score| **score <= nes)
                                    .count();
                                let n_reference = reference_neg.len() as f64;
                                if n_reference.gt(&0.0) {
                                    (reference_count as f64).div(n_reference)
                                } else {
                                    0.0f64
                                }
                            })
                            .collect();
                    }
                    reference_props
                        .iter()
                        .sum::<f64>()
                        .div(reference_props.len() as f64)
                        .div(obs_prop)
                        .min(1.0f64)
                })
            })
            .collect()
    }
}

#[cfg(test)]
mod test_gsea_input {
    use super::*;

    fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn test_norm_pos_and_neg() {
        let weights = vec![1.0, 2.0, 3.0, 4.0];
        let hits = vec![1, 3];
        let gsea = GSEAInput::new(&weights, &hits);

        let norm_pos = gsea.norm_pos(None);
        let norm_neg = gsea.norm_neg(None);

        assert!(approx_eq(norm_pos, 1.0 / (2.0 + 4.0), 1e-10));
        assert!(approx_eq(norm_neg, 1.0 / 2.0, 1e-10));
    }

    #[test]
    fn test_score() {
        let weights = vec![0.1, 0.2, 0.3, 0.4];
        let hits = vec![1, 3];
        let gsea = GSEAInput::new(&weights, &hits);

        let result = gsea.score(None, None, None);
        assert_eq!(result.running.len(), weights.len());
        assert!(result.es_pos < weights.len());
        assert!(result.es.abs() > 0.0);
    }

    #[test]
    fn test_es_consistency_with_score() {
        let weights = vec![0.1, 0.2, 0.3, 0.4];
        let hits = vec![1, 3];
        let gsea = GSEAInput::new(&weights, &hits);

        let es_score = gsea.score(None, None, None).es;
        let es_fast = gsea.es(None, None, None);
        assert!(approx_eq(es_score, es_fast, 1e-10));
    }

    #[test]
    fn test_switch_hits_override() {
        let weights = vec![1.0, 1.0, 1.0];
        let default_hits = vec![0];
        let override_hits = vec![1];
        let gsea = GSEAInput::new(&weights, &default_hits);

        let h1 = gsea.switch_hits(None);
        let h2 = gsea.switch_hits(Some(&override_hits));
        assert_eq!(h1, default_hits);
        assert_eq!(h2, override_hits);
    }
}

// Need nightly channel, always add `+nightly` when running `cargo +nightly bench`
#[cfg(all(feature = "bench", test))]
mod bench_gsea_input {
    extern crate test;
    use test::Bencher;

    use super::*;

    #[bench]
    fn bench_score(b: &mut Bencher) {
        let weights: Vec<f64> =
            (0 .. 20000).map(|i| (i as f64 + 1.0) / 20000.0).collect();
        let hits: Vec<usize> = (0 .. 20000).step_by(100).collect();
        let gsea = GSEAInput::new(&weights, &hits);
        b.iter(|| gsea.score(None, None, None));
    }

    #[bench]
    fn bench_es(b: &mut Bencher) {
        let weights: Vec<f64> =
            (0 .. 20000).map(|i| (i as f64 + 1.0) / 20000.0).collect();
        let hits: Vec<usize> = (0 .. 20000).step_by(100).collect();
        let gsea = GSEAInput::new(&weights, &hits);
        b.iter(|| gsea.es(None, None, None));
    }
}
