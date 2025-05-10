use core::f64;
use std::{
    collections::HashSet,
    ops::{Div, Not},
};

use extendr_api::{list, List, Robj};
use rayon::prelude::*;

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

pub struct GSEAPermutateList(Vec<GSEAPermutate>);

impl FromIterator<GSEAPermutate> for GSEAPermutateList {
    fn from_iter<I: IntoIterator<Item = GSEAPermutate>>(iter: I) -> Self {
        Self::new(iter.into_iter().collect::<Vec<GSEAPermutate>>())
    }
}

impl FromParallelIterator<GSEAPermutate> for GSEAPermutateList {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: IntoParallelIterator<Item = GSEAPermutate>,
    {
        Self::new(par_iter.into_par_iter().collect::<Vec<GSEAPermutate>>())
    }
}

impl GSEAPermutateList {
    fn new(value: Vec<GSEAPermutate>) -> Self {
        Self(value)
    }

    pub fn normalize_and_test(self) -> GSEAOutput {
        let (nes, (pvalue, null_nes_list)): (
            Vec<f64>,
            (Vec<f64>, Vec<Vec<f64>>),
        ) = self.0.par_iter().map(|x| x.normalize_and_test()).unzip();
        let fdr = self.fdr(&nes, &null_nes_list);

        // let null_list: List = value.null.into_iter().collect();
        let n = self.0.len();
        let mut running_es = Vec::with_capacity(n);
        let mut es: Vec<f64> = Vec::with_capacity(n);
        let mut es_pos = Vec::with_capacity(n);
        for v in self.0 {
            running_es.push(v.running);
            es.push(v.es);
            es_pos.push(v.es_pos);
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

    pub fn fdr(
        &self,
        nes_list: &[f64],
        null_nes_list: &[Vec<f64>],
    ) -> Vec<f64> {
        // we remove empty null, and normalize the null enrichment scores
        let null_list: Vec<&Vec<f64>> = null_nes_list
            .par_iter()
            .filter(|null| null.len().gt(&0usize))
            .collect();

        // early exit when all genesets are failed
        if null_list.len().le(&0usize) {
            return vec![f64::NAN; self.0.len()];
        }

        // early exit when no permutations, but this shouldn't occur
        let nperm = unsafe { null_list.get_unchecked(0).len() };
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
                if !x.success {
                    return f64::NAN;
                }
                let reference_props: Vec<f64>; // count.col
                let obs_prop: f64; // obs.count.col
                if *nes >= 0.0 {
                    let obs_count =
                        obs_pos.iter().filter(|score| **score >= nes).count();
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
                    let obs_count =
                        obs_neg.iter().filter(|score| **score <= nes).count();
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
            .collect()
    }
}

pub struct GSEAPermutate {
    pub success: bool,
    pub running: Vec<f64>,
    pub es_pos: usize,
    pub es: f64,
    pub null: Vec<f64>,
}

impl GSEAPermutate {
    pub fn empty() -> Self {
        Self {
            success: false,
            running: Vec::new(),
            es_pos: 0,
            es: f64::NAN,
            null: Vec::new(),
        }
    }

    pub fn new(
        running: Vec<f64>,
        es_pos: usize,
        es: f64,
        null: Vec<f64>,
    ) -> Self {
        GSEAPermutate {
            success: true,
            running,
            es_pos,
            es,
            null,
        }
    }

    fn split_null(&self) -> (Vec<&f64>, Vec<&f64>) {
        let pos_null: Vec<&f64> =
            self.null.iter().filter(|es| **es >= 0.0f64).collect();
        let neg_null: Vec<&f64> =
            self.null.iter().filter(|es| **es <= 0.0f64).collect();
        (pos_null, neg_null)
    }

    // Get nes and pvalue
    fn nes(&self, pos_mean: f64, neg_mean: f64) -> f64 {
        if self.es.ge(&0.0) {
            self.es.div(pos_mean.abs())
        } else {
            self.es.div(neg_mean.abs())
        }
    }

    fn test(&self, pos_null: Vec<&f64>, neg_null: Vec<&f64>) -> f64 {
        let total: f64;
        let n: f64;
        if self.es.ge(&0.0) {
            total = pos_null.len() as f64;
            n = pos_null.iter().filter(|score| ***score >= self.es).count()
                as f64;
        } else {
            total = neg_null.len() as f64;
            n = neg_null.iter().filter(|score| ***score <= self.es).count()
                as f64;
        }
        n / total
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

    pub fn normalize_and_test(&self) -> (f64, (f64, Vec<f64>)) {
        if !self.success {
            return (f64::NAN, (f64::NAN, Vec::new()));
        }
        let (pos_null, neg_null) = self.split_null();
        let pos_mean = pos_null
            .iter()
            .copied()
            .sum::<f64>()
            .div(pos_null.len() as f64);
        let neg_mean = neg_null
            .iter()
            .copied()
            .sum::<f64>()
            .div(neg_null.len() as f64);
        let nes = self.nes(pos_mean, neg_mean);
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

pub fn gsea_hits(geneset: &HashSet<&str>, identifiers: &[&str]) -> Vec<bool> {
    identifiers
        .par_iter()
        .map(|id| geneset.contains(id))
        .collect::<Vec<bool>>()
}

pub fn gsea_running_es(
    weights: &[f64],
    hits: &[bool],
    norm_pos: f64,
    norm_neg: f64,
) -> (Vec<f64>, usize, f64) {
    let running: Vec<f64> = hits
        .iter()
        .zip(weights)
        .map(
            |(hit, weight)| {
                if *hit {
                    weight * norm_pos
                } else {
                    -norm_neg
                }
            },
        )
        .scan(0.0, |acc, score| {
            *acc += score;
            Some(*acc)
        })
        .collect();
    let (es_pos, &es) =
        running.iter().enumerate().fold((0usize, &0.0f64), |x, y| {
            if x.1.abs() > y.1.abs() {
                x
            } else {
                y
            }
        });
    (running, es_pos, es)
}

pub fn gsea_norm_pos(weights: &[f64], hits: &[bool]) -> f64 {
    let sum_pos: f64 = hits
        .iter()
        .zip(weights)
        .filter_map(|(hit, weight)| hit.then_some(weight))
        .sum();
    1.0 / sum_pos
}

pub fn gsea_norm_neg(hits: &[bool]) -> f64 {
    1.0 / (hits
        .iter()
        // Only keep elements not in the geneset
        .filter(|hit| hit.not())
        .count() as f64)
}

pub fn gsea_es(
    weights: &[f64],
    hits: &[bool],
    norm_pos: f64,
    norm_neg: f64,
) -> f64 {
    hits.iter()
        .zip(weights)
        .fold(
            (0.0f64, 0.0f64),
            |(score, es), (hit, metric)| -> (f64, f64) {
                let cumsum;
                if *hit {
                    cumsum = score + (metric * norm_pos)
                } else {
                    cumsum = score - norm_neg
                };
                if es.abs() > cumsum.abs() {
                    (cumsum, es)
                } else {
                    (cumsum, cumsum)
                }
            },
        )
        .1 // Get the `es`
}
