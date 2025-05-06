use core::f64;
use std::{collections::HashSet, ops::Not};

use extendr_api::{list, List, Robj};
use rayon::prelude::*;

pub struct GSEAOutput {
    pub running_es: Vec<Vec<f64>>,
    pub es: Vec<f64>,
    pub es_pos: Vec<usize>,
    pub nes: Vec<f64>,
    pub pvalue: Vec<f64>,
}

impl GSEAOutput {
    pub fn new(n: usize) -> Self {
        Self {
            running_es: Vec::with_capacity(n),
            es: Vec::with_capacity(n),
            es_pos: Vec::with_capacity(n),
            nes: Vec::with_capacity(n),
            pvalue: Vec::with_capacity(n),
        }
    }

    pub fn add_output(
        &mut self,
        running: Vec<f64>,
        es_pos: usize,
        es: f64,
        null: &[f64],
    ) {
        let mean = null.iter().sum::<f64>() / (null.len() as f64);
        self.pvalue.push(gsea_pvalue(es, null));
        self.nes.push(gsea_nes(es, mean));
        self.es.push(es);
        self.es_pos.push(es_pos);
        self.running_es.push(running);
    }

    pub fn add_empty(&mut self) {
        self.pvalue.push(1.0);
        self.nes.push(f64::NAN);
        self.es.push(f64::NAN);
        self.es_pos.push(0);
        self.running_es.push(vec![]);
    }
}

impl From<GSEAOutput> for Robj {
    fn from(value: GSEAOutput) -> Self {
        let running_es_list: List = value.running_es.into_iter().collect();
        list!(
            running_es = running_es_list,
            es = value.es,
            es_pos = value.es_pos,
            nes = value.nes,
            pvalue = value.pvalue
        )
        .into()
    }
}

pub fn gsea_par_prerank<'a>(
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

pub fn gsea_prerank<'a>(
    identifiers: &[&'a str],
    metrics: &[f64],
    exponent: f64,
) -> (Vec<usize>, Vec<&'a str>, Vec<f64>) {
    // Create a list of indices
    let mut indices: Vec<usize> = (0 .. metrics.len()).collect();

    // Sort indices by descending metrics
    indices.sort_unstable_by(|a, b| f64::total_cmp(&metrics[*b], &metrics[*a]));

    // Define the weights
    let weights: Vec<f64> =
        metrics.iter().map(|x| x.abs().powf(exponent)).collect();

    // re-order the identifiers and weights
    let identifiers = indices
        .iter()
        .map(|&i| identifiers[i])
        .collect::<Vec<&str>>();
    let weights = indices.iter().map(|&i| weights[i]).collect();
    (indices, identifiers, weights)
}

pub fn gsea_par_hits(
    geneset: &HashSet<&str>,
    identifiers: &[&str],
) -> Vec<bool> {
    identifiers
        .par_iter()
        .map(|id| geneset.contains(id))
        .collect::<Vec<bool>>()
}

pub fn gsea_hits(geneset: &HashSet<&str>, identifiers: &[&str]) -> Vec<bool> {
    identifiers
        .iter()
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
        .fold(0.0f64, |old, (hit, metric)| -> f64 {
            let score = if *hit { metric * norm_pos } else { -norm_neg };
            let new = old + score;
            if old.abs() > new.abs() {
                old
            } else {
                new
            }
        })
}

pub fn gsea_nes(es: f64, mean: f64) -> f64 {
    es / mean.abs()
}

pub fn gsea_pvalue(es: f64, null: &[f64]) -> f64 {
    let total: f64 = null.len() as f64;
    let n: f64;
    if es >= 0.0 {
        n = null.iter().filter(|score| **score >= es).count() as f64;
    } else {
        n = null.iter().filter(|score| **score <= es).count() as f64;
    }
    n / total
}

pub fn gsea_fdr(null: &[f64], mean: f64) {}
