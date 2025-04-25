use extendr_api::{prelude::Robj, RobjItertools};

// we should provide a method to convert this object into Robj
struct GSEAResult {
    es: Vec<f64>,
    running_es: Vec<f64>,
    nes: Vec<f64>,
    pvalue: Vec<f64>,
    qvalue: Vec<f64>,
}

fn gsea() {}

struct GSEAInput {
    metrics: Vec<f64>,
    genesets: Vec<String>,
    exponent: f64,
    nperm: usize,
    threads: usize,
    seed: usize,
}

fn running_es(metrics: &[f64], hits: &[bool]) -> Vec<f64> {
    let scores_pos: f64 = hits
        .iter()
        .zip(metrics.iter())
        .filter_map(|(&hit, metric)| if hit { Some(metric) } else { None })
        .sum();
    let norm_pos = 1.0 / scores_pos;
    let norm_neg = 1.0 / ((metrics.len() - hits.len()) as f64);
    hits.iter()
        .zip(metrics.iter())
        .map(|(&hit, metric)| if hit { metric * norm_pos } else { -norm_neg })
        .scan(0.0, |acc, score| {
            *acc += score;
            Some(*acc)
        })
        .collect()
}

fn es(running_es: &[f64]) -> f64 {}

fn gene_permutate(metrics: &[f64], hit_list: &[&[bool]], weight: f64) {
    let mut metrics: Vec<f64> = metrics.iter().map(|x| x.abs().powf(weight)).collect();
    metrics.sort_unstable_by(|a, b| f64::total_cmp(b, a));
    let mut es_list: Vec<f64> = Vec::with_capacity(hit_list.len());
    let mut running_es_list: Vec<Vec<f64>> = Vec::with_capacity(hit_list.len());
    for hits in hit_list {
        let running_escores = running_es(&metrics, hits);
        es_list.push(es(&running_escores));
        running_es_list.push(running_escores);
    }

    for hits in hit_list {}
}
