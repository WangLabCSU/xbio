use extendr_api::Robj;

#[allow(unused)]
fn check_metrics(metrics: &Robj) -> Result<&[f64], String> {
    metrics
        .as_real_slice()
        .ok_or("expected a numeric vector".to_string())
}
