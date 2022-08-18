use diploid_contam_estimator::wrapper;

fn main() {
    // parse cli argumnets
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    let result = wrapper();
    if let Err(e) = result {
        println!("Error: {}", e)
    }
}
