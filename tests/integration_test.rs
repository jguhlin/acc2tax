extern crate acc2tax;

use std::time::{Instant};

#[test]
fn process_taxonomy() {
    let now = Instant::now();
    println!("Parsing (or loading) taxonomy related files...");
    acc2tax::init(64, 
                "/mnt/data/nt/nucl_gb.accession2taxid.gz".to_string(), 
                "/mnt/data/nt/taxdmp/nodes.dmp".to_string(),
                "/mnt/data/nt/taxdmp/names.dmp".to_string());
    println!("Time to parse and save file... {}", now.elapsed().as_secs());

    assert_eq!(acc2tax::check_taxon("X59856.2".to_string()), 9913);

    
    
}
