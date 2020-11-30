extern crate acc2tax;

use std::time::Instant;

#[test]
fn process_taxonomy() {
    let now = Instant::now();
    println!("Parsing (or loading) taxonomy related files...");
    acc2tax::init(
        64,
        "/mnt/data/nt/nucl_gb.accession2taxid.gz".to_string(),
        "/mnt/data/nt/taxdmp/nodes.dmp".to_string(),
        "/mnt/data/nt/taxdmp/names.dmp".to_string(),
    );
    println!("Time to parse and save file... {}", now.elapsed().as_secs());

    assert_eq!(
        acc2tax::get_taxon("X59856.2".to_string()),
        9913,
        "Get taxon error"
    );

    assert_eq!(
        acc2tax::get_complete_taxonomy(9913),
        [
            9913, 9903, 27592, 9895, 35500, 9845, 91561, 314145, 1437010, 9347, 32525, 40674,
            32524, 32523, 1338369, 8287, 117571, 117570, 7776, 7742, 89593, 7711, 33511, 33213,
            6072, 33208, 33154, 2759, 131567, 1
        ],
        "Get complete taxonomy error"
    );
}
