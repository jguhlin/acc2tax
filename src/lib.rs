extern crate flate2;
extern crate itertools;
extern crate bincode;
extern crate mimalloc;
extern crate serde;
extern crate bytelines;
extern crate snap;

use mimalloc::MiMalloc;

#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

mod parser;

use std::str;
use std::sync::Arc;
use std::io::{BufReader, BufRead};
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::collections::HashMap;
use std::thread::Builder;
use std::path::Path;

use once_cell::sync::OnceCell;

use twox_hash::XxHash;

use itertools::Itertools;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crossbeam::utils::Backoff;
use crossbeam::atomic::AtomicCell;
use crossbeam::queue::{ArrayQueue, PushError, SegQueue};

use bytelines::*;

pub type Acc2TaxInner = HashMap<Vec<u8>, u32, BuildHasherDefault<XxHash>>;
pub type Acc2Tax = HashMap<u32, Acc2TaxInner, BuildHasherDefault<XxHash>>;
pub type Result = (u32, Vec<u8>, u32);

static NAMES: OnceCell<Vec<String>>       = OnceCell::new();
static TAXON2PARENT: OnceCell<Vec<usize>> = OnceCell::new();
static TAXON_RANK: OnceCell<Vec<String>>  = OnceCell::new();

fn load_taxon(filename: &str) -> Acc2TaxInner {
    let fh = BufReader::with_capacity(64 * 1024 * 1024, File::open(filename).unwrap());
    bincode::deserialize_from(snap::Reader::new(fh)).expect("Unable to read file...")
}

#[pyfunction]
pub fn get_taxon(accession: String) -> u32 {

    let accession = accession.as_bytes().to_vec();
    let short: u32 = parser::shorten(&accession[0..4]);
    
    let map = load_taxon(&format!("acc2tax_db/{}.bc", &short.to_string()));

    *map.get(&accession).unwrap()
    
}

fn load_existing() -> (Option<Acc2Tax>, Vec<String>, Vec<usize>, Vec<String>) {
    let names_fh = BufReader::with_capacity(64 * 1024 * 1024, File::open("names.bc").unwrap());
    let names = bincode::deserialize_from(snap::Reader::new(names_fh)).expect("Unable to read file...");

    let t2p_fh = BufReader::with_capacity(64 * 1024 * 1024, File::open("t2p.bc").unwrap());
    let taxon_to_parent = bincode::deserialize_from(snap::Reader::new(t2p_fh)).expect("Unable to read file...");

    let taxon_rank_fh = BufReader::with_capacity(64 * 1024 * 1024, File::open("taxon_rank.bc").unwrap());
    let taxon_rank = bincode::deserialize_from(snap::Reader::new(taxon_rank_fh)).expect("Unable to read file...");

    (None, names, taxon_to_parent, taxon_rank)

}

#[pyfunction]
pub fn init(num_threads: usize, acc2tax_filename: String, nodes_filename: String, names_filename: String) {
// Initializes the database

    let (data, names, taxon_to_parent, taxon_rank);

    if Path::new("names.bc").exists() {
        data = load_existing();
    } else {
        println!("Binary files do not exist, generating... This can take up to 60 minutes the first time...");
        data = parser::read_taxonomy(num_threads, acc2tax_filename, nodes_filename, names_filename);

        for (short, all) in data.0.unwrap().iter() {
            let mut acc2tax_fh = snap::Writer::new(File::create(format!("acc2tax_db/{}.bc", short.to_string())).unwrap());
            bincode::serialize_into(&mut acc2tax_fh, &all).expect("Unable to write to bincode file");
        }

        // let mut acc2tax_fh = snap::Writer::new(File::create("acc2tax.bc").unwrap());
        // bincode::serialize_into(&mut acc2tax_fh, &data.0).expect("Unable to write to bincode file");
        // serde_json::to_writer(acc2tax_fh, &acc2tax).expect("Unable to write JSON file...");

        let mut names_fh = snap::Writer::new(File::create("names.bc").unwrap());
        bincode::serialize_into(&mut names_fh, &data.1).expect("Unable to write to bincode file");

        let mut taxon_to_parent_fh = snap::Writer::new(File::create("t2p.bc").unwrap());
        bincode::serialize_into(&mut taxon_to_parent_fh, &data.2).expect("Unable to write to bincode file");

        let mut taxon_rank_fh = snap::Writer::new(File::create("taxon_rank.bc").unwrap());
        bincode::serialize_into(&mut taxon_rank_fh, &data.3).expect("Unable to write to bincode file");

    }

    names = data.1;
    taxon_to_parent = data.2;
    taxon_rank = data.3;

    NAMES.set(names).expect("Unable to set. Already initialized?");
    TAXON2PARENT.set(taxon_to_parent).expect("Unable to set. Already initialized?");
    TAXON_RANK.set(taxon_rank).expect("Unable to set. Already initialized?");
    println!("Loaded taxonomy databases");
}

#[pymodule]
fn acc2tax(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(init))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon))?;
    // Need a function to get the rank of a taxon
    // Need a function to get the parents
    // Need a function to get the parents & ranks given a taxon (or accession)
    Ok(())
}
