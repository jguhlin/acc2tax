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
use std::collections::HashSet;

use once_cell::sync::OnceCell;

use twox_hash::XxHash;

use itertools::Itertools;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crossbeam::utils::Backoff;
use crossbeam::atomic::AtomicCell;
use crossbeam::queue::{ArrayQueue, PushError, SegQueue};

use bytelines::*;
use rayon::prelude::*;

use thincollections::thin_vec::ThinVec;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

pub type Acc2TaxInner = HashMap<Vec<u8>, u32, BuildHasherDefault<XxHash>>;
pub type Acc2Tax = HashMap<u32, Acc2TaxInner, BuildHasherDefault<XxHash>>;
pub type TaxonLevels2AccInner = HashMap<u32, Vec<ThinVec<u8>>, BuildHasherDefault<XxHash>>;
pub type TaxonLevels2Acc = HashMap<u32, TaxonLevels2AccInner, BuildHasherDefault<XxHash>>;
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
pub fn get_complete_taxonomy (taxon: usize) -> Vec<usize> {
    let mut complete_taxon: Vec<usize> = Vec::with_capacity(20);
    
    let mut cur_taxon = taxon;
    let taxon_to_parent = TAXON2PARENT.get().expect("Data not initialized");

    complete_taxon.push(cur_taxon);

    while cur_taxon != 1 {
        cur_taxon = *taxon_to_parent.get(cur_taxon).or(Some(&1)).unwrap();
        complete_taxon.push(cur_taxon);
    }

    complete_taxon.shrink_to_fit();
    complete_taxon
}

#[pyfunction]
pub fn init(num_threads: usize, acc2tax_filename: String, nodes_filename: String, names_filename: String) {
// Initializes the database

    let (data, names, taxon_to_parent, taxon_rank);
    let mut new = false;
    let mut acc2tax: Acc2Tax = Default::default();


    if Path::new("taxon_level_to_acc.bc").exists() {
        data = load_existing();
    } else {
        new = true;
        println!("Binary files do not exist, generating... This can take up to 60 minutes the first time...");
        data = parser::read_taxonomy(num_threads, acc2tax_filename, nodes_filename, names_filename);

        acc2tax = data.0.unwrap();

        // TODO: This part should be serialized, actually...
        acc2tax.par_iter().for_each(|(short, all)| {
            let mut acc2tax_fh = snap::Writer::new(File::create(format!("acc2tax_db/{}.bc", short.to_string())).unwrap());
            bincode::serialize_into(&mut acc2tax_fh, &all).expect("Unable to write to bincode file");
        });

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

    if new {
        let mut taxon_level_to_acc: TaxonLevels2Acc = Default::default();
        taxon_level_to_acc.reserve(100_000);

        let taxon_rank = TAXON_RANK.get().expect("Taxon rank not properly loaded...");
        let taxon_to_parent = TAXON2PARENT.get().expect("Taxon2Parent not properly loaded");

        // let mut accessions: Vec<ThinVec<u8>> = Vec::with_capacity(300_000_000); // 260_424_480
                                                                 
        /* for (_, mut all) in acc2tax.drain() {
            for (acc, _) in all.drain() {
                let mut accthin = ThinVec::with_capacity(acc.len());
                accthin.extend_from_slice(&acc[..]);
                accessions.push(accthin);
            }
        }

        drop(acc2tax);

        for acc in accessions {
            let tax_id = get_taxon(String::from_utf8(acc.to_vec()).unwrap());
            let mut work_tax_id = tax_id;
            loop {
                if work_tax_id == 0 { break; }
                if work_tax_id == 1 { break; }

                let parent_tax_id: u32 = taxon_to_parent[work_tax_id as usize] as u32;
                
                let entry = taxon_level_to_acc
                                .entry(parent_tax_id)
                                .or_insert_with(Vec::new);

                entry.push( (acc.clone(), work_tax_id) );

                // Iteratively move up...
                work_tax_id = parent_tax_id;
            }
        } */

        let pb = ProgressBar::new(acc2tax.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta} {msg}")
        .progress_chars("█▇▆▅▄▃▂▁  "));

        let mut i = 0;
        
        for (_, mut all) in pb.wrap_iter(acc2tax.drain()) {
            i += 1;
            let mut alllen = all.len();
            pb.set_message(&format!("{}", alllen));
            for (acc, tax_id) in all.drain() {
                
                let mut work_tax_id = tax_id;

                let mut accthin = ThinVec::with_capacity(acc.len());
                accthin.extend_from_slice(&acc[..]);

                loop {
                    if work_tax_id == 0 { break; }
                    if work_tax_id == 1 { break; }

                    let parent_tax_id: u32 = taxon_to_parent[work_tax_id as usize] as u32;
                    
                    let entry = taxon_level_to_acc
                                    .entry(parent_tax_id)
                                    .or_insert_with(Default::default);

                    let entry_inner = entry.entry(work_tax_id).or_insert_with(Vec::new);
                    entry_inner.push(accthin.clone());

                    // Should be another hashmap that is work_tax_id => accthin
                    // So we aren't duplicating all of this all the time...

                    // entry.push( (accthin, work_tax_id) );

                    // Iteratively move up...
                    assert_ne!(work_tax_id, parent_tax_id);

                    if i >= 7898 {
                        println!("{}", i);
                        println!("{} {}", work_tax_id, parent_tax_id);
                    }

                    work_tax_id = parent_tax_id;
                
                }
            }
        }

        // let mut taxon_level_to_acc_fh = snap::Writer::new(File::create("taxon_level_to_acc.bc").unwrap());
        // bincode::serialize_into(&mut taxon_level_to_acc_fh, &taxon_level_to_acc).expect("Unable to write to bincode file");
    }

    println!("Loaded taxonomy databases");
}

#[pymodule]
fn acc2tax(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(init))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon))?;
    m.add_wrapped(wrap_pyfunction!(get_complete_taxonomy))?;
    // Need a function to get the rank of a taxon
    // Need a function to get the parents
    // Need a function to get the parents & ranks given a taxon (or accession)
    Ok(())
}
