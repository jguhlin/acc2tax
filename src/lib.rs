#![feature(shrink_to)]

extern crate bincode;
extern crate flate2;
extern crate itertools;
extern crate bytelines;
extern crate byteorder;
extern crate rand;
extern crate serde;
extern crate snap;
extern crate zerocopy;
extern crate mimalloc;
extern crate rocksdb;
extern crate rayon;

use byteorder::BigEndian;
use byteorder::ByteOrder;
use zerocopy::{U32};
use rocksdb::{DB, Options};
use rayon::prelude::*;

use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

mod parser;
// mod fasta;

use hashbrown::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::str;
use std::sync::Arc;
use std::thread::Builder;

use once_cell::sync::OnceCell;

use itertools::Itertools;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crossbeam::atomic::AtomicCell;
use crossbeam::queue::ArrayQueue;
use crossbeam::utils::Backoff;

use bytelines::*;

use thincollections::thin_vec::ThinVec;

pub type TaxonLevels2AccInner = HashMap<u32, Vec<ThinVec<u8>>>;
pub type TaxonLevels2Acc = HashMap<u32, TaxonLevels2AccInner>;
pub type Result = (String, U32<BigEndian>);

static NAMES: OnceCell<Vec<String>> = OnceCell::new();
static TAXON2PARENT: OnceCell<Vec<usize>> = OnceCell::new();
static TAXON_RANK: OnceCell<Vec<String>> = OnceCell::new();
static ACC2TAX: OnceCell<Arc<DB>> = OnceCell::new();
static TAXIDS: OnceCell<Vec<u32>> = OnceCell::new();

// TODO: Update this and fasta annotator stuff...
/*fn load_taxon(filename: &str) -> Option<Acc2TaxInner> {
    let file = match File::open(filename) {
        Ok(file) => file,
        Err(_x)   => return None
    };

    let fh = BufReader::with_capacity(256 * 1024, file);
    Some(bincode::deserialize_from(snap::read::FrameDecoder::new(fh)).expect("Unable to read file..."))
}*/

#[pyfunction]
pub fn get_taxon(accession: String) -> u32 {
    let accession = accession
        .split_ascii_whitespace()
        .take(1)
        .collect::<String>();
    let acc2tax = ACC2TAX.get().expect("Data not initialized");
    let x = match acc2tax.get(accession) {
        Ok(x) => x,
        Err(y) => panic!("DB Error: {}", y),
    };

    match x {
        None => 0,
        Some(x) => BigEndian::read_u32(&x),
    }
}

fn load_existing() -> (Vec<String>, Vec<u32>, Vec<usize>, Vec<String>) {
    let names_fh = BufReader::with_capacity(2 * 1024 * 1024, File::open("names.bc").unwrap());
    let names = bincode::deserialize_from(snap::read::FrameDecoder::new(names_fh))
        .expect("Unable to read file...");

    let taxids_fh = BufReader::with_capacity(2 * 1024 * 1024, File::open("taxids.bc").unwrap());
    let taxids = bincode::deserialize_from(snap::read::FrameDecoder::new(taxids_fh))
            .expect("Unable to read file...");

    let t2p_fh = BufReader::with_capacity(2 * 1024 * 1024, File::open("t2p.bc").unwrap());
    let taxon_to_parent = bincode::deserialize_from(snap::read::FrameDecoder::new(t2p_fh))
        .expect("Unable to read file...");

    let taxon_rank_fh =
        BufReader::with_capacity(2 * 1024 * 1024, File::open("taxon_rank.bc").unwrap());
    let taxon_rank = bincode::deserialize_from(snap::read::FrameDecoder::new(taxon_rank_fh))
        .expect("Unable to read file...");

    (names, taxids, taxon_to_parent, taxon_rank)
}

#[pyfunction]
fn get_complete_taxonomy(taxon: usize) -> Vec<usize> {
    let mut complete_taxon: Vec<usize> = Vec::with_capacity(20);

    let mut cur_taxon = taxon;
    let taxon_to_parent = TAXON2PARENT.get().expect("Data not initialized");

    complete_taxon.push(cur_taxon);

    while cur_taxon != 1 && cur_taxon != 0 {
        cur_taxon = *taxon_to_parent.get(cur_taxon).or(Some(&1)).unwrap();
        complete_taxon.push(cur_taxon);
    }

    complete_taxon.shrink_to_fit();
    complete_taxon
}

#[pyfunction]
pub fn get_complete_taxonomy_dict(taxon: usize) -> HashMap<String, usize> {
    let mut complete_taxon: HashMap<String, usize> = HashMap::with_capacity(20);

    let mut cur_taxon = taxon;
    let taxon_to_parent = TAXON2PARENT.get().expect("Data not initialized");
    let taxon_rank = TAXON_RANK.get().expect("Taxon Rank not initialized");

    complete_taxon.insert(taxon_rank[cur_taxon].clone(), cur_taxon);

    while cur_taxon != 1 && cur_taxon != 0 {
        cur_taxon = *taxon_to_parent.get(cur_taxon).or(Some(&1)).unwrap();
        complete_taxon.insert(taxon_rank[cur_taxon].clone(), cur_taxon);
    }

    complete_taxon.shrink_to_fit();
    complete_taxon
}

#[pyfunction]
pub fn get_complete_taxonomy_names_dict(taxon: usize) -> HashMap<String, String> {
    let mut complete_taxon: HashMap<String, String> = HashMap::with_capacity(20);

    let mut cur_taxon = taxon;
    let taxon_to_parent = TAXON2PARENT.get().expect("Data not initialized");
    let taxon_rank = TAXON_RANK.get().expect("Taxon Rank not initialized");
    let names = NAMES.get().expect("Names not initialized");

    complete_taxon.insert(taxon_rank[cur_taxon].clone(), names[cur_taxon].clone());

    let mut loops: usize = 0;

    while cur_taxon != 1 && cur_taxon != 0 {
        loops += 1;
        if loops > 100 {
            println!(
                "Stuck in some sort of a loop... {} {}",
                cur_taxon,
                *taxon_to_parent.get(cur_taxon).unwrap()
            );
        }
        cur_taxon = *taxon_to_parent.get(cur_taxon).or(Some(&1)).unwrap();
        complete_taxon.insert(taxon_rank[cur_taxon].clone(), names[cur_taxon].clone());
    }
    complete_taxon.shrink_to_fit();
    complete_taxon
}

#[pyfunction]
pub fn init(
    num_threads: usize,
    acc2tax_filename: String,
    nodes_filename: String,
    names_filename: String,
) {
    if TAXON_RANK.get().is_some() { return }

    let (data, names, taxon_to_parent, taxon_rank, taxids);
    // let mut new = false;

    let mut rockoptions = Options::default();
    rockoptions.increase_parallelism(num_threads as i32);
    rockoptions.set_level_compaction_dynamic_level_bytes(true);
    rockoptions.set_compaction_readahead_size(8 * 1024 * 1024);
    
    // let a2tdb = Arc::new(sledconfig.open().expect("Unable to open database path, delete existing acc2tax.db and re-initialize. Check for free space on device."));
    // let a2tdb: Arc<sled::Db> = Arc::new(sled::open("acc2tax.db").expect("Unable to open database path, delete existing acc2tax.db and re-initialize. Check for free space on device."));
    if Path::new("taxon_rank.bc").exists() {
    // if a2tdb.len() > 0 {
        data = load_existing();
    } else {
        let mut writeoptions = rockoptions.clone();
        writeoptions.create_if_missing(true);
        writeoptions.set_bytes_per_sync(8 * 1024 * 1024);
        writeoptions.set_db_write_buffer_size(128 * 1024 * 1024);
        writeoptions.set_write_buffer_size(256 * 1024 * 1024);
        writeoptions.prepare_for_bulk_load();

        let a2tdb = Arc::new(DB::open(&writeoptions, "acc2tax.db").expect("Unable to open db for writing"));

        println!("Binary files do not exist, generating... This can take up to 60 minutes the first time...");
        println!("Processing acc2tax file");
        data = parser::read_taxonomy(
            num_threads,
            Arc::clone(&a2tdb),
            acc2tax_filename,
            nodes_filename,
            names_filename,
        );

        println!("Processing names file");

        // let mut acc2tax_fh = snap::Writer::new(File::create("acc2tax.bc").unwrap());
        // bincode::serialize_into(&mut acc2tax_fh, &data.0).expect("Unable to write to bincode file");
        // serde_json::to_writer(acc2tax_fh, &acc2tax).expect("Unable to write JSON file...");

        let mut names_fh = snap::write::FrameEncoder::new(File::create("names.bc").unwrap());
        bincode::serialize_into(&mut names_fh, &data.0).expect("Unable to write to bincode file");

        let mut taxids_fh = snap::write::FrameEncoder::new(File::create("taxids.bc").unwrap());
        bincode::serialize_into(&mut taxids_fh, &data.1).expect("Unable to write to bincode file");

        println!("Processing taxon to parent file");
        let mut taxon_to_parent_fh =
            snap::write::FrameEncoder::new(File::create("t2p.bc").unwrap());
        bincode::serialize_into(&mut taxon_to_parent_fh, &data.2)
            .expect("Unable to write to bincode file");

        println!("Processing taxon rank file");
        let mut taxon_rank_fh =
            snap::write::FrameEncoder::new(File::create("taxon_rank.bc").unwrap());
        bincode::serialize_into(&mut taxon_rank_fh, &data.3)
            .expect("Unable to write to bincode file");

        println!("Finished... Flushing and compacting database...");
        a2tdb.flush().expect("Unable to flush");
        a2tdb.compact_range(None::<&[u8]>, None::<&[u8]>);
    }

    let a2tdb = Arc::new(DB::open(&rockoptions, "acc2tax.db").expect("Unable to open database path, delete existing acc2tax.db and re-initialize"));

    names = data.0;
    taxids = data.1;
    taxon_to_parent = data.2;
    taxon_rank = data.3;

    NAMES
        .set(names)
        .expect("Unable to set. Already initialized?");
    TAXIDS
        .set(taxids)
        .expect("Unable to set. Already initialized?");
    TAXON2PARENT
        .set(taxon_to_parent)
        .expect("Unable to set. Already initialized?");
    TAXON_RANK
        .set(taxon_rank)
        .expect("Unable to set. Already initialized?");
    ACC2TAX
        .set(a2tdb)
        .expect("Unable to set. Already initialized?");

    println!("Loaded taxonomy databases");
}

#[pyfunction]
fn get_taxons_count() -> usize {
    let names = NAMES.get().expect("Names not initialized");
    names.len()
}

#[pyfunction]
fn get_child_taxons(parent_taxon: usize) -> Vec<usize> {
    let mut child_taxons: Vec<usize>; // = Vec::with_capacity(1000);

    let taxon_to_parent = TAXON2PARENT.get().expect("Data not initialized!");
    
    // child_taxons = taxon_to_parent.into_par_iter().enumerate().filter(|(_i, x)| **x == parent_taxon).map(|(i, _x)| i).collect();
    // child_taxons.par_sort();

    child_taxons = taxon_to_parent.iter().enumerate().filter(|(_i, x)| **x == parent_taxon).map(|(i, _x)| i).collect();

/*    for (x, item) in taxon_to_parent.iter().enumerate() {
        if *item == parent_taxon {
            child_taxons.push(x)
        }
    } */

    // child_taxons.shrink_to_fit();
    child_taxons.sort_unstable();
    child_taxons
}

#[pyfunction]
fn get_parent_taxons(taxon: usize) -> Vec<usize> {
    let mut parent_taxons: Vec<usize> = Vec::with_capacity(100);

    let taxon_to_parent = TAXON2PARENT.get().expect("Data not initialized!");

    let mut cur_tax_id = taxon;
    while cur_tax_id > 0 {
        cur_tax_id = taxon_to_parent[cur_tax_id];
        parent_taxons.push(cur_tax_id);
        if cur_tax_id == 1 {
            break;
        }
    }

    parent_taxons.shrink_to_fit();
    parent_taxons
}

#[pyfunction]
fn get_parent_taxons_names(taxon: usize) -> Vec<(usize, String)> {
    let parent_taxons = get_parent_taxons(taxon);
    let mut parent_taxons_names = Vec::with_capacity(parent_taxons.len());

    let names = NAMES.get().expect("Names not initialized");

    for parent in parent_taxons {
        parent_taxons_names.push((parent, names.get(parent).unwrap().to_string()))
    }

    parent_taxons_names
}

#[pyfunction]
fn get_child_taxons_names(parent_taxon: usize) -> Vec<(usize, String)> {
    let child_taxons = get_child_taxons(parent_taxon);
    let mut child_taxons_names = Vec::with_capacity(child_taxons.len());

    let names = NAMES.get().expect("Names not initialized");

    for child in child_taxons {
        child_taxons_names.push((child, names.get(child).unwrap().to_string()))
    }

    child_taxons_names
}

#[pyfunction]
fn get_taxon_rank(taxon: usize) -> String {
    let taxon_rank = TAXON_RANK.get().expect("Taxon Rank not initialized");
    taxon_rank.get(taxon).unwrap().to_string()
}

#[pyfunction]
fn get_taxon_name(taxon: usize) -> String {
    let taxon_names = NAMES.get().expect("Taxon Names not initialized");
    taxon_names
        .get(taxon)
        .unwrap_or(&"Not Found".to_string())
        .to_string()
}

#[pyfunction]
fn get_all_taxids() -> Vec<u32> {
    return TAXIDS.get().expect("Taxids not initialized").to_vec()
}

#[pymodule]
fn acc2tax(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pyfunction!(get_all_taxids))?;
    m.add_wrapped(wrap_pyfunction!(init))?;
    m.add_wrapped(wrap_pyfunction!(get_taxons_count))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon))?;
    m.add_wrapped(wrap_pyfunction!(get_complete_taxonomy))?;
    m.add_wrapped(wrap_pyfunction!(get_complete_taxonomy_dict))?;
    m.add_wrapped(wrap_pyfunction!(get_complete_taxonomy_names_dict))?;
    //m.add_wrapped(wrap_pyfunction!(convert_ntfasta_file))?;
    //m.add_wrapped(wrap_pyfunction!(filter_annotated_file))?;
    //m.add_wrapped(wrap_pyfunction!(filter_annotated_file_singlethreaded))?;
    m.add_wrapped(wrap_pyfunction!(get_parent_taxons))?;
    m.add_wrapped(wrap_pyfunction!(get_parent_taxons_names))?;
    m.add_wrapped(wrap_pyfunction!(get_child_taxons))?;
    m.add_wrapped(wrap_pyfunction!(get_child_taxons_names))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon_rank))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon_name))?;
    //m.add_wrapped(wrap_pyfunction!(child_taxon_seqlengths))?;
    //m.add_wrapped(wrap_pyfunction!(split_train_test_validation))?;
    //m.add_wrapped(wrap_pyfunction!(chunk_file))?;
    //m.add_wrapped(wrap_pyfunction!(shuffle_file))?;
    //m.add_wrapped(wrap_pyfunction!(balance_sequences))?;

    // Need a function to get the parents
    // Need a function to get the parents & ranks given a taxon (or accession)
    Ok(())
}

/*#[pyfunction]
fn split_train_test_validation(filename: String, output_prefix: String, test_pct: f32, validation_pct: f32)
{
    fasta::split_train_test_validation(filename, output_prefix, test_pct, validation_pct);
}

#[pyfunction]
fn chunk_file(filename: String, output_filename: String, size: usize)
{
    fasta::chunk_file(filename, output_filename, size);
}

#[pyfunction]
fn shuffle_file(filename: String, output_filename: String, bins: usize)
{
    fasta::shuffle_file(filename, output_filename, bins);
}

#[pyfunction]
fn child_taxon_seqlengths(filename: String, parent_tax_id: usize) {
    fasta::child_taxon_seqlengths(filename, parent_tax_id);
}


#[pyfunction]
fn filter_annotated_file(filename: String,
                         tax_id: usize,
                         num_threads: usize,)
{
    fasta::filter_annotated_file(filename, tax_id, num_threads);
}

#[pyfunction]
fn filter_annotated_file_singlethreaded(
                        filename: String,
                        tax_id: usize,)
{
    fasta::filter_annotated_file_singlethreaded(filename, tax_id);
}

#[pyfunction]
fn convert_ntfasta_file(filename: String,
                    output: String,
                    num_threads: usize,)
{
    fasta::convert_ntfasta_file(filename, output, num_threads);
}

#[pyfunction]
fn balance_sequences(filename: String,
    output_filename: String,
    parent_tax_id: usize, // Have to balance relative to a parent taxon
    minimum_length: usize)
{
    fasta::balance_sequences(filename, output_filename, parent_tax_id, minimum_length);
} */
