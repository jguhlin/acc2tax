use byteorder::{BigEndian, LittleEndian};
use byteorder::ByteOrder;
use rocksdb::{Options, DB};

use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

mod parser;

use std::collections::{HashMap, HashSet};
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
pub type Result = (String, u32);

static NAMES: OnceCell<Vec<String>> = OnceCell::new();
static TAXON2PARENT: OnceCell<Vec<usize>> = OnceCell::new();
static TAXON_RANK: OnceCell<Vec<String>> = OnceCell::new();
static ACC2TAX: OnceCell<Arc<DB>> = OnceCell::new();
static TAXIDS: OnceCell<Vec<u32>> = OnceCell::new();

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
        Some(x) => LittleEndian::read_u32(&x),
    }
}

fn load_existing() -> (Vec<String>, Vec<u32>, Vec<usize>, Vec<String>) {
    let names_fh = BufReader::with_capacity(2 * 1024 * 1024, File::open("names.bc").unwrap());
    let bincoded: Vec<u8> = bincode::deserialize_from(names_fh).expect("unable to read names.bc");
    let decompressed = zstd::stream::decode_all(&bincoded[..]).expect("unable to read names.bc");
    let names: Vec<String>  = bincode::deserialize_from(&decompressed[..]).expect("Unable to read file...");

    let taxids_fh = BufReader::with_capacity(2 * 1024 * 1024, File::open("taxids.bc").unwrap());
    let bincoded: Vec<u8> = bincode::deserialize_from(taxids_fh).expect("unable to read taxids.bc");
    let decompressed = zstd::stream::decode_all(&bincoded[..]).expect("unable to read taxids.bc");
    let taxids: Vec<u32>  = bincode::deserialize_from(&decompressed[..]).expect("Unable to read file...");

    let t2p_fh = BufReader::with_capacity(2 * 1024 * 1024, File::open("t2p.bc").unwrap());
    let bincoded: Vec<u8> = bincode::deserialize_from(t2p_fh).expect("unable to read t2p.bc");
    let decompressed = zstd::stream::decode_all(&bincoded[..]).expect("unable to read t2p.bc");
    let taxon_to_parent: Vec<usize>  = bincode::deserialize_from(&decompressed[..]).expect("Unable to read file...");

    let taxon_rank_fh =
        BufReader::with_capacity(2 * 1024 * 1024, File::open("taxon_rank.bc").unwrap());
    let bincoded: Vec<u8> = bincode::deserialize_from(taxon_rank_fh).expect("unable to read taxon_rank.bc");
    let decompressed = zstd::stream::decode_all(&bincoded[..]).expect("unable to read taxon_rank.bc");
    let taxon_rank: Vec<String>  = bincode::deserialize_from(&decompressed[..]).expect("Unable to read file...");

    (names, taxids, taxon_to_parent, taxon_rank)
}

#[pyfunction]
pub fn get_complete_taxonomy(taxon: usize) -> Vec<usize> {
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
    if TAXON_RANK.get().is_some() {
        return;
    }

    let (data, names, taxon_to_parent, taxon_rank, taxids);

    let mut rockoptions = Options::default();
    rockoptions.increase_parallelism(num_threads as i32);
    rockoptions.set_level_compaction_dynamic_level_bytes(true);

    // let a2tdb = Arc::new(sledconfig.open().expect("Unable to open database path, delete existing acc2tax.db and re-initialize. Check for free space on device."));
    // let a2tdb: Arc<sled::Db> = Arc::new(sled::open("acc2tax.db").expect("Unable to open database path, delete existing acc2tax.db and re-initialize. Check for free space on device."));
    if Path::new("taxon_rank.bc").exists() {
        // if a2tdb.len() > 0 {
        data = load_existing();
    } else {
        let mut writeoptions = rockoptions.clone();
        writeoptions.create_if_missing(true);
        writeoptions.prepare_for_bulk_load();

        let a2tdb =
            Arc::new(DB::open(&writeoptions, "acc2tax.db").expect("Unable to open db for writing"));

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

        let mut names_fh = File::create("names.bc").unwrap();
        let bincoded = bincode::serialize(&data.0).expect("Unable to write to bincode file");
        let compressed = zstd::encode_all(&bincoded[..], -3).expect("Unable to compress");
        bincode::serialize_into(&mut names_fh, &compressed).expect("Unable to write to bincode file");
        drop(names_fh);

        let mut taxids_fh = File::create("taxids.bc").unwrap();
        let bincoded = bincode::serialize(&data.1).expect("Unable to write to bincode file");
        let compressed = zstd::encode_all(&bincoded[..], -3).expect("Unable to compress");
        bincode::serialize_into(&mut taxids_fh, &compressed).expect("Unable to write to bincode file");
        drop(taxids_fh);
        
        println!("Processing taxon to parent file");
        let mut taxon_to_parent_fh = File::create("t2p.bc").unwrap();
        let bincoded = bincode::serialize(&data.2).expect("Unable to write to bincode file");
        let compressed = zstd::encode_all(&bincoded[..], -3).expect("Unable to compress");
        bincode::serialize_into(&mut taxon_to_parent_fh, &compressed).expect("Unable to write to bincode file");
        drop(taxon_to_parent_fh);

        println!("Processing taxon rank file");
        let mut taxon_rank_fh = File::create("taxon_rank.bc").unwrap();
        let bincoded = bincode::serialize(&data.3).expect("Unable to write to bincode file");
        let compressed = zstd::encode_all(&bincoded[..], -3).expect("Unable to compress");
        bincode::serialize_into(&mut taxon_rank_fh, &compressed).expect("Unable to write to bincode file");
        drop(taxon_rank_fh);

        println!("Finished... Flushing and compacting database...");
        a2tdb.flush().expect("Unable to flush");
        a2tdb.compact_range(None::<&[u8]>, None::<&[u8]>);
    }

    let a2tdb = Arc::new(
        DB::open(&rockoptions, "acc2tax.db")
            .expect("Unable to open database path, delete existing acc2tax.db and re-initialize"),
    );

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

    child_taxons = taxon_to_parent
        .iter()
        .enumerate()
        .filter(|(_i, x)| **x == parent_taxon)
        .map(|(i, _x)| i)
        .collect();

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
pub fn get_taxon_rank(taxon: usize) -> String {
    let taxon_rank = TAXON_RANK.get().expect("Taxon Rank not initialized");
    taxon_rank.get(taxon).unwrap().to_string()
}

#[pyfunction]
pub fn get_taxon_ranks() -> Vec<String> {
    let taxon_rank = TAXON_RANK.get().expect("Taxon Rank not initialized");
    taxon_rank.into_iter().collect::<HashSet<&String>>().into_iter().map(|x| x.to_string()).collect()
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
    return TAXIDS.get().expect("Taxids not initialized").to_vec();
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
    m.add_wrapped(wrap_pyfunction!(get_parent_taxons))?;
    m.add_wrapped(wrap_pyfunction!(get_parent_taxons_names))?;
    m.add_wrapped(wrap_pyfunction!(get_child_taxons))?;
    m.add_wrapped(wrap_pyfunction!(get_child_taxons_names))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon_rank))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon_ranks))?;
    m.add_wrapped(wrap_pyfunction!(get_taxon_name))?;
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
