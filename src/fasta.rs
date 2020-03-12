use std::thread::{Builder, Thread, JoinHandle};
use std::thread;

use std::io::prelude::*;
use std::io::BufWriter;

use std::sync::{Arc, RwLock};

use std::fs::File;
use std::io::{BufReader, Read, BufRead};
use std::time::Duration;

use crate::parser;

use std::collections::HashMap;

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

use pyo3::prelude::*;

use rand::Rng;
use rand::prelude::*;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use super::{get_taxon, get_complete_taxonomy};

use serde::{Serialize, Deserialize};

#[derive(PartialEq)]
pub enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Sequence> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> Sequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

impl ThreadCommand<FastaSequence> {
    // Consumes the ThreadCommand, which is just fine...
    pub fn unwrap(self) -> FastaSequence {
        match self {
            ThreadCommand::Work(x)   => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

#[derive(PartialEq, Serialize, Deserialize)]
pub struct Sequence {
    pub seq: Vec<u8>,
    pub id:  String,
    pub taxons: Vec<usize>,
    pub taxon: usize,
}

#[derive(PartialEq, Serialize, Deserialize)]
pub struct FastaSequence {
    pub seq: Vec<u8>,
    pub id:  String,
}

pub struct ThreadPool {
    pub generator_done: RwLock<bool>,
    pub children: RwLock<Vec<JoinHandle<()>>>,
    pub num_threads: usize,
    pub main_thread_handle: Thread,
}

impl ThreadPool {
    pub fn new(num_threads: usize) -> Self {
        ThreadPool { 
            generator_done: RwLock::new(false),
            children: RwLock::new(Vec::new()),
            num_threads: num_threads,
            main_thread_handle: thread::current()
        }
    }
}

pub const STACKSIZE: usize = 16 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

/*#[pyfunction]
pub fn filter_annotated_file(
    filename: String, 
    tax_id: usize, 
    num_threads: usize)
{
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 128));
    let output_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(256));

    let generator_done = Arc::new(RwLock::new(false));

    let generator;
    let mut children = Vec::new();

    let main_thread_handle = thread::current();

    { // Explicit lifetime
        let generator_done = Arc::clone(&generator_done);
        let seq_queue = Arc::clone(&seq_queue);
        let output_queue = Arc::clone(&output_queue);
        let num_threads = num_threads.clone();

        generator =  match Builder::new()
                            .name("Generator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move|| 
                                fasta_generator(
                                    filename, 
                                    seq_queue, 
                                    generator_done, 
                                    num_threads,
                                    main_thread_handle,
                                ))
                                {
                                Ok(x) => x,
                                Err(y) => panic!("{}", y)
                            };
        
   }

   // Children to filter the FASTA sequences
   for _ in 0..num_threads {
        let seq_queue = Arc::clone(&seq_queue);
        let output_queue = Arc::clone(&output_queue);
        let tax_id = tax_id.clone();

        let child = match Builder::new()
                        .name("SequenceFilterWorker".into())
                        .stack_size(STACKSIZE)
                        .spawn(move || filter_sequence_worker(tax_id, seq_queue, output_queue)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

   // A single thread for the output

    { // explicit lifetime
        let output_queue = Arc::clone(&output_queue);
        let filename = format!("filtered_taxonomy_level_{}", tax_id.to_string());

        let outputchild = match Builder::new()
                            .name("OutputWorker".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move|| process_output(filename, output_queue)) {
                                Ok(x) => x,
                                Err(y) => panic!("{}", y)
                            };

        children.push(outputchild);
    }
   
    thread::park();

    let backoff = Backoff::new();

    while !*generator_done.read().unwrap() {
        thread::park();
        backoff.snooze();
    }

    generator.join().expect("Unable to join generator thread...");

    while !seq_queue.is_empty() {
        backoff.snooze();
    }

    output_queue.push(ThreadCommand::Terminate).expect("Unable to send terminate command to output buffer");

    for child in children {
        child.join().expect("Unable to join child thread");
    }

    println!("Finish filtering file!");

}*/

pub fn filter_annotated_file_singlethreaded(filename: String, tax_id: usize) {
    let f_output = File::create(format!("filtered_taxonomy_{}.sz", tax_id)).expect("Unable to write to file!");
    let out_buffer = BufWriter::with_capacity(16 * 1024 * 1024, f_output);
    let out_fh = snap::write::FrameEncoder::new(out_buffer);
    let mut out_fh = BufWriter::with_capacity(8 * 1024 * 1024, out_fh);

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(32 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(16 * 1024 * 1024, fasta);

    loop {
        let seq: Sequence = match bincode::deserialize_from(&mut reader) {
            Ok(x)   => x,
            Err(_)  => break
        };

        if seq.taxons.contains(&tax_id) {
            bincode::serialize_into(&mut out_fh, &seq).expect("Unable to write to bincode file");

        }
    }
}

fn open_file_with_progress_bar(filename: String) -> (Box<dyn Read>, ProgressBar)
{
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(64 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);
    return (Box::new(reader), pb)
}

fn snappy_output(filename: String) -> Box<dyn Write> {
    let buffer = BufWriter::with_capacity(64 * 1024 * 1024, 
        File::create(filename).expect("Unable to write to file!"));
    Box::new(BufWriter::with_capacity(16 * 1024 * 1024, snap::write::FrameEncoder::new(buffer)))
}

pub fn split_train_test_validation(
        filename: String,
        output_prefix: String, 
        test: f32, 
        validation: f32) 
{
    let mut train_out = snappy_output(format!("{}.train.sz", output_prefix));
    let mut test_out = snappy_output(format!("{}.test.sz", output_prefix));
    let mut validation_out = snappy_output(format!("{}.validation.sz", output_prefix));
    let (mut reader, mut pb) = open_file_with_progress_bar(filename);
    
    let mut rng = rand::thread_rng();

    loop {
        let seq: Sequence = match bincode::deserialize_from(&mut reader) {
            Ok(x)   => x,
            Err(_)  => break
        };

        if rng.gen::<f32>() < validation { 
            bincode::serialize_into(&mut validation_out, &seq).expect("Unable to write to bincode file");
            continue
        }
        if rng.gen::<f32>() < test { 
            bincode::serialize_into(&mut test_out, &seq).expect("Unable to write to bincode file");
            continue
        }

        bincode::serialize_into(&mut train_out, &seq).expect("Unable to write to bincode file");
    }
}

pub fn chunk_file(filename: String, output_filename: String, chunk_size: usize) {
    let f_output = File::create(format!("{}.sz", output_filename)).expect("Unable to write to file!");
    let out_buffer = BufWriter::with_capacity(64 * 1024 * 1024, f_output);
    let out_fh = snap::write::FrameEncoder::new(out_buffer);
    let mut out_fh = BufWriter::with_capacity(16 * 1024 * 1024, out_fh);

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(64 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);

    loop {
        let seq: Sequence = match bincode::deserialize_from(&mut reader) {
            Ok(x)   => x,
            Err(_)  => break
        };

        for chunk in seq.seq.chunks(chunk_size) {
            let newseq = Sequence { 
                id: seq.id.clone(),
                seq: chunk.to_vec(),
                taxons: seq.taxons.clone(),
                taxon: seq.taxon.clone()
            };
            bincode::serialize_into(&mut out_fh, &newseq).expect("Unable to write to bincode file");
        }
    }
}

pub fn shuffle_file(filename: String, output_filename: String, bins: usize) {

    let mut rng = rand::thread_rng();
    let x: u16 = rng.gen();
    let mut files: Vec<String> = Vec::with_capacity(bins);

    let mut shuffle_bins: Vec<_> = Vec::new();
    for i in 0..bins {
        let bin_filename = format!("temp_{}_shuffle_file_{}.bc", x, i);
        files.push(bin_filename.clone());
        let f_output = File::create(bin_filename).expect("Unable to write to file!");
        let out_buffer = BufWriter::with_capacity(4 * 1024 * 1024, f_output);
        shuffle_bins.push(out_buffer);
    }
    
    // Open file to be shuffled...
    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(64 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);

    loop {
        let seq: Sequence = match bincode::deserialize_from(&mut reader) {
            Ok(x)   => x,
            Err(_)  => break
        };

        let x = rng.gen_range(0, shuffle_bins.len());
                
        bincode::serialize_into(&mut shuffle_bins[x], &seq).expect("Unable to write to bincode file");
    }

    drop(shuffle_bins);
    files.shuffle(&mut rng);

    let out_buffer = BufWriter::with_capacity(64 * 1024 * 1024, 
        File::create(format!("{}.sz", output_filename)).expect("Unable to write to file!"));
    let mut shuffled_out = BufWriter::with_capacity(16 * 1024 * 1024, snap::write::FrameEncoder::new(out_buffer));

    for input_bin in files {
        let file = match File::open(&input_bin) {
            Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
            Ok(file) => file,
        };
        let mut reader = BufReader::with_capacity(32 * 1024 * 1024, file);

        loop {
            let seq: Sequence = match bincode::deserialize_from(&mut reader) {
                Ok(x)   => x,
                Err(_)  => break
            };
            bincode::serialize_into(&mut shuffled_out, &seq).expect("Unable to write to bincode file");
        }
        drop(reader);
        std::fs::remove_file(input_bin).expect("Unable to delete temp file");
    }
}

pub fn child_taxon_seqlengths(filename: String, parent_tax_id: usize) {

    let mut seqcounts: HashMap<usize, usize> = HashMap::with_capacity(64);
    let mut child_taxons: Vec<usize> = Vec::new();
    let mut names: HashMap<usize, String> = HashMap::new();

    let child_taxon_names = super::get_child_taxons_names(parent_tax_id);

    for (tax_id, tax_name) in child_taxon_names {
        let tax_name_lc = tax_name.to_lowercase();
        if tax_name_lc.contains("candidatus") ||
            tax_name_lc.contains("unclassified") ||
            tax_name_lc.contains("environmental") {
                continue;
        } else {
            // println!("{}: {}", tax_name, tax_id);
            names.insert(tax_id, tax_name);
            child_taxons.push(tax_id);
        }
    }

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(32 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(16 * 1024 * 1024, fasta);

    loop {
        let seq: Sequence = match bincode::deserialize_from(&mut reader) {
            Ok(x)   => x,
            Err(_)  => break
        };

        let mut foundone: bool = false;

        for child_taxon in &child_taxons {
            if seq.taxons.contains(&child_taxon) || seq.taxon == *child_taxon {
                if foundone {
                    panic!("Found one! But found another!");
                }
                let count = seqcounts.entry(child_taxon.clone()).or_insert(0);
                *count += seq.seq.len();
                foundone = true;
            }
        }
    }

    let f_output = File::create(format!("seqlength_count_parent_taxid_{}.sz", parent_tax_id)).expect("Unable to write to file!");
    let out_buffer = BufWriter::with_capacity(1 * 1024 * 1024, f_output);
    let out_fh = snap::write::FrameEncoder::new(out_buffer);
    let mut out_fh = BufWriter::with_capacity(1 * 1024 * 1024, out_fh);
    bincode::serialize_into(&mut out_fh, &seqcounts).expect("Unable to write to bincode file");

    pb.finish();

    for (key, val) in seqcounts {
        let name = names.get(&key).unwrap();
        println!("{}: {}", name, val);
    }
}


pub fn filter_annotated_file(filename: String, tax_id: usize, num_threads: usize) {
    let seq_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 128));
    let output_queue = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 1024));

    let threadpool = Arc::new(ThreadPool::new(num_threads));

    let generator;

    { // Explicit lifetime
        let seq_queue = Arc::clone(&seq_queue);
        let threadpool = Arc::clone(&threadpool);

        generator =  match Builder::new()
                            .name("SequenceGenerator".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move|| 
                                sequence_generator(
                                    filename, 
                                    seq_queue, 
                                    threadpool))
                                {
                                Ok(x) => x,
                                Err(y) => panic!("{}", y)
                            };
   }

   threadpool.children.write().unwrap().push(generator);

   let in_queue = seq_queue;
   let out_queue = output_queue;

   for _ in 0..num_threads {
        let work_fn = filter_annotated_sequences_fn(tax_id);
        let in_queue  = Arc::clone(&in_queue);
        let out_queue = Arc::clone(&out_queue);
        let child = match Builder::new().name("ThreadPoolWorker".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move|| {
                                queue_worker(work_fn, in_queue, out_queue);
                            }) {
                                Ok(x)  => x,
                                Err(y) => panic!("Unable to create thread: {}", y),
                            };
        threadpool.children.write().unwrap().push(child);
    }

    let writer_fn = write_output_fn_generator::<Sequence>(format!("filtered_taxonomy_{}", tax_id));

    let in_queue = Arc::clone(&out_queue);
    let child = match Builder::new().name("ThreadPoolWorker".to_string())
                        .stack_size(STACKSIZE)
                        .spawn(move|| {
                            output_worker(writer_fn, in_queue);
                        }) {
                            Ok(x)  => x,
                            Err(y) => panic!("Unable to create thread: {}", y),
                        };
    threadpool.children.write().unwrap().push(child);

    while !*threadpool.generator_done.read().unwrap() {
        std::thread::sleep(Duration::from_millis(1500));
        for child in &*threadpool.children.read().unwrap() {
            child.thread().unpark();
        }
    }

    terminate_queue(out_queue, 1);

    let threadpool = match Arc::try_unwrap(threadpool) {
        Ok(x) => x,
        Err(_) => panic!("Error"),
    };

    for child in threadpool.children.into_inner().unwrap() {
        child.join().expect("Unable to join worker thread!");
    }

    println!("Finished filtering file!");
}

fn process_output(
    filename: String, 
    output_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>) 
{
    let backoff = Backoff::new();

    let f_train = File::create(format!("{}.sz", filename)).expect("Unable to write to file");
    let out_buffer = BufWriter::with_capacity(16 * 1024 * 1024, f_train);
    let out_fh = snap::write::FrameEncoder::new(out_buffer);
    let mut out_fh = BufWriter::with_capacity(64 * 1024 * 1024, out_fh);

    loop {
        if let Ok(command) = output_queue.pop() {
            if let ThreadCommand::Terminate = command {
                return;
            }

            let seq = command.unwrap();
            bincode::serialize_into(&mut out_fh, &seq).expect("Unable to write to bincode file");
        } else {
            backoff.snooze();
        }
    }
}

fn split_process_output(
    filename: String, 
    output_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    validation: f32,
    test: f32,
) {
    let backoff = Backoff::new();

    let f_train = File::create(format!("{}.train.snappy", filename)).expect("Unable to write to file");
    let out_train_buffer = BufWriter::with_capacity(32 * 1024 * 1024, f_train);
    let mut out_train = snap::write::FrameEncoder::new(out_train_buffer);

    let f_validation = File::create(format!("{}.validation.snappy", filename)).expect("Unable to write to file");
    let out_validation_buffer = BufWriter::with_capacity(32 * 1024 * 1024, f_validation);
    let mut out_validation = snap::write::FrameEncoder::new(out_validation_buffer);

    let f_test = File::create(format!("{}.test.snappy", filename)).expect("Unable to write to file");
    let out_test_buffer = BufWriter::with_capacity(32 * 1024 * 1024, f_test);
    let mut out_test = snap::write::FrameEncoder::new(out_test_buffer);

    let mut rng = rand::thread_rng();

    loop {
        if let Ok(command) = output_queue.pop() {
            if let ThreadCommand::Terminate = command {
                return;
            }

            let seq = command.unwrap();

            if rng.gen::<f32>() < validation {
                out_validation.write_all(format!(">{}\n", &seq.id).as_bytes()).expect("Unable to write to file!");
                out_validation.write_all(&seq.seq).expect("Unable to write to file!");
                out_validation.write(b"\n").expect("Unable to write to file!");
            } else if rng.gen::<f32>() < test {
                out_test.write_all(format!(">{}\n", &seq.id).as_bytes()).expect("Unable to write to file!");
                out_test.write_all(&seq.seq).expect("Unable to write to file!");
                out_test.write(b"\n").expect("Unable to write to file!");
            } else {
                out_train.write_all(format!(">{}\n", &seq.id).as_bytes()).expect("Unable to write to file!");
                out_train.write_all(&seq.seq).expect("Unable to write to file!");
                out_train.write(b"\n").expect("Unable to write to file!");
            }
        } else {
            backoff.snooze();
            backoff.snooze();
        }
    }
}
/*
fn filter_sequence_worker(
        filter_tax_id: usize, 
        seq_queue: Arc<ArrayQueue<ThreadCommand<FastaSequence>>>,
        output_queue: Arc<ArrayQueue<ThreadCommand<FastaSequence>>>,)
{
    let mut contains_cache: HashMap<u32, bool> = HashMap::with_capacity(1024 * 8);
    let mut taxon_cache: HashMap<String, Option<super::Acc2TaxInner>> = HashMap::with_capacity(1024 * 8);
    let backoff = Backoff::new();

    loop {
        if let Ok(command) = seq_queue.pop() {
            if let ThreadCommand::Terminate = command {
                return;
            }

            let seq = command.unwrap();

            let short: String = parser::shorten(&seq.id);

            let map = match taxon_cache
                        .entry(short.clone())
                        .or_insert_with(|| 
                            super::load_taxon(&format!("acc2tax_db/{}.bc", short)))
                        {
                            Some(x) => x,
                            None    => continue
                        };

            // let tax_id = get_taxon(seq.id.clone());

            // println!("{}", seq.id);
            let accession = seq.id.split_ascii_whitespace().take(1).collect::<String>();

            let tax_id = match map.get(&accession) {
                Some(x) => *x,
                None    => continue
            };

            // println!("\n\n{}\n\n", tax_id);

            let contains = 
                contains_cache.entry(tax_id).or_insert_with(|| 
                    get_complete_taxonomy(tax_id as usize).contains(&filter_tax_id));

            //if complete.contains(&filter_tax_id) {
            if *contains { //|| rng.gen::<f32>() < other {
                let wp = ThreadCommand::Work(seq);

                let mut result = output_queue.push(wp);
                while let Err(PushError(wp)) = result {
                    backoff.snooze();
                    result = output_queue.push(wp);
                }
            }
        
            if contains_cache.len() > 1024 * 8 {
                contains_cache.clear();
		        contains_cache.shrink_to(1024 * 8);
            }

            if taxon_cache.len() > 1024 * 8 {
                taxon_cache.clear();
		        taxon_cache.shrink_to(1024 * 8);
            }

        } else {
            backoff.snooze();
            backoff.reset();
        }
    }
}*/

fn submit_to_queue<T>(queue: &Arc<ArrayQueue<ThreadCommand<T>>>, 
                      wp: ThreadCommand<T>) 
{
    let mut result = queue.push(wp);
    while let Err(PushError(wp)) = result {
        // println!("Queue is full: {:#?}", std::thread::current().name());
        result = queue.push(wp);
    }
}

fn queue_worker<T, W>(work_fn:   impl Fn(T) -> Option<ThreadCommand<W>>,
                      in_queue:  Arc<ArrayQueue<ThreadCommand<T>>>,
                      out_queue: Arc<ArrayQueue<ThreadCommand<W>>>,) 
{
    let backoff = Backoff::new();

    loop {
        if let Ok(command) = in_queue.pop() {
            match command {
                ThreadCommand::Terminate => return,
                ThreadCommand::Work(x)   => {
                    let output = work_fn(x);
                    match output {
                        Some(x) => submit_to_queue(&out_queue, x),
                        None    => ()
                    };
                }
            }
        } else {
            backoff.snooze();
/*            eventual_sleep += 1;
            if eventual_sleep == 1024 * 1024 * 1024 {
                thread::park();
                eventual_sleep = 0;
            } */
        }
    }
}

fn output_worker<T>(mut work_fn:   impl FnMut(T),
                    in_queue:  Arc<ArrayQueue<ThreadCommand<T>>>) 
{
    let backoff = Backoff::new();

    loop {
        if let Ok(command) = in_queue.pop() {
            match command {
                ThreadCommand::Terminate => return,
                ThreadCommand::Work(x)   => {
                    backoff.reset();
                    work_fn(x);
                }
            }
        } else {
            // We are the final destination, never park
            backoff.snooze();
        }
    }
}

fn filter_annotated_sequences_fn(tax_id: usize) 
        -> impl Fn(Sequence) -> Option<ThreadCommand<Sequence>> {
    move |seq: Sequence| {
        if seq.taxons.contains(&tax_id) {
            return Some(ThreadCommand::Work(seq))
        } 
        return None
    }
}

fn write_output_fn_generator<T>(filename: String)
        -> impl FnMut(T)
    where T: serde::Serialize
{
    let f_output = File::create(format!("{}.sz", filename)).expect("Unable to write to file!");
    let out_buffer = BufWriter::with_capacity(64 * 1024 * 1024, f_output);
    let out_fh = snap::write::FrameEncoder::new(out_buffer);
    let mut out_fh = BufWriter::with_capacity(64 * 1024 * 1024, out_fh);

    move |seq: T| {
        bincode::serialize_into(&mut out_fh, &seq).expect("Unable to write to bincode file");
    }
}

fn terminate_queue<T>(queue: Arc<ArrayQueue<ThreadCommand<T>>>, 
                      count: usize) 
{
    for _ in 0..count {
        let mut result = queue.push(ThreadCommand::<T>::Terminate);

        while let Err(PushError(wp)) = result {
            result = queue.push(wp);
        }
    }
}

#[pyfunction]
pub fn convert_ntfasta_file(
    filename: String, 
    output: String,
    num_threads: usize)
// Convert file to bincode/snappy for faster processing
// Stores accession/taxon information inside the Sequence struct
{
    let seq_queue      = Arc::new(ArrayQueue::<ThreadCommand<FastaSequence>>::new(64 * 1024));
    let output_queue   = Arc::new(ArrayQueue::<ThreadCommand<Sequence>>::new(1024 * 1024));

    let generator_done = Arc::new(RwLock::new(false));

    let generator;
    let mut children = Vec::new();

    let main_thread_handle = thread::current();

    { // Explicit lifetime
        let generator_done = Arc::clone(&generator_done);
        let seq_queue = Arc::clone(&seq_queue);
        let num_threads = num_threads.clone();

        generator =  match Builder::new()
                            .name("FastaGeneratorWorker".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move|| 
                                fasta_generator(
                                    filename, 
                                    seq_queue, 
                                    generator_done, 
                                    num_threads,
                                    main_thread_handle,
                                ))
                                {
                                Ok(x) => x,
                                Err(y) => panic!("{}", y)
                            };
        
   }

   // Children to filter the FASTA sequences
   for _ in 0..num_threads {
        let seq_queue = Arc::clone(&seq_queue);
        let output_queue = Arc::clone(&output_queue);

        let child = match Builder::new()
                        .name("AnnotateSequenceWorker".into())
                        .stack_size(STACKSIZE)
                        .spawn(move || annotate_sequences(seq_queue, output_queue)) {
                            Ok(x)  => x,
                            Err(y) => panic!("{}", y)
                        };
        
        children.push(child);
    }

   // A single thread for the output

    { // explicit lifetime
        let output_queue = Arc::clone(&output_queue);

        let outputchild = match Builder::new()
                            .name("OutputWorker".to_string())
                            .stack_size(STACKSIZE)
                            .spawn(move|| process_output(output, output_queue)) {
                                Ok(x) => x,
                                Err(y) => panic!("{}", y)
                            };

        children.push(outputchild);
    }
   
    thread::park();

    let backoff = Backoff::new();

    while !*generator_done.read().unwrap() {
        thread::park();
        backoff.snooze();
    }

    generator.join().expect("Unable to join generator thread...");

    while !seq_queue.is_empty() {
        backoff.snooze();
    }

    output_queue.push(ThreadCommand::Terminate).expect("Unable to send terminate command to output buffer");

    for child in children {
        child.join().expect("Unable to join child thread");
    }

    println!("Finish filtering file!");

}

pub fn annotate_sequences(
    in_queue:  Arc<ArrayQueue<ThreadCommand<FastaSequence>>>, 
    out_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,) 
{
    let backoff = Backoff::new();

    let mut taxon_cache: HashMap<String, Option<super::Acc2TaxInner>> = HashMap::with_capacity(1024 * 16);

    loop {
        if let Ok(command) = in_queue.pop() {
            if let ThreadCommand::Terminate = command {
                return;
            }

            let seq = command.unwrap();

            let short: String = parser::shorten(&seq.id);

            let map = match taxon_cache
                        .entry(short.clone())
                        .or_insert_with(|| 
                            super::load_taxon(&format!("acc2tax_db/{}.bc", short)))
                        {
                            Some(x) => x,
                            None    => continue
                        };

            let accession = seq.id.split_ascii_whitespace().take(1).collect::<String>();

            let tax_id = match map.get(&accession) {
                Some(x) => *x as usize,
                None    => continue
            };

            let complete = get_complete_taxonomy(tax_id as usize);

            let wp = ThreadCommand::Work(
                Sequence { seq: seq.seq, id: seq.id, taxons: complete, taxon: tax_id });

            submit_to_queue(&out_queue, wp);

            if taxon_cache.len() > 1024 * 16 {
                taxon_cache.clear();
		        taxon_cache.shrink_to(1024 * 16);
            }

        } else {
            backoff.snooze();
        }
    }
}

pub fn fasta_generator(
    filename: String, 
    out_queue: Arc<ArrayQueue<ThreadCommand<FastaSequence>>>, 
    generator_done: Arc<RwLock<bool>>,
    num_threads: usize,
    main_thread_handle: Thread,) 
{

let file = match File::open(&filename) {
    Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
    Ok(file) => file,
};

let pb = ProgressBar::new(file.metadata().unwrap().len());
pb.set_style(ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                .progress_chars("█▇▆▅▄▃▂▁  "));

let mut buffer: Vec<u8> = Vec::with_capacity(1024);
let mut id: String = String::from("INVALID_ID_FIRST_ENTRY_YOU_SHOULD_NOT_SEE_THIS");
let mut seqbuffer: Vec<u8> = Vec::with_capacity(16 * 1024 * 1024); // 8 Mb to start, will likely increase...
let mut seqlen: usize = 0;

let file = BufReader::with_capacity(32 * 1024 * 1024, pb.wrap_read(file));

let fasta: Box<dyn Read> = if filename.ends_with("gz") {
    Box::new(flate2::read::GzDecoder::new(file))
} else if filename.ends_with("snappy") || filename.ends_with("sz") {
    Box::new(snap::read::FrameDecoder::new(file))
} else {
    Box::new(file)
};

let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);

let backoff = Backoff::new();

while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {

    if bytes_read == 0 { 
        // No more reads, thus no more data...
        // Submit the last sequence (or in the case of some genomes, the entire sequence)
        let wp = ThreadCommand::Work(FastaSequence { seq: seqbuffer[..seqlen].to_vec(), id: id });
        seqbuffer.clear();
        submit_to_queue(&out_queue, wp);

        // Slight delay before shutting down...
        for _ in 0..5 {
            backoff.snooze();
        }

        break;
    }

    match buffer[0] {
        // 62 is a > meaning we have a new sequence id.
        62 => {
            let wp = ThreadCommand::Work(FastaSequence { seq: seqbuffer[..seqlen].to_vec(), id: id });
            seqbuffer.clear();
            seqlen = 0;

            submit_to_queue(&out_queue, wp);

            let slice_end = bytes_read.saturating_sub(1);
            id = String::from_utf8(buffer[1..slice_end].to_vec()).expect("Invalid UTF-8 encoding...");
            // pb.set_message(&format!("Seq Queue: {} Output: {}", seq_queue.len(), output_queue.len()));
        },
        _  => {
            let slice_end = bytes_read.saturating_sub(1);
            seqbuffer.extend_from_slice(&buffer[0..slice_end]);
            seqlen = seqlen.saturating_add(slice_end);
        }
    }

    buffer.clear();
    }

    *generator_done.write().unwrap() = true;
    main_thread_handle.unpark();

    terminate_queue(out_queue, num_threads);
}

pub fn sequence_generator(
    filename: String,
    out_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    threadpool: Arc<ThreadPool>)
{

    let file = match File::open(&filename) {
        Err(why) => panic!("Couldn't open {}: {}", filename, why.to_string()),
        Ok(file) => file,
    };

    let pb = ProgressBar::new(file.metadata().unwrap().len());
    pb.set_style(ProgressStyle::default_bar()
                    .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>5}/{len:5} {eta_precise} {msg}")
                    .progress_chars("█▇▆▅▄▃▂▁  "));

    let file = BufReader::with_capacity(32 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else if filename.ends_with("snappy") || filename.ends_with("sz") {
        Box::new(snap::read::FrameDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(32 * 1024 * 1024, fasta);

    while let Ok(seq) = bincode::deserialize_from(&mut reader) {
        let wp = ThreadCommand::Work(seq);
        submit_to_queue(&out_queue, wp);
        // pb.set_message(&format!("Queue Size: {}", out_queue.len()));
    }

    *threadpool.generator_done.write().unwrap() = true;
    threadpool.main_thread_handle.unpark();

    terminate_queue(out_queue, threadpool.num_threads);
}