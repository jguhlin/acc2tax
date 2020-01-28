use std::thread::{Builder, Thread};
use std::thread;

use std::io::prelude::*;
use std::io::BufWriter;

use std::sync::{Arc, RwLock};

use std::fs::File;
use std::io::{BufReader, Read, BufRead};

use crate::parser;

use std::collections::HashMap;

use crossbeam::queue::{ArrayQueue, PushError};
use crossbeam::utils::Backoff;

use pyo3::prelude::*;

use rand::Rng;

use indicatif::ProgressBar;
use indicatif::ProgressStyle;

use super::{get_taxon, get_complete_taxonomy};

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

#[derive(PartialEq)]
pub struct Sequence {
    pub seq: Vec<u8>,
    pub id:     String,
}

const STACKSIZE: usize = 16 * 1024 * 1024;  // Stack size (needs to be > BUFSIZE + SEQBUFSIZE)

// And unpark it from the generator...
#[pyfunction]
pub fn filter_fasta_file(
    filename: String, 
    tax_id: usize, 
    num_threads: usize, 
    validation: f32, // Fraction to put in the validation file
    test: f32, // Fraction to put in the test file
    other: f32) // Chance for a non-matching sequence to end up in one of the files...
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
                                sequence_generator(
                                    filename, 
                                    seq_queue, 
                                    output_queue,
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
                        .spawn(move || filter_sequence_worker(tax_id, seq_queue, output_queue, other)) {
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
                            .spawn(move|| process_output(filename, output_queue, validation, test)) {
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

fn process_output(
    filename: String, 
    output_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
    validation: f32,
    test: f32,
) {
    let backoff = Backoff::new();

    let f_train = File::create(format!("{}.train.snappy", filename)).expect("Unable to write to file");
    let out_train_buffer = BufWriter::with_capacity(32 * 1024 * 1024, f_train);
    let mut out_train = snap::Writer::new(out_train_buffer);

    let f_validation = File::create(format!("{}.validation.snappy", filename)).expect("Unable to write to file");
    let out_validation_buffer = BufWriter::with_capacity(32 * 1024 * 1024, f_validation);
    let mut out_validation = snap::Writer::new(out_validation_buffer);

    let f_test = File::create(format!("{}.test.snappy", filename)).expect("Unable to write to file");
    let out_test_buffer = BufWriter::with_capacity(32 * 1024 * 1024, f_test);
    let mut out_test = snap::Writer::new(out_test_buffer);

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

fn filter_sequence_worker(
        filter_tax_id: usize, 
        seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
        output_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>,
        other: f32)
{
    let mut contains_cache: HashMap<u32, bool> = HashMap::with_capacity(1024 * 8);
    let mut taxon_cache: HashMap<String, Option<super::Acc2TaxInner>> = HashMap::with_capacity(1024 * 8);
    let backoff = Backoff::new();

    let mut rng = rand::thread_rng();

    
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
            if *contains || rng.gen::<f32>() < other {
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
}

fn sequence_generator(
        filename: String, 
        seq_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>, 
        output_queue: Arc<ArrayQueue<ThreadCommand<Sequence>>>, 
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

    let file = BufReader::with_capacity(128 * 1024 * 1024, pb.wrap_read(file));

    let fasta: Box<dyn Read> = if filename.ends_with("gz") {
        Box::new(flate2::read::GzDecoder::new(file))
    } else {
        Box::new(file)
    };

    let mut reader = BufReader::with_capacity(64 * 1024 * 1024, fasta);

    let backoff = Backoff::new();

    while let Ok(bytes_read) = reader.read_until(b'\n', &mut buffer) {

        if bytes_read == 0 { 
            // No more reads, thus no more data...
            // Submit the last sequence (or in the case of some genomes, the entire sequence)
            let wp = ThreadCommand::Work(Sequence { seq: seqbuffer[..seqlen].to_vec(), id: id });
            seqbuffer.clear();

            let mut result = seq_queue.push(wp);
            while let Err(PushError(wp)) = result {
                result = seq_queue.push(wp);
            }

            backoff.spin();
            backoff.spin();
            backoff.spin();
            backoff.spin(); // Very slight delay before leaving...

            break;
        }

        match buffer[0] {
            // 62 is a > meaning we have a new sequence id.
            62 => {
                let wp = ThreadCommand::Work(Sequence { seq: seqbuffer[..seqlen].to_vec(), id: id });
                seqbuffer.clear();
                seqlen = 0;

                let mut result = seq_queue.push(wp);
                while let Err(PushError(wp)) = result {
                    result = seq_queue.push(wp);
                }

                let slice_end = bytes_read.saturating_sub(1);
                id = String::from_utf8(buffer[1..slice_end].to_vec()).expect("Invalid UTF-8 encoding...");
                pb.set_message(&format!("Seq Queue: {} Output: {}", seq_queue.len(), output_queue.len()));
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

    for _ in 0..num_threads {
        let mut result = seq_queue.push(ThreadCommand::Terminate);

        while let Err(PushError(wp)) = result {
            result = seq_queue.push(wp);
        }
    }
}
