// In machine learning balanced classes are very important
// as well as shuffled input data
// This does both.

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

use crate::fasta::{ThreadCommand, Sequence, sequence_generator, STACKSIZE, filter_sequence_worker};

// And unpark it from the generator...
#[pyfunction]
pub fn balance_and_shuffle(
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