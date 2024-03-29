use super::*;
use indicatif::ProgressBar;
use indicatif::ProgressStyle;
use zerocopy::{AsBytes, U32};
use redb::{Database, Error, ReadableTable, TableDefinition, Table};

use std::collections::HashSet;

enum ThreadCommand<T> {
    Work(T),
    Terminate,
}

impl ThreadCommand<Vec<Vec<u8>>> {
    // Consumes the ThreadCommand, which is just fine...
    fn unwrap(self) -> Vec<Vec<u8>> {
        match self {
            ThreadCommand::Work(x) => x,
            ThreadCommand::Terminate => panic!("Unable to unwrap terminate command"),
        }
    }
}

// const SHORT_LENGTH: usize = 8;

/* pub fn shorten(acc: &str) -> String {
    let length = acc.len();
    let max = if length < SHORT_LENGTH {
        length
    } else {
        SHORT_LENGTH
    };

    acc[..max].to_string()
} */

fn parse_line(line: &[u8]) -> Result {
    let data = unsafe {
        std::str::from_utf8_unchecked(line)
            .split_ascii_whitespace()
            .skip(1)
            .take(2)
            .collect::<Vec<&str>>()
    };
    let taxon = data[1].parse::<u32>().unwrap();
    // let acc: Vec<u8> = data[0].as_bytes().to_vec();
    let acc: String = data[0].to_string();

    (acc, U32::<BigEndian>::new(taxon).get())
}

pub fn read_taxonomy(
    num_threads: usize,
    a2tdb: Arc<Database>,
    acc2tax_filename: String,
    nodes_filename: String,
    names_filename: String,
) -> (Vec<String>, Vec<u32>, Vec<usize>, Vec<String>) {
    let names_child = match Builder::new()
        .name("ParseNames".into())
        .spawn(move || parse_names(names_filename))
    {
        Ok(x) => x,
        Err(y) => panic!("{}", y),
    };

    let nodes_child = match Builder::new()
        .name("ParseNodes".into())
        .spawn(move || parse_nodes(nodes_filename))
    {
        Ok(x) => x,
        Err(y) => panic!("{}", y),
    };

    let gb2accession_fh = File::open(acc2tax_filename).unwrap();

    let pb = ProgressBar::new(gb2accession_fh.metadata().unwrap().len());
    pb.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {eta} {msg}")
            .expect("Unable to set progress bar style")
            .progress_chars("█▇▆▅▄▃▂▁  "),
    );

    let gb2accession_fh = BufReader::new(pb.wrap_read(gb2accession_fh));
    let gb2accession = flate2::read::GzDecoder::new(gb2accession_fh);
    let gb2accession = BufReader::new(gb2accession);
    let backoff = Backoff::new();

    // let mut children = Vec::new();
    //let queue = Arc::new(ArrayQueue::<ThreadCommand<Vec<Vec<u8>>>>::new(2048));

    let jobs = Arc::new(AtomicCell::new(0_usize));

/*
    for _ in 0..num_threads {
        let queue = Arc::clone(&queue);
        let jobs = Arc::clone(&jobs);
        let a2tdb = Arc::clone(&a2tdb);

        let child = match Builder::new()
            .name("TaxonReader".into())
            .spawn(move || _worker_thread(queue, a2tdb, jobs))
        {
            Ok(x) => x,
            Err(y) => panic!("{}", y),
        };
        children.push(child);
    } */

    jobs.fetch_add(1); // So the merger thread doesn't stop right away...

    let mut lines = 0;

    let write_txn = a2tdb.begin_write().expect("Unable to begin read transaction");
    let mut table = write_txn.open_table(TABLE).expect("Unable to open table");

    for (id, val) in gb2accession
        .byte_lines()
        .into_iter()
        .skip(1)
        // .chunks(512 * 1024)
        .map(|x| parse_line(&x.unwrap()))
    {

        // lines += 1;
        // pb.set_message(format!("{} lines", lines));

        // jobs.fetch_add(1);

        table.insert(id.as_str(), val).expect("Unable to insert");
    }

    drop(table);
    write_txn.commit().expect("Unable to commit transaction");

    pb.set_message("Completed processing all lines...");

    jobs.fetch_sub(1); // Merger thread extra...

    while jobs.load() > 0 {
        pb.set_message(format!("{} jobs remaining", jobs.load()));
        backoff.snooze();
    }

    pb.set_message("All jobs finished...");

    let names = names_child
        .join()
        .expect("Unable to join taxonomy names thread");
    let (taxon_to_parent, taxon_rank) = nodes_child
        .join()
        .expect("Unable to join taxonomy nodes thread");

    pb.finish_with_message("Complete");

    (names.0, names.1, taxon_to_parent, taxon_rank)
}

fn _worker_thread(
    queue: Arc<ArrayQueue<ThreadCommand<Vec<Vec<u8>>>>>,
    a2tdb: Arc<Database>,
    jobs: Arc<AtomicCell<usize>>,
) {
    let backoff = Backoff::new();

    loop {
        if let Some(command) = queue.pop() {
            if let ThreadCommand::Terminate = command {
                break;
            }

            let lines = command.unwrap();

            let result = lines
                .into_iter()
                .map(|x| parse_line(&x))
                .collect::<Vec<Result>>();

            let write_txn = a2tdb.begin_write().expect("Unable to begin read transaction");
            let mut table = write_txn.open_table(TABLE).expect("Unable to open table");

            for (x, y) in result {
                table.insert(x.as_str(), y).expect("Unable to insert");
            }
            drop(table);
            write_txn.commit().expect("Unable to commit transaction");

            jobs.fetch_sub(1);
        } else {
            backoff.snooze();
        }
    }
}

pub fn parse_names(filename: String) -> (Vec<String>, Vec<u32>) {
    let mut names: Vec<String> = Vec::with_capacity(3_006_098);

    let reader = BufReader::new(File::open(filename).expect("Unable to open taxonomy names file"));

    let lines = reader.lines();
    let mut taxids = HashSet::new();

    for line in lines {
        let split = line
            .expect("Error reading line")
            .split('|')
            .map(|x| x.trim().to_string())
            .collect::<Vec<String>>();

        let id: usize = split[0].parse().expect("Error converting to number");
        let name: &str = &split[1];
        let class: &str = &split[3];

        match names.get(id) {
            None => {
                names.resize(id + 1, "".to_string());
                names[id] = name.into();
                taxids.insert(id as u32);
            }
            Some(_) => {
                if class == "scientific name" {
                    names[id] = name.into();
                }
            }
        };
    }

    let mut taxids = taxids.into_iter().collect::<Vec<u32>>();
    taxids.sort_unstable();

    (names, taxids)
}

pub fn parse_nodes(filename: String) -> (Vec<usize>, Vec<String>) {
    let mut taxon_to_parent: Vec<usize> = Vec::with_capacity(3_000_000);
    let mut taxon_rank: Vec<String> = Vec::with_capacity(3_000_000);

    let reader = BufReader::new(File::open(filename).expect("Unable to open taxonomy names file"));

    let lines = reader.lines();

    for line in lines {
        let split = line
            .expect("Error reading line")
            .split('|')
            .map(|x| x.trim().to_string())
            .collect::<Vec<String>>();

        let tax_id: usize = split[0].parse().expect("Error converting to number");
        let parent_id: usize = split[1].parse().expect("Error converting to number");
        let rank: &str = &split[2];

        if taxon_to_parent.len() < tax_id + 1 {
            taxon_to_parent.resize(tax_id + 100, 0);
            taxon_rank.resize(tax_id + 100, "".to_string());
        }

        taxon_to_parent[tax_id] = parent_id;
        taxon_rank[tax_id] = rank.into();
    }

    taxon_to_parent.shrink_to_fit();
    taxon_rank.shrink_to_fit();

    (taxon_to_parent, taxon_rank)
}

// Unit tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_line() {
        let line: Vec<u8> = "X59856  X59856.2        9913    109659794"
            .as_bytes()
            .to_vec();
        assert_eq!(parse_line(&line), ("X59856.2".to_string(), 9913_u32));
    }
}
