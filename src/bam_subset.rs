use rust_htslib::bam::{Record};
use std::collections::HashSet;
use std::collections::HashMap;
use rust_htslib::bam::record::Aux;

use std::thread;
use rayon::prelude::*;
use rayon::slice::ParallelSlice;


pub struct BamSubset {
    records: HashMap<String, Vec::<Record>>,
    cells: HashSet<String>,
    num_threads: usize,
}

impl BamSubset{
    
    pub fn new(num_threads: usize) -> Self {
        Self{
            records: HashMap::new(),
            cells: HashSet::<String>::new(),
            num_threads,
        }
    }

    pub fn push(&mut self, record:Record ) {

        if let Some(tag_value) = record.aux(b"CB").ok() {
            if let Aux::String(barcode) = tag_value {
                // Convert the barcode to a String and do something with it
                let record_barcode = barcode.to_string();
                // Now you can use record_barcode
                self.records.entry(record_barcode.clone())
                    .or_insert_with(Vec::new)
                    .push(record);
                self.cells.insert(record_barcode);
            }
        }
    }

    pub fn len( &self) -> usize{
        self.cells.len()
    }

    pub fn assemble_contigs( &self ) -> Vec<String> {

        // Collect keys from the HashMap into a Vec
        let keys: Vec<String> = self.records.keys().cloned().collect();

        // Use par_chunks to process the keys in parallel
        let consensus_sequence = keys.par_chunks(self.num_threads) // Adjust chunk size as needed
            .map(|key_chunk| {
                let mut local_consensus = Vec::new();
                
                for key in key_chunk {
                    // Use the key to access records from the HashMap
                    if let Some(records) = self.records.get(key) {
                        let mut contig_reads = Vec::new();
                        let mut current_contig_start = 0;
                        let mut current_contig_end = 0;
                        let mut id = 1;

                        for record in records {
                            let start_pos = record.pos();
                            let end_pos = record.cigar().end_pos();

                            if start_pos <= current_contig_end {
                                contig_reads.push(record.clone());
                                if end_pos > current_contig_end {
                                    current_contig_end = end_pos;
                                }
                            } else {
                                if !contig_reads.is_empty() {
                                    let consensus = self.build_consensus(&contig_reads);
                                    local_consensus.push(format!(">{}|transcript_{}\n{}", key, id, consensus));
                                    id += 1;
                                }

                                contig_reads.clear();
                                contig_reads.push(record.clone());
                                current_contig_start = start_pos;
                                current_contig_end = end_pos;
                            }
                        }

                        if !contig_reads.is_empty() {
                            let consensus = self.build_consensus(&contig_reads);
                            local_consensus.push(format!(">{}|transcript_{}\n{}", key, id, consensus));
                        }
                    }
                }
                
                local_consensus
            })
           .flatten() // Flatten the results into a single iterator
            .collect(); // Collect into the final Vec

        consensus_sequence
    }

    // Helper function to build consensus from a set of overlapping reads
    pub fn build_consensus( &self, contig_reads:&Vec<Record> ) -> String {
        // 1. Initialize a vector to hold consensus bases for each position
        let mut consensus_bases: Vec<char> = Vec::new();

        // 2. Determine the range of positions covered by the reads
        let start_pos = contig_reads[0].pos();
        let end_pos = contig_reads.iter().map(|r| r.cigar().end_pos()).max().unwrap();

        // 3. Iterate through each position within this range
        for pos in start_pos..end_pos {
            let mut base_counts = std::collections::HashMap::new();

            // 4. For each read that covers this position, count the occurrences of each base
            for read in contig_reads {

                if pos >= read.pos() && pos < read.cigar().end_pos() {
                    // Convert position to index within the read's sequence
                    let idx = (pos - read.pos()) as usize;
                    // Access the base at this position
                    let array = read.seq().as_bytes();
                    if array.len() > idx{
                        let base = read.seq().as_bytes()[idx] as char;
                        *base_counts.entry(base).or_insert(0) += 1;
                    }
                }
            }

            // 5. Determine the most common base at this position (consensus)
            if let Some((&consensus_base, _)) = base_counts.iter().max_by_key(|&(_, count)| count) {
                consensus_bases.push(consensus_base);
            }
        }
        // 6. Convert the vector of consensus bases to a string and return it
        consensus_bases.iter().collect()
    }
}

