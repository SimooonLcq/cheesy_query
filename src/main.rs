use std::collections::HashSet;
use bio::io::fasta::Reader;
use clap::Parser;
use serde::Deserialize;
use serde_json;
use std::fs;

#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Args {
    ///Kmer length
    #[arg(short)]
    k: usize,

    ///Query Fasta filename
    #[arg(short = 'q', long = "query")]
    query: String,

    ///Index Json filename
    #[arg(short = 'i', long = "index")]
    index: String,
}

#[derive(Deserialize)]
struct HSet{
    x: HashSet<Vec<u8>>,
}

fn main() {
    let parser = Args::parse();

    let qfile = Reader::from_file(parser.query).expect("Error while reading query");
    let json = fs::read_to_string(parser.index).expect("Error while reading index");
    let deserialized: HSet = serde_json::from_str(&json).unwrap();
    let kmers: HashSet<Vec<u8>> = deserialized.x;

    let mut count_pos = 0;
    let mut count_neg = 0;
    for record in qfile.records(){
        let seq = record.expect("Error during fasta record parsing").seq().to_vec();
        for i in 0..seq.len()-parser.k+1{
            let kmer = (&seq[i..i+parser.k]).to_vec();
            let rc = revcomp(&kmer);
            if kmers.contains(&kmer) || kmers.contains(&rc) {
                count_pos += 1;
            }
            else {
                count_neg += 1;
            }
        }
    }
    println!("{:?} positive kmers", count_pos);
    println!("{:?} negative kmers", count_neg);
}

fn revcomp(kmer:&Vec<u8>)->Vec<u8> {
    let mut rc:Vec<u8> = Vec::new();
    for nuc in kmer{
        match nuc{
            65 => rc.push(84),
            84 => rc.push(65),
            67 => rc.push(71),
            71 => rc.push(67),
            _ => println!("Wrong character in a kmer, correct it and start again"),
        }
    }
    rc.reverse();
    return rc;
}