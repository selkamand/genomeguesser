use clap::{command, Arg};
use core::fmt;
use std::string::String;
use colored::Colorize;
use csv::StringRecord;

//type VarRecord = (String, usize, String);
type VarRecord = std::collections::HashMap<String, String>;

#[derive(Clone)]
enum Coordinate {
    Zero,   // Zero-Based coordinate system
    One,    // One-based coordinate system
    Unsure, // Unsure
}

const SHRUGGING_MAN: &str = "¯\\_(ツ)_/¯";
const EQUALS_HEADER: &str = "================================";

impl fmt::Display for Coordinate {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Coordinate::Zero => write!(f, "{}", "0-based ✅".bold().green()),
            Coordinate::One => write!(f, "{}", "1-based ✅".bold().green()),
            Coordinate::Unsure => write!(f, "{}" ,format!("{} {}",SHRUGGING_MAN.bold(), "(Unclear)").yellow()),
        }
    }
}

#[derive(Clone)]
enum RefGenomeMatches {
    Yes,    // Zero-Based coordinate system
    No,     // One-based coordinate system
    Unsure, // Unsure
}

impl fmt::Display for RefGenomeMatches {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            RefGenomeMatches::Yes => write!(f, "{}", "Perfectly aligns ✅".bold().green()),
            RefGenomeMatches::No => write!(f, "{}", "No, Ref mismatches found ❌".bold().red()),
            RefGenomeMatches::Unsure => write!(f, "{}", format!("{} {}",SHRUGGING_MAN.bold(), "(Unclear)").yellow()),
        }
    }
}

fn assess_ref_genomes_support(path_genomes: Vec<String>, path_variants: &str) -> String {

    
    let results: Vec<(String, RefGenomeMatches, Coordinate)> = path_genomes.iter().map(|genome| { assess_ref_genome_support(path_variants, genome) }).collect();
    
    let genome_matches: Vec<RefGenomeMatches> = results.iter().map(|tup| tup.1.clone() ).collect();
    let results_coord: Vec<Coordinate> = results.iter().map(|tup: &(String, RefGenomeMatches, Coordinate)| tup.2.clone() ).collect();


    let ref_genome_matched: Vec<bool> = genome_matches.iter().map(|ref_match| match ref_match{
        RefGenomeMatches::Yes => true,
        RefGenomeMatches::No => false,
        RefGenomeMatches::Unsure => false}
    ).collect();

    let num_matches = ref_genome_matched.iter().filter(|x| **x).count();
    let mut final_result = "Unsure";

    if num_matches == 0{
        final_result = "no matches";
    }
    else if num_matches == 1{
        final_result = "1 match";
    }
    else if num_matches > 1{
        final_result = "Ambiguous, multiple matches";
    }

    // Collect Genome Names
    let matched_genome_name: String;
    let matched_coord: String;
    let genome_names: Vec<String> = path_genomes.iter().map(|genome_path| {
        std::path::Path::new(genome_path.to_string().as_str())
            .file_name().map(|name| name.to_string_lossy().into_owned())
            .unwrap_or_else(|| panic!("Failed to parse genome path as filepath"))
    }).collect();
    
    // 
    if final_result == "1 match"{
        let index_of_matched = ref_genome_matched.iter().enumerate()
        .find_map(|(i, &value)| if value { Some(i) } else { None })
        .expect("No match found, but expected one match");

    // Access the matched genome name safely
    matched_genome_name = genome_names.get(index_of_matched).unwrap().to_string();
    matched_coord = results_coord.get(index_of_matched).unwrap().to_string();
    

    }
    else if final_result == "Ambiguous, multiple matches" {
        matched_genome_name = r"¯\_(ツ)_/¯ (Variant file works equally well for multiple variants)".to_string();
        matched_coord = Coordinate::Unsure.to_string();
    }
    else if final_result == "no matches" {
        matched_genome_name = r"¯\_(ツ)_/¯ (No unambiguous matches to any reference genome)".to_string();
        matched_coord = Coordinate::Unsure.to_string();
    }
    else  {
        panic!("Should never see this. If you are developer has added 'final_result' possibilities and not accounted for them")
    }

    //Print Details
    println!("\n\n{EQUALS_HEADER}\n{}\n{EQUALS_HEADER}", "Final Summary".bold().blue());
    println!("Matched reference genome: {} ✅", matched_genome_name.bold().bright_magenta());
    println!("Coordinate system genome: {}", matched_coord.bold().bright_magenta());

    // Return Result
    matched_genome_name.to_string()

}

fn main() {
    let matches = command!()
        .about(
            "Checks variant file aligns with reference genome. Also tells you if it is 0 or 1based?",
        )
        .arg(
            Arg::new("variants")
                .required(true)
                .help("TSV file with Variants"),
        )
        .arg(
            Arg::new("genome")
                .required(true)
                .value_delimiter(',')
                .help("FASTA file with genome to test"),
        )
        .get_matches();

    let path_variants = matches.get_one::<String>("variants").unwrap();
    let path_genomes: Vec<String> = matches.get_many::<String>("genome").unwrap().cloned().collect();

    assess_ref_genomes_support(path_genomes, path_variants);

}  

fn assess_ref_genome_support(path_variants: &str, path_genome: &str) -> (String, RefGenomeMatches, Coordinate) {
    
    // Create Reader
    let reader = rust_htslib::faidx::Reader::from_path(path_genome).unwrap_or_else(|_err| {
        panic!("Failed to read genome file ({}). Missing file", path_genome)
    });

    // Pull all valid chromosomes
    let valid_chromosomes = get_all_chromosomes(&reader);
    
    // Read TSV file
    // Create Reader
    let mut var_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path_variants)
        .unwrap_or_else(|_err| {
            panic!("Failed to read TSV file ({}). Missing file?", path_variants)
        });

    // Allow users to add column names (Will need to take from user input)
    // Note by default values will be None and the default valid names above will be used
    let user_specified_chrom = None;
    let user_specified_pos = None;
    let user_specified_ref = None;

    // Use user supplied colnames unless they are None
    let chrom_column_names = get_column_names(
        user_specified_chrom,
        vec!["chrom", "Chrom", "CHROM", "chromosome", "Chromosome"],
    );
    let pos_column_names = get_column_names(
        user_specified_pos,
        vec!["pos", "Pos", "POS", "position", "Position"],
    );
    let ref_column_names = get_column_names(
        user_specified_ref,
        vec!["ref", "REF", "Ref", "Reference", "Reference_Allele"],
    );

    // Ensure File Has Headers
    if !var_reader.has_headers() {
        panic!("genomeguesser only works if the variant file has a header row. Please add a header and try again")
    }

    // Figure Out Column Names
    let headers = var_reader.headers().expect("Failed to find headers");
    let col_chrom = find_colname(headers, &chrom_column_names);
    let col_pos = find_colname(headers, &pos_column_names);
    let col_ref = find_colname(headers, &ref_column_names);

    // println!("\n-------- TSV Col Names ---------");
    // println!("Chromosome: {}", col_chrom);
    // println!("Position: {}", col_pos);
    // println!("Reference: {}", col_ref);

    let valid_bases = ["A", "C", "T", "G"];

    // Initialise some counts
    let mut n_variants_supporting_0base = 0;
    let mut n_variants_supporting_1base = 0;
    let mut n_variants_with_chromosome_name_not_in_genome: i32 = 0;
    let mut n_variants_tested = 0;
    // Iterate Dataset
    for result in var_reader.deserialize() {
        // Parse Record
        let record: VarRecord = result.expect("Failed to parse tsv file record");

        let pos: usize = record[&col_pos]
            .parse()
            .expect("Failed to read position as numeric");

        let pos_if_1base = pos - 1;

        // Pull the ref base described in the variant file
        let ref_in_varfile = &record[&col_ref].to_uppercase();

        // Only consider variants with 1 base, that is one of valid_bases (A, C, G, or T) - No Ns allowed
        let n_refbases = ref_in_varfile.chars().count();
        if n_refbases != 1 || !valid_bases.contains(&ref_in_varfile.as_str()) {
            continue;
        }

        // Check if Chromosome Is in Ref genome 
        // (We explicitly do this check before querying fasta as otherwise htslib will just error 
        // when we try and fetch the sequence)
        if !valid_chromosomes.contains(&record[&col_chrom].to_string()) {
            n_variants_tested += 1;
            n_variants_with_chromosome_name_not_in_genome += 1;
            continue;
        }

        //Collect Sequence Information 
        let ref_base = fetch_seq(&reader, record[&col_chrom].to_string(), pos, pos, true);
        let ref_base_if_1base = fetch_seq(
            &reader,
            record[&col_chrom].to_string(),
            pos_if_1base,
            pos_if_1base,
            true,
        );

        // Check If Retrieved Base is the same as the base in variant file
        if ref_base == *ref_in_varfile {
            n_variants_supporting_0base += 1;
        }

        if ref_base_if_1base == *ref_in_varfile {
            n_variants_supporting_1base += 1;
        }

        n_variants_tested += 1;

        // println!("Reference Base for record {:?} is {}", record, ref_base)
    }

    let mut coord = Coordinate::Unsure;
    let mut ref_genome_matches = RefGenomeMatches::Unsure;

    if n_variants_supporting_0base != n_variants_supporting_1base {
        if n_variants_supporting_1base == n_variants_tested {
            coord = Coordinate::One;
            ref_genome_matches = RefGenomeMatches::Yes;
        }
        if n_variants_supporting_0base == n_variants_tested {
            coord = Coordinate::Zero;
            ref_genome_matches = RefGenomeMatches::Yes;
        }
    } else {
        let n_unexplained_variants = n_variants_tested
            - std::cmp::max(n_variants_supporting_1base, n_variants_supporting_0base);
        if n_unexplained_variants != 0 {
            ref_genome_matches = RefGenomeMatches::No
        }
    }

    let genome_filename = std::path::Path::new(path_genome)
        .file_name()
        .unwrap_or_else(|| panic!("Failed to parse genome as filepath"));

    // Print Results
    println!(
        "\n\n{EQUALS_HEADER}\nGenome: {}\n{EQUALS_HEADER}",
        genome_filename.to_owned().to_string_lossy().blue().bold()
    );

    println!("\n-------- Details --------");
    println!(
        "0-based matches to ref genome: {}/{}",
        n_variants_supporting_0base, n_variants_tested
    );
    println!(
        "1-based matches to ref genome: {}/{}",
        n_variants_supporting_1base, n_variants_tested
    );
    println!(
        "Chromosome names missing from ref genome: {}/{}",
        n_variants_with_chromosome_name_not_in_genome, n_variants_tested
    );
    
    println!("\n-------- Summary --------");
    println!("Coordinate System: {}", coord);
    println!("Reference Genome Matches: {}", ref_genome_matches);

    let result_summary = (genome_filename.to_str().unwrap_or_default().to_string(), ref_genome_matches, coord);
    result_summary
}

fn get_column_names<'a>(
    user_specified: Option<&'a str>,
    default_names: Vec<&'a str>,
) -> Vec<&'a str> {
    match user_specified {
        Some(name) => vec![name], // Now correctly associates the lifetime 'a with the input &str
        None => default_names, // And ensures that the default_names Vec<&'a str> matches the return type
    }
}

fn find_colname(headers: &StringRecord, valid_column_names: &[&str]) -> String {
    headers
        .into_iter()
        .find(|&col| valid_column_names.contains(&col))
        .unwrap_or_else(|| panic!("Couldn't find required column"))
        .to_string()
}

fn fetch_seq(
    reader: &rust_htslib::faidx::Reader,
    seqname: String,
    begin: usize,
    end: usize,
    upper: bool,
) -> String {
    let seq_string = reader
        .fetch_seq_string(seqname, begin, end)
        .expect("Failed to fetch sequencec string");

    // Convert to Uppercase
    if upper {
        return seq_string.to_uppercase();
    }

    // Return String
    seq_string
}


fn get_all_chromosomes(reader: &rust_htslib::faidx::Reader) -> Vec<String>{
    
    // Identify all chromosomes 
    let valid_chromosomes: Vec<String> = (0..reader.n_seqs()).map(|i| reader.seq_name(i as i32).unwrap()).collect();

    // Return
    valid_chromosomes
}