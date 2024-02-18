use std::string::String;

use csv::StringRecord;

// struct VarRecord {
//     chrom: String,
//     pos: usize,
//     refbase: String,
// }

//type VarRecord = (String, usize, String);
type VarRecord = std::collections::HashMap<String, String>;

// #[derive(Debug, Deserialize)]
// #[serde(rename_all = "PascalCase")]
// struct VarRecord {
//     chrom: String,
//     pos: usize,
//     refbase: String,
// }

fn main() {
    // Setup Paths
    let path_genome = "inst/genome11.fasta";
    let path_tsv = "inst/variants_1based";

    // Create Reader
    let reader_hg38 = rust_htslib::faidx::Reader::from_path(path_genome)
        .expect("Failed to read hg38 genome file");

    // Read TSV file
    // Create Reader
    let mut var_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path_tsv)
        .expect("Failed to read TSV file");

    //Setup Valid Chrom column Names
    let mut chrom_column_names = vec!["chrom", "Chrom", "CHROM", "chromosome", "Chromosome"];
    let mut pos_column_names = vec!["pos", "Pos", "POS", "position", "Position"];
    let mut ref_column_names = vec!["ref", "REF", "Ref", "Reference", "Reference_Allele"];

    // Allow users to add column names (Will need to take from user input)
    // Note by default values will be None and the default valid names above will be used
    let user_specified_chrom = None;
    let user_specified_pos = None;
    let user_specified_ref = None;

    if user_specified_chrom.is_some() {
        chrom_column_names = vec![user_specified_chrom.unwrap()];
    }
    if user_specified_pos.is_some() {
        pos_column_names = vec![user_specified_pos.unwrap()];
    }
    if user_specified_ref.is_some() {
        ref_column_names = vec![user_specified_ref.unwrap()];
    }

    // Ensure File Has Headers
    if !var_reader.has_headers() {
        panic!("genomeguesser only works if the variant file has a header row. Please add a header and try again")
    }

    // Figure Out Column Names
    let col_chrom;
    let col_pos;
    let col_ref;
    let headers = var_reader.headers().expect("Failed to find headers");
    col_chrom = find_colname(headers, &chrom_column_names);
    col_pos = find_colname(headers, &pos_column_names);
    col_ref = find_colname(headers, &ref_column_names);

    println!("Chromosome Column Name: {}", col_chrom);
    println!("Position Column Name: {}", col_pos);
    println!("Reference Column Name: {}", col_ref);

    let valid_bases = ["A", "C", "T", "G"];

    let mut n_variants_supporting_0base = 0;
    let mut n_variants_supporting_1base = 0;
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

        //Collect Sequence Information ()
        let ref_base = fetch_seq(&reader_hg38, record[&col_chrom].to_string(), pos, pos, true);
        let ref_base_if_1base = fetch_seq(
            &reader_hg38,
            record[&col_chrom].to_string(),
            pos_if_1base,
            pos_if_1base,
            true,
        );

        // Check If Retrieved Base is the same as the base in variant file
        if ref_base == ref_in_varfile.to_string() {
            n_variants_supporting_0base += 1;
        }

        if ref_base_if_1base == ref_in_varfile.to_string() {
            n_variants_supporting_1base += 1;
        }

        n_variants_tested += 1;

        // println!("Reference Base for record {:?} is {}", record, ref_base)
    }

    // Print Results
    println!(
        "\nSupport for Genomes
        0-based matches to ref genome: {}/{}, 
        1-based matches to ref genome: {}/{}",
        n_variants_supporting_0base,
        n_variants_tested,
        n_variants_supporting_1base,
        n_variants_tested
    )
}

fn find_colname(headers: &StringRecord, valid_column_names: &Vec<&str>) -> String {
    let colname = headers
        .into_iter()
        .find(|&col| valid_column_names.contains(&col))
        .expect("Couldnt Find Chrom Col");

    colname.to_string()
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
