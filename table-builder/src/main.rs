use fixed::types::U0F64;
use std::fs::File;
use std::io::Write;
use std::{env, process};

fn build_cordic_atan_table(num_entries: u8) -> Vec<u64> {
    let mut table = Vec::with_capacity(num_entries as usize);
    let mut angle = 1.0f64;

    for _ in 0..num_entries {
        let atan = angle.atan();
        let ival = U0F64::from_num(atan);
        table.push(ival.to_bits());
        angle = angle * 0.5;
    }

    table
}

// NOTE: we generate exp(x) - 1.0 instead of exp(x) directly
// so that the results can fit into a U0F64 like all our other
// CORDIC lookup tables.
fn build_cordic_exp_minus_one_table(num_entries: u8) -> Vec<u64> {
    let mut table = Vec::with_capacity(num_entries as usize);
    let mut angle = 0.5f64;

    for _ in 0..num_entries {
        let exp = angle.exp();
        let ival = U0F64::from_num(exp - 1.0);
        table.push(ival.to_bits());
        angle = angle * 0.5;
    }

    table
}

fn vec_u64_to_le_bytes(data: Vec<u64>) -> Vec<u8> {
    let mut bytes = vec![];
    data.into_iter()
        .for_each(|d| bytes.extend(&d.to_le_bytes()));
    bytes
}

fn write_table(table_filepath: &str, num_bits: Option<u8>, data: Vec<u8>) {
    let mut table_file = File::create(table_filepath).expect("unable to open file");
    if let Some(num_bits) = num_bits {
        table_file
            .write_all(&[num_bits])
            .expect("unable to write to file");
    }
    table_file
        .write_all(&data)
        .expect("unable to write to file");
}

fn print_usage(error: Option<&str>) {
    if let Some(error) = error {
        eprintln!("Error: {}", error);
    }
    println!(
        "usage: table-builder {{cordic-atan|cordic-exp-minus-one|cordic-cumprod}} num_bits table_filepath"
    );
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() == 1 {
        print_usage(None);
        process::exit(0);
    }
    if args.len() != 4 {
        print_usage(Some("Incorrect number of parameters"));
        process::exit(1);
    }
    let num = args[2].parse::<u8>().expect("unable to parse num_bits");
    match args[1].as_str() {
        "cordic-atan" => {
            let table = build_cordic_atan_table(num);
            let data = vec_u64_to_le_bytes(table);
            write_table(args[3].as_str(), None, data);
        }
        "cordic-exp-minus-one" => {
            let table = build_cordic_exp_minus_one_table(num);
            let data = vec_u64_to_le_bytes(table);
            write_table(args[3].as_str(), None, data);
        }
        _ => {
            print_usage(Some("Invalid table type"));
            process::exit(1)
        }
    }
}
