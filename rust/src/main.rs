pub mod lib;
extern crate clap;
use clap::{App, Arg};

fn main() {
    let matches = App::new("Gaussian internal")
        .author("Aidan Goldfarb <agoldfa7@u.rochester.edu>")
        .about("Gausian elimination seq and parallel using channels")
        .arg(
            Arg::with_name("SIZE")
                .short("s")
                .long("size")
                .help("Sets matrix dimensions")
                .required(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("VERBOSITY")
                .short("v")
                .long("verbose")
                .help("Sets the verbosity")
                .takes_value(false)
                .required(false),
        )
        .arg(
            Arg::with_name("NUM_THREADS")
                .short("n")
                .long("num_threads")
                .help("Sets the number of threads to spawn")
                .takes_value(true)
                .default_value("1")
                .required(false),
        )
        .get_matches();

    let size: usize = matches.value_of("SIZE").unwrap().parse::<usize>().unwrap();
    let mut verbose = false;
    let num_of_threads = matches
        .value_of("NUM_THREADS")
        .unwrap()
        .parse::<usize>()
        .unwrap();

    //println!("size: {}, verbose: {}, num threads {}", size, verbose, num_threads);

    if matches.is_present("VERBOSITY") {
        verbose = true;
    }

    let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(size);
    let mut b: Vec<f64> = Vec::with_capacity(size);
    let mut c: Vec<f64> = Vec::with_capacity(size);
    let mut v: Vec<f64> = Vec::with_capacity(size);
    let mut swap: Vec<u64> = Vec::with_capacity(size);

    let mut data = lib::Data {
        nsize: size,
        matrix: &mut matrix,
        b: &mut b,
        c: &mut c,
        v: &mut v,
        swap: &mut swap,
        num_threads: num_of_threads,
    };

    lib::init(&mut data);
    if verbose {
        lib::print(&data);
    }
    // lib::compute_gauss(&mut data);
    // if verbose {
    //     lib::print(&data);
    // }
    // lib::solve_gauss(&mut data);
    // lib::verify(&data);
    lib::compute_gauss_p(&mut data);
    if verbose {
        lib::print(&data);
    }
    lib::verify(&data);
}
