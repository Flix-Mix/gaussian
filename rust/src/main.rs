pub mod lib;
extern crate clap;
use clap::{App, Arg};
use std::sync::{Arc};// Mutex};
use std::time::Instant;

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
                .default_value("0")
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

    let matrix: Vec<Vec<f64>> = Vec::with_capacity(size);
    let b: Vec<f64> = Vec::with_capacity(size);
    let c: Vec<f64> = Vec::with_capacity(size);
    let v: Vec<f64> = Vec::with_capacity(size);
    let swap: Vec<u64> = Vec::with_capacity(size);

    let data = lib::Data {
        nsize: size,
        matrix,
        b,
        c,
        v,
        swap,
        num_threads: num_of_threads,
    };

    // lib::init(&mut data);
    // if verbose {
    //     lib::print(&data);
    // }
    // // if num_of_threads > 0 {
    // //     data = lib::compute_gauss_p(data);
    // // } else {
    // //     lib::compute_gauss(&mut data)
    // // }
    // lib::compute_gauss(&mut data);
    // if verbose {
    //     lib::print(&data);
    // }
    let now = Instant::now();
    let data_arc = lib::initp(data);
    let time = now.elapsed().as_nanos();
    println!("Program finished in {} sec", (time as f64)/10e8);
    if verbose {
        let guard = Arc::try_unwrap(data_arc).unwrap();
        let inner = guard.lock().unwrap(); 
        let t  = &*inner;
        lib::print(t);
    }
    // let guard = Arc::try_unwrap(data_arc).unwrap();
    // let inner = guard.lock().unwrap(); 
    // let mut t  = &*inner;
    // lib::solve_gauss(t);
    // lib::verify(&data);
    // ib::compute_gauss(&mut data);
    // if verbose {
    //     lib::print(&data);
    // }
    // lib::verify(&data);
}
