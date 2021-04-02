use std::env;
pub mod lib;

fn main() {
    //let args: Vec<String> = env::args().collect();

    let size: usize = env::args()
        .nth(1)
        .expect("no size given")
        .parse::<usize>()
        .unwrap();
    if size < 2 {
        panic!("size must be > 1");
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
    };

    lib::init(&mut data);
    lib::print(&data);
    lib::compute_gauss(&mut data);
    lib::print(&data);
    lib::solve_gauss(&mut data);
    lib::verify(&mut data);
}
