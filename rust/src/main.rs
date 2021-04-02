use std::env;
pub mod lib;

fn main() {
    let args: Vec<String> = env::args().collect();

    let size: usize = args[1].parse::<usize>().unwrap();
    let mut matrix: Vec<Vec<f32>> = Vec::with_capacity(size);
    let mut b: Vec<f32> = Vec::with_capacity(size);
    let mut c: Vec<f32> = Vec::with_capacity(size);
    let mut v: Vec<f32> = Vec::with_capacity(size);
    let mut swap: Vec<u32> = Vec::with_capacity(size);

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
}
