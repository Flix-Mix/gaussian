use std::env;
pub mod lib;

fn main(){
	let args: Vec<String> = env::args().collect();
	
	let nsize: usize = args[1].parse::<usize>().unwrap();
	let mut matrix: Vec<Vec<f32>> = Vec::new();
	let mut b: Vec<f32> = Vec::new();
	let mut c: Vec<f32> = Vec::new();
	let mut v: Vec<f32> = Vec::new();
	let mut swap: Vec<u32> = Vec::new();


	lib::init(nsize,&mut matrix,&mut b,&mut swap, &mut c, &mut v);
	lib::print(&matrix);
	lib::compute_gauss(nsize,&mut matrix,&mut b,&mut swap);
	lib::print(&matrix);
}

