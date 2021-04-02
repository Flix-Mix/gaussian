use std::env;
pub mod lib;

fn main(){
	let args: Vec<String> = env::args().collect();
	
	let nsize: usize = args[1].parse::<usize>().unwrap();
	
	let mut matrix: Vec<Vec<f32>> = Vec::with_capacity(nsize);
	let mut b: Vec<f32> = Vec::with_capacity(nsize);
	let mut c: Vec<f32> = Vec::with_capacity(nsize);
	let mut v: Vec<f32> = Vec::with_capacity(nsize);
	let mut swap: Vec<u32> = Vec::with_capacity(nsize);


	lib::init(nsize,&mut matrix,&mut b,&mut swap, &mut c, &mut v);
	lib::print(&matrix);
	lib::compute_gauss(nsize,&mut matrix,&mut b,&mut swap);
	lib::print(&matrix);
}