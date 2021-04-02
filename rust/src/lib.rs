use std::process;

pub fn init(nsize: usize, matrix: &mut Vec<Vec<f32>>, b: &mut Vec<f32>, swap: &mut Vec<u32>, c: &mut Vec<f32>, v: &mut Vec<f32>){
	for i in 0..nsize{
		matrix.push(Vec::new());
		for _ in 0..nsize{
			matrix[i].push(0.0);
		}
		b.push(i as f32);
		swap.push(i as u32);
		c.push(0.0);
		v.push(0.0);
	}

	for i in 0..matrix.len(){
		for j in 0..matrix.len(){
			let ii: f32 = i as f32;
			let jj: f32 = j as f32;
			matrix[i][j] = if jj < ii {2.0*(jj+1.0)} else {2.0*(ii+1.0)};
		}
	}
}

//#[allow(dead_code,unused_variables)]
pub fn compute_gauss(nsize: usize, matrix: &mut Vec<Vec<f32>>, b: &mut Vec<f32>, swap: &mut Vec<u32>){
	for i in 0..nsize{
		pivot(matrix,nsize,i,b,swap);

		let mut pivot_val;// = matrix[i][i];
		for j in i+1..nsize{
			pivot_val = matrix[j][i];
			matrix[j][i] = 0.0;
			for k in i+1..nsize{
				matrix[j][k] -= pivot_val * matrix[i][k];
			}
			b[j] -= pivot_val * b[i];
		}
	}
}

//#[allow(dead_code,unused_variables)]
fn pivot(matrix: &mut Vec<Vec<f32>>, nsize: usize, currow: usize, b: &mut Vec<f32>, swap: &mut Vec<u32>){
	let (mut irow, mut big, mut tmp);

	big = matrix[currow][currow];
	irow = currow;

	if big == 0.0 {
		for i in currow..nsize{
			tmp = matrix[i][currow];
			if tmp != 0.0 {
				big = tmp;
				irow = i;
				break;
			}
		}
	}

	if big == 0.0 {
		println!("The matrix is singular");
		process::exit(-1);
	}

	if irow != currow {
		for i in currow..nsize {
			let tmp = matrix[irow][i]; //make this a macro
			matrix[irow][i] = matrix[currow][i];
			matrix[currow][i] = tmp;
		}
		let tmp = b[irow];
		b[irow] = b[currow];
		b[currow] = tmp;

		let tmp = swap[irow];
		swap[irow] = swap[currow];
		swap[currow] = tmp;
	}

	{
		let pivot_val = matrix[currow][currow];
		if pivot_val != 1.0{
			matrix[currow][currow] = 1.0;
			for i in currow+1..nsize{
				matrix[currow][i] /= pivot_val;
			}
			b[currow] /= pivot_val;
		}
	}
}

#[allow(dead_code,unused_variables)]
fn solve_gauss(matrix: &mut Vec<Vec<f32>>, nsize: usize, b: &Vec<f32>, c: &mut Vec<f32>, v: &mut Vec<f32>){
	v[nsize-1] = b[nsize-1];
	for i in (0..nsize-2).rev(){
		v[i] = b[i];
		for j in (i+1..nsize-1).rev(){
			v[i] -= matrix[i][j] * v[j];
		}
	}
	for i in 0..nsize{
		c[i] = v[i];
	}
}


pub fn print(matrix: &Vec<Vec<f32>>){
	print!("{{\n");
	for i in 0..matrix.len()-1{
		println!(" {:?}",matrix[i]);
	}
	print!(" {:?}",matrix[matrix.len()-1]);
	println!("\n}}");
}


// #[macro_export]
// macro_rules! SWAP {
// 	($a:expr,$b:expr) => {
// 		{
// 			let t = $a;
// 			$a = $b;
// 			$b = $t;
// 		}
// 	}
// }

#[cfg(test)]
mod tests {	
	use super::*;

    #[test]
    fn init_matrix_test() {
    	let nsize = 3;
        let mut matrix: Vec<Vec<f32>> = Vec::new();
        let mut b: Vec<f32> = Vec::new();
        let mut c: Vec<f32> = Vec::new();
		let mut v: Vec<f32> = Vec::new();
		let mut swap: Vec<u32> = Vec::new();
        init(nsize,&mut matrix,&mut b,&mut swap, &mut c, &mut v);
        assert_eq!(matrix, [[2.0,2.0,2.0],[2.0,4.0,4.0],[2.0,4.0,6.0]]);
    }

    #[test]
    fn compute_gauss_test() {
    	let nsize = 3;
        let mut matrix: Vec<Vec<f32>> = Vec::new();
        let mut b: Vec<f32> = Vec::new();
        let mut c: Vec<f32> = Vec::new();
		let mut v: Vec<f32> = Vec::new();
		let mut swap: Vec<u32> = Vec::new();
        init(nsize,&mut matrix,&mut b,&mut swap, &mut c, &mut v);
        compute_gauss(nsize,&mut matrix,&mut b,&mut swap);
        assert_eq!(matrix, [[1.0,1.0,1.0],[0.0,1.0,1.0],[0.0,0.0,1.0]]);
    }
} 
 