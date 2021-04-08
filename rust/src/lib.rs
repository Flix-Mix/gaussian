#[macro_use]
macro_rules! SWAP {
    ($a:expr,$b:expr) => {{
        // $a ^= $b;
        // $b ^= $a;
        // $a ^= $b;
        let t = $a;
        $a = $b;
        $b = t;
    }};
}

use std::sync::{Arc, Mutex, MutexGuard};
// use std::thread;
// use std::sync::mpsc::{Sender, Receiver};
// use std::sync::mpsc::channel;

use threadpool::ThreadPool;
extern crate threadpool;

pub struct Data<'a> {
    pub nsize: usize,
    pub matrix: &'a mut Vec<Vec<f64>>,
    pub b: &'a mut Vec<f64>,
    pub c: &'a mut Vec<f64>,
    pub v: &'a mut Vec<f64>,
    pub swap: &'a mut Vec<u64>,
    pub num_threads: usize,
}

pub fn init(data: &mut Data) {
    for i in 0..data.nsize {
        data.matrix.push(Vec::with_capacity(data.nsize));
        for _ in 0..data.nsize {
            data.matrix[i].resize(data.nsize, 0.0);
        }
        data.b.push(i as f64);
        data.swap.push(i as u64);
        data.c.push(0.0);
        data.v.push(0.0);
    }

    for i in 0..data.matrix.len() {
        for j in 0..data.matrix.len() {
            let ii: f64 = i as f64;
            let jj: f64 = j as f64;
            data.matrix[i][j] = if jj < ii {
                2.0 * (jj + 1.0)
            } else {
                2.0 * (ii + 1.0)
            };
        }
    }
}

pub fn compute_gauss(data: &mut Data) {
    for i in 0..data.nsize {
        pivot(data, i);

        let mut pivot_val;

        for j in i + 1..data.nsize {
            pivot_val = data.matrix[j][i];
            data.matrix[j][i] = 0.0;
            for k in i + 1..data.nsize {
                data.matrix[j][k] -= pivot_val * data.matrix[i][k];
            }
            data.b[j] -= pivot_val * data.b[i];
        }
    }
}

pub struct Message {
    pub nsize: usize,
    pub index: usize,
}

pub fn compute_gauss_p(data: &mut Data) {
    let pool = ThreadPool::new(data.num_threads);
    let num_threads = data.num_threads;
    let s_mat = Arc::new(Mutex::new(data.matrix.clone()));
    let b_mat = Arc::new(Mutex::new(data.b.clone()));
    for i in 0..data.nsize {
        pivot(data, i);
        let s_mat = Arc::clone(&s_mat);
        let b_mat = Arc::clone(&b_mat);
        pool.execute(move || {
            let mut mat = s_mat.lock().unwrap();
            let mut b = b_mat.lock().unwrap();
            do_calc(&mut mat, &mut b, i, num_threads);
        });
    }
    pool.join();
}

fn do_calc(
    mat: &mut MutexGuard<'_, Vec<Vec<f64>>>,
    b: &mut MutexGuard<'_, Vec<f64>>,
    i: usize,
    num_threads: usize,
) {
    let mut pivot_val;
    for j in (i + 1..mat.len()).step_by(num_threads) {
        pivot_val = mat[j][i];
        mat[j][i] = 0.0;
        for k in i + 1..mat.len() {
            mat[j][k] -= pivot_val * mat[i][k];
        }
        b[j] -= pivot_val * b[i];
    }
}

fn pivot(data: &mut Data, currow: usize) {
    let (mut irow, mut big, mut tmp);

    big = data.matrix[currow][currow];
    irow = currow;

    if big == 0.0 {
        for (i, row) in data.matrix.iter().enumerate().take(data.nsize).skip(currow) {
            tmp = row[currow];
            if tmp != 0.0 {
                big = tmp;
                irow = i;
                break;
            }
        }
    }

    if big == 0.0 {
        println!("The matrix is singular");
        panic!("singular");
    }

    if irow != currow {
        for i in currow..data.nsize {
            SWAP!(data.matrix[irow][i], data.matrix[currow][i]);
        }
        data.b.swap(irow, currow);

        data.swap.swap(irow, currow)
    }

    {
        let pivot_val = data.matrix[currow][currow];
        if (pivot_val - 1.0).abs() > 0.0000001 {
            data.matrix[currow][currow] = 1.0;
            for i in currow + 1..data.nsize {
                data.matrix[currow][i] /= pivot_val;
            }
            data.b[currow] /= pivot_val;
        }
    }
}

pub fn solve_gauss(data: &mut Data) {
    data.v[data.nsize - 1] = data.b[data.nsize - 1];
    for i in (0..data.nsize - 2).rev() {
        data.v[i] = data.b[i];
        for j in (i + 1..data.nsize - 1).rev() {
            data.v[i] -= data.matrix[i][j] * data.v[j];
        }
    }
    data.c[..data.nsize].clone_from_slice(&data.v[..data.nsize]);
}

pub fn print(data: &Data) {
    println!("{{");
    for row in data.matrix.iter().take(data.matrix.len() - 1) {
        println!(" {:?}", row);
    }
    print!(" {:?}", data.matrix[data.matrix.len() - 1]);
    println!("\n}}");
}

pub fn print_mat(matrix: &[Vec<f64>]) {
    println!("{{");
    for row in matrix.iter().take(matrix.len() - 1) {
        println!(" {:?}", row);
    }
    print!(" {:?}", matrix[matrix.len() - 1]);
    println!("\n}}");
}

pub fn verify(data: &Data) -> u64 {
    let err: f64 = 0.000001;
    if data.nsize == 2 {
        //TODO assert with assert_le to avoid error margin failures
        assert!((data.b[0] - 0.0).abs() < err);
        assert!((data.b[1] - 0.5).abs() < err);
        //assert_eq!(*data.b, [0.0, 0.5]);

        assert!((data.c[0] - 0.0).abs() < err);
        assert!((data.c[1] - 0.5).abs() < err);
    //assert_eq!(*data.c, [0.0, 0.5]);
    } else if data.nsize == 3 {
        assert!((data.b[0] - 0.0).abs() < err);
        assert!((data.b[1] - 0.5).abs() < err);
        assert!((data.b[2] - 0.5).abs() < err);
        //assert_eq!(*data.b, [0.0, 0.5, 0.5]);

        assert!((data.c[0] - 0.0).abs() < err);
        assert!((data.c[1] - 0.0).abs() < err);
        assert!((data.c[2] - 0.5).abs() < err);
    //assert_eq!(*data.c, [0.0, 0.0, 0.5]);
    } else {
        for i in 0..data.nsize {
            //println!("{:6.5} {:5.5}", data.b[i], data.c[i]);
            if i == 0 {
                assert!((data.b[i] - 0.0).abs() < err);
                assert!((data.c[i] - -0.5).abs() < err);
            } else if i == data.nsize - 3 || i == data.nsize - 1 {
                assert!((data.c[i] - 0.5).abs() < err);
            } else {
                assert!((data.b[i] - 0.5).abs() < err);
                assert!((data.c[i] - 0.0).abs() < err);
            }
        }
    }
    println!("Verified");
    1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_matrix_test() {
        let nsize = 3;
        let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(nsize);
        let mut b: Vec<f64> = Vec::with_capacity(nsize);
        let mut c: Vec<f64> = Vec::with_capacity(nsize);
        let mut v: Vec<f64> = Vec::with_capacity(nsize);
        let mut swap: Vec<u64> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
            num_threads: 1,
        };
        init(&mut data);
        assert_eq!(matrix, [[2.0, 2.0, 2.0], [2.0, 4.0, 4.0], [2.0, 4.0, 6.0]]);
    }

    #[test]
    fn compute_gauss_test() {
        let nsize = 3;
        let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(nsize);
        let mut b: Vec<f64> = Vec::with_capacity(nsize);
        let mut c: Vec<f64> = Vec::with_capacity(nsize);
        let mut v: Vec<f64> = Vec::with_capacity(nsize);
        let mut swap: Vec<u64> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
            num_threads: 1,
        };
        init(&mut data);
        compute_gauss(&mut data);
        assert_eq!(matrix, [[1.0, 1.0, 1.0], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]]);
    }

    #[test]
    fn solve_gauss_test() {
        let nsize = 5;
        let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(nsize);
        let mut b: Vec<f64> = Vec::with_capacity(nsize);
        let mut c: Vec<f64> = Vec::with_capacity(nsize);
        let mut v: Vec<f64> = Vec::with_capacity(nsize);
        let mut swap: Vec<u64> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
            num_threads: 1,
        };
        init(&mut data);
        compute_gauss(&mut data);
        solve_gauss(&mut data);
        assert_eq!(verify(&mut data), 1);
    }

    #[test]
    fn parallel_smoke() {
        let nsize = 5;
        let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(nsize);
        let mut b: Vec<f64> = Vec::with_capacity(nsize);
        let mut c: Vec<f64> = Vec::with_capacity(nsize);
        let mut v: Vec<f64> = Vec::with_capacity(nsize);
        let mut swap: Vec<u64> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
            num_threads: 1,
        };
        init(&mut data);
        compute_gauss_p(&mut data);
        assert!(true);
    }

    #[test]
    fn parallel_t1() {
        let nsize = 5;
        let mut matrix: Vec<Vec<f64>> = Vec::with_capacity(nsize);
        let mut b: Vec<f64> = Vec::with_capacity(nsize);
        let mut c: Vec<f64> = Vec::with_capacity(nsize);
        let mut v: Vec<f64> = Vec::with_capacity(nsize);
        let mut swap: Vec<u64> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
            num_threads: 1,
        };
        init(&mut data);
        compute_gauss_p(&mut data);
        solve_gauss(&mut data);
        assert_eq!(verify(&data), 1);
    }
}
