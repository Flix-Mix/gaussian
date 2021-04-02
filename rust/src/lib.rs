use std::process;

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

pub struct Data<'a> {
    pub nsize: usize,
    pub matrix: &'a mut Vec<Vec<f32>>,
    pub b: &'a mut Vec<f32>,
    pub c: &'a mut Vec<f32>,
    pub v: &'a mut Vec<f32>,
    pub swap: &'a mut Vec<u32>,
}

pub fn init(data: &mut Data) {
    for i in 0..data.nsize {
        data.matrix.push(Vec::with_capacity(data.nsize));
        for _ in 0..data.nsize {
            data.matrix[i].resize(data.nsize, 0.0);
        }
        data.b.push(i as f32);
        data.swap.push(i as u32);
        data.c.push(0.0);
        data.v.push(0.0);
    }

    for i in 0..data.matrix.len() {
        for j in 0..data.matrix.len() {
            let ii: f32 = i as f32;
            let jj: f32 = j as f32;
            data.matrix[i][j] = if jj < ii {
                2.0 * (jj + 1.0)
            } else {
                2.0 * (ii + 1.0)
            };
        }
    }
}

//#[allow(dead_code,unused_variables)]
pub fn compute_gauss(data: &mut Data) {
    for i in 0..data.nsize {
        pivot(data, i);

        let mut pivot_val; // = matrix[i][i];
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

//#[allow(dead_code,unused_variables)]
fn pivot(data: &mut Data, currow: usize) {
    let (mut irow, mut big, mut tmp);

    big = data.matrix[currow][currow];
    irow = currow;

    if big == 0.0 {
        //for i in currow..nsize {
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
        process::exit(-1);
    }

    if irow != currow {
        for i in currow..data.nsize {
            SWAP!(data.matrix[irow][i], data.matrix[currow][i]);
            // let tmp = matrix[irow][i]; //make this a macro
            // matrix[irow][i] = matrix[currow][i];
            // matrix[currow][i] = tmp;
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

#[allow(dead_code, unused_variables)]
fn solve_gauss(data: Data) {
    data.v[data.nsize - 1] = data.b[data.nsize - 1];
    for i in (0..data.nsize - 2).rev() {
        data.v[i] = data.b[i];
        for j in (i + 1..data.nsize - 1).rev() {
            data.v[i] -= data.matrix[i][j] * data.v[j];
        }
    }
    // for i in 0..nsize {
    //     c[i] = v[i];
    // }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn init_matrix_test() {
        let nsize = 3;
        let mut matrix: Vec<Vec<f32>> = Vec::with_capacity(nsize);
        let mut b: Vec<f32> = Vec::with_capacity(nsize);
        let mut c: Vec<f32> = Vec::with_capacity(nsize);
        let mut v: Vec<f32> = Vec::with_capacity(nsize);
        let mut swap: Vec<u32> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
        };
        init(&mut data);
        assert_eq!(matrix, [[2.0, 2.0, 2.0], [2.0, 4.0, 4.0], [2.0, 4.0, 6.0]]);
    }

    #[test]
    fn compute_gauss_test() {
        let nsize = 3;
        let mut matrix: Vec<Vec<f32>> = Vec::with_capacity(nsize);
        let mut b: Vec<f32> = Vec::with_capacity(nsize);
        let mut c: Vec<f32> = Vec::with_capacity(nsize);
        let mut v: Vec<f32> = Vec::with_capacity(nsize);
        let mut swap: Vec<u32> = Vec::with_capacity(nsize);
        let mut data = Data {
            nsize: nsize,
            matrix: &mut matrix,
            b: &mut b,
            c: &mut c,
            v: &mut v,
            swap: &mut swap,
        };
        init(&mut data);
        compute_gauss(&mut data);
        assert_eq!(matrix, [[1.0, 1.0, 1.0], [0.0, 1.0, 1.0], [0.0, 0.0, 1.0]]);
    }
}
