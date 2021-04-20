type alias Context = {
    matrix: Array[Array[Float64]], 
    b: Array[Float64], 
    v: Array[Float64], 
    c: Array[Float64],
    swap: Array[Int32]
}

def buildContext(matrixSize: Int32): Context & Impure = 
    {
        matrix = [[0.0; matrixSize]; matrixSize],
        b = [0.0; matrixSize],
        v = [0.0; matrixSize],
        c = [0.0; matrixSize],
        swap = [0; matrixSize]
    }  

def initMatrixhelper(context: Context, array: Array[Float64], i: Int32): () & Impure = 
    Array.foreachWithIndex((element, j) -> ref element := Int32.toFloat64(if (j < i) 2 * (j + 1) else 2 * (i + 1)), array);
    Array.put(context.b, i, Int32.toFloat64(i));
    let _ = Array.put(context.swap, i, i);

def initMatrix(context: Context): () & Impure = 
    Array.foreachWithIndex((array, i) -> initMatrixhelper(context, array, i), context.matrix)        

def main(): String & Impure =
    // The following variables would be command line arguments, but Flix doesn't
    // have the ability to recieve command line arguments
    let matrixSize = 1024;
    let verify = true;
    let numThreads = 4;
    
    let context = buildContext(matrixSize);
    initMatrix(context);

    "Hey Flix!"

    