type alias Context = {
    matrix: Array[Array[Float64]], 
    b: Array[Float64], 
    v: Array[Float64], 
    c: Array[Float64],
    swap: Array[Int32]
}

def formatWith[a, e](f: a -> String & e, sep: String, a: Array[a]): String & Impure =
        let sb = StringBuilder.new();
        let step = (x,i) ->
            if (i == 0)
                StringBuilder.appendString!(sb, f(x))
            else {
                StringBuilder.appendString!(sb, sep);
                StringBuilder.appendString!(sb, f(x))
            };
        Array.foreachWithIndex(step, a);
        StringBuilder.toString(sb)

def arrayToString[a](toString: a -> String, array: Array[a]): String & Impure =
    let stringBuilder = StringBuilder.new();
    StringBuilder.appendString!(stringBuilder, "[");
    StringBuilder.appendString!(stringBuilder, formatWith((element) -> toString(element), ", ", array));
    StringBuilder.appendString!(stringBuilder, "]");
    StringBuilder.toString(stringBuilder)

def matrixToString(matrix: Array[Array[Float64]]): String & Impure =
    let stringBuilder = StringBuilder.new();
    StringBuilder.appendLine!(stringBuilder, "matrix: [");
    StringBuilder.appendString!(
        stringBuilder, 
        formatWith((array) -> arrayToString(Float64.toString, array), " -- ", matrix));
    StringBuilder.appendLine!(stringBuilder, "]");
    StringBuilder.toString(stringBuilder)

/// def matrixPrint(matrix: Array[Array[Float64]]): () & Impure =
///     Console.printLine(matrixToString(matrix))

def contextToString(context: Context): String & Impure =
    let stringBuilder = StringBuilder.new();
    StringBuilder.appendLine!(stringBuilder, "Context: ");
    StringBuilder.appendLine!(stringBuilder, matrixToString(context.matrix));
    StringBuilder.appendString!(stringBuilder, "b: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, context.b));
    StringBuilder.appendString!(stringBuilder, "v: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, context.v));
    StringBuilder.appendString!(stringBuilder, "c: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, context.c));
    StringBuilder.appendString!(stringBuilder, "swap: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Int32.toString, context.swap));
    StringBuilder.toString(stringBuilder)

def buildContext(matrixSize: Int32): Context & Impure = 
    {
        matrix = [[0.0; matrixSize]; matrixSize],
        b = [0.0; matrixSize],
        v = [0.0; matrixSize],
        c = [0.0; matrixSize],
        swap = [0; matrixSize]
    }  

def initMatrixhelper(context: Context, row: Array[Float64], i: Int32): () & Impure = 
    Array.foreachWithIndex((_, j) -> (context.matrix)[i][j] = Int32.toFloat64(if (j < i) 2 * (j + 1) else 2 * (i + 1)), row);
    (context.b)[i] = Int32.toFloat64(i);
    (context.swap)[i] = i

def initMatrix(context: Context): () & Impure = 
    Array.foreachWithIndex((row, i) -> initMatrixhelper(context, row, i), context.matrix)        


def swapDoubles(a: Ref[Float64], b: Ref[Float64]): () & Impure =
    let tmp = a;
    ref a := b;
    ref b := tmp;
    ()

def swapInts(a: Ref[Int32], b: Ref[Int32]): () & Impure =
    let tmp = a;
    ref a := b;
    ref b := tmp;
    ()

def iterateRows(context: Context, i: Int32, matrixSize: Int32, currow: Int32, big: Float64, irow: Int32): () & Impure =
    if (i == matrixSize) () 
    else {
        let tmp = (context.matrix)[i][currow];
        if (tmp != 0.0) {
            ref big := tmp;
            ref irow := i
        } else {
            iterateRows(context, i+1, matrixSize, currow, big, irow)
        }
    }
    
def swapRows(context: Context, i: Int32, matrixSize: Int32, currow: Int32, big: Float64, irow: Int32): () & Impure =
    if (i == matrixSize) ()
    else {
	    swapDoubles(ref (context.matrix)[irow][i], ref (context.matrix)[currow][i]);
        swapRows(context, i+1, matrixSize, currow, big, irow)
    }
    
	/// 		for(i = currow + 1; i < nsize; i++){
	/// 			matrix[currow][i] /= pivotVal;
	/// 		}
def divideElementsByPivot(context: Context, i: Int32, matrixSize: Int32, currow: Int32, pivotVal: Float64): () & Impure =
    if (i == matrixSize) ()
    else {
        ref (context.matrix)[currow][i] := pivotVal;
        divideElementsByPivot(context, i+1, matrixSize, currow, pivotVal)
    }

def getPivot(context: Context, currow: Int32, matrixSize: Int32): () & Impure =
	let big = (context.matrix)[currow][currow];
	let irow = currow;

	if (big == 0.0) {
        iterateRows(context, currow, matrixSize, currow, big, irow)
	} else ();			

	if (big == 0.0) {
        Console.printLine("The matrix is singular\n");
        import java.lang.Runtime:getRuntime();
        import java.lang.Runtime.exit(Int32);
        let rt = getRuntime();
        rt.exit(0)
	} else ();

	if (irow != currow){
        swapRows(context, currow, matrixSize, currow, big, irow);
        swapDoubles(ref (context.b)[irow], ref (context.b)[currow]);
        swapInts(ref (context.swap)[irow], ref (context.swap)[currow])
    }
    else ();

    let pivotVal = (context.matrix)[currow][currow];
    if (pivotVal != 1.0) {
        ref (context.matrix)[currow][currow] := 1.0;
        divideElementsByPivot(context, currow+1, matrixSize, currow, pivotVal);
		ref (context.b)[currow] := (context.b)[currow]/pivotVal
    } else ()


def main(): String & Impure =
    // The following variables would be command line arguments, but Flix doesn't
    // have the ability to recieve command line arguments
    let matrixSize = 16;
    /// let verify = true;
    /// let numThreads = 4;
    
    let context = buildContext(matrixSize);
    initMatrix(context);
    getPivot(context, 1, matrixSize);
	/// computeGauss(context, matrixSize);

    contextToString(context)

    