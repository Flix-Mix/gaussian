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
    Array.foreach(
        (row) -> StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, row)),
        matrix);
    StringBuilder.appendLine!(stringBuilder, "]");
    StringBuilder.toString(stringBuilder)

/// def matrixPrint(matrix: Array[Array[Float64]]): () & Impure =
///     Console.printLine(matrixToString(matrix))

def contextToString(context: Context): String & Impure =
    let stringBuilder = StringBuilder.new();
    StringBuilder.appendLine!(stringBuilder, "Context: ");
    StringBuilder.appendString!(stringBuilder, matrixToString(context.matrix));
    StringBuilder.appendString!(stringBuilder, "b: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, context.b));
    StringBuilder.appendString!(stringBuilder, "v: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, context.v));
    StringBuilder.appendString!(stringBuilder, "c: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, context.c));
    StringBuilder.appendString!(stringBuilder, "swap: ");
    StringBuilder.appendLine!(stringBuilder, arrayToString(Int32.toString, context.swap));
    StringBuilder.toString(stringBuilder)

def buildContextJLoop(context: Context, matrixSize: Int32, i: Int32, j: Int32): Context & Impure = 
    if (j == matrixSize) context
    else {
        if (j < i) {
            (context.matrix)[i][j] = Int32.toFloat64(2 * (j + 1));
            ()
        } else {
            ((context.matrix)[i])[j] = Int32.toFloat64(2 * (i + 1));
            ()
        };
        buildContextJLoop(context, matrixSize, i, j+1)
    }


def buildContextILoop(context: Context, matrixSize: Int32, i: Int32): Context & Impure = 
    if (i == matrixSize) context
    else {
        buildContextJLoop(context, matrixSize, i, 0);
        (context.b)[i] = Int32.toFloat64(i);
        (context.swap)[i] = i;
        buildContextILoop(context, matrixSize, i+1)
    }


def buildMatrix(matrix: Array[Array[Float64]], i: Int32, matrixSize: Int32): Array[Array[Float64]] & Impure = 
    if (i == matrixSize) matrix
    else {
        matrix[i] = [0.0; matrixSize];
        buildMatrix(matrix, i+1, matrixSize)
    } 

def buildContext(matrixSize: Int32): Context & Impure = 
    let context = {
        matrix = [[0.0; matrixSize]; matrixSize],
        b = [0.0; matrixSize],
        v = [0.0; matrixSize],
        c = [0.0; matrixSize],
        swap = [0; matrixSize]
    };
    ref (context.matrix) := buildMatrix(context.matrix, 0, matrixSize);
    ref context := buildContextILoop(context, matrixSize, 0);
    context


def computeGaussKLoop(context: Context, pivotValue: Float64, i: Int32, j: Int32, k: Int32): () & Impure =
    if (k < Array.length((context.matrix)[i])) {
        (context.matrix)[j][k] = (pivotValue * (context.matrix)[i][k]) - (context.matrix)[j][k];
        computeGaussKLoop(context, pivotValue, i, j, k + 1)
    } else ()

def computeGaussJLoop(context: Context, i: Int32, j: Int32): () & Impure =
    if (j < Array.length((context.matrix)[i])) {
        let pivotValue = (context.matrix)[j][i];
        (context.matrix)[j][i] = 0.0;
        computeGaussKLoop(context, pivotValue, i, j, i + 1);
        (context.b)[j] = (pivotValue * (context.b)[i]) - (context.b)[j];
        computeGaussJLoop(context, i, j + 1)
    } else ()

def computeGauss(context: Context, matrixSize: Int32): () & Impure =
    let computeGaussILoop = (row: Array[Float64], i: Int32) -> {
        getPivot(context, i, matrixSize);    
        Array.foreachWithIndex((_, _) -> computeGaussJLoop(context, i, i + 1), row)
    };
    Array.foreachWithIndex((row, i) -> computeGaussILoop(row, i), context.matrix)    


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
        // Swap matrix[irow] and matrix[currow]
        let tmp = (context.matrix)[irow][i];
        (context.matrix)[irow][i] = (context.matrix)[currow][i];
        (context.matrix)[currow][i] = tmp;

        swapRows(context, i+1, matrixSize, currow, big, irow)
    }
    
def divideElementsByPivot(context: Context, i: Int32, matrixSize: Int32, currow: Int32, pivotVal: Float64): () & Impure =
    if (i == matrixSize) ()
    else {
        (context.matrix)[currow][i] = (context.matrix)[currow][i]/pivotVal;
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

        // Swap b[irow] and b[currow]
        let tmp = (context.b)[irow];
        (context.b)[irow] = (context.b)[currow];
        (context.b)[currow] = tmp;

        // Swap swap[irow] and swap[currow]
        let tmpInt = (context.swap)[irow];
        (context.swap)[irow] = (context.swap)[currow];
        (context.swap)[currow] = tmpInt
    }
    else ();

    let pivotVal = (context.matrix)[currow][currow];
    if (pivotVal != 1.0) {
        (context.matrix)[currow][currow] = 1.0;
        divideElementsByPivot(context, currow+1, matrixSize, currow, pivotVal);
		(context.b)[currow] = (context.b)[currow]/pivotVal
    } else ()

def solveGaussVJLoop(context: Context, i: Int32, j: Int32, matrixSize: Int32): Context & Impure =
    if (j <= i) context
    else {
	    (context.v)[i] = (context.v)[i] - ((context.matrix)[i][j] * (context.v)[j]);
        solveGaussVJLoop(context, i, j-1, matrixSize)
    }

def solveGaussVILoop(context: Context, i: Int32, matrixSize: Int32): Context & Impure =
    if (i < 0) context
    else {
        (context.v)[i] = (context.b)[i];
        solveGaussVJLoop(context, i, matrixSize-1, matrixSize);
        solveGaussVILoop(context, i-1, matrixSize)
    }
    
def solveGaussCILoop(context: Context, i: Int32, matrixSize: Int32): Context & Impure =
    if (i >= matrixSize) context
    else {
        (context.c)[i] = (context.v)[i];
        solveGaussCILoop(context, i+1, matrixSize)
    }

def solveGauss(context: Context, matrixSize: Int32): Context & Impure =
	(context.v)[matrixSize-1] = (context.b)[matrixSize-1];
    ref context := solveGaussVILoop(context, matrixSize-2, matrixSize);
    solveGaussCILoop(context, 0, matrixSize)


def verifyLoop(context: Context, i: Int32, matrixSize: Int32): Unit & Impure =
    if (i == matrixSize) ()
    else {
        let stringBuilder = StringBuilder.new();
        let bVal = Float64.toString((context.b)[i]);
        let cVal = Float64.toString((context.c)[i]);
        StringBuilder.appendString!(stringBuilder, bVal);
        StringBuilder.appendString!(stringBuilder, " ");
        StringBuilder.appendString!(stringBuilder, cVal);
        let row = StringBuilder.toString(stringBuilder);
        Console.printLine(row);
        verifyLoop(context, i+1, matrixSize)
    }

def main(): String & Impure =
    // The following variables would be command line arguments, but Flix doesn't
    // have the ability to recieve command line arguments
    let matrixSize = 128;
    let verify = true;
    let printContext = false;
    /// let numThreads = 4;
    let context = buildContext(matrixSize);

    computeGauss(context, matrixSize);

	if (verify) ref context := solveGauss(context, matrixSize)
	else ();

    if (printContext) Console.printLine(contextToString(context))
	else ();

	if (verify) {
		verifyLoop(context, 0, matrixSize)
	} else ();
    
    "end"