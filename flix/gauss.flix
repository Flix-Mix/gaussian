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
        (array) -> StringBuilder.appendLine!(stringBuilder, arrayToString(Float64.toString, array)),
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

def buildContext(matrixSize: Int32): Context & Impure = 
    let context = {
        matrix = [[0.0; matrixSize]; matrixSize],
        b = [0.0; matrixSize],
        v = [0.0; matrixSize],
        c = [0.0; matrixSize],
        swap = [0; matrixSize]
    };

    {
        let initMatrixHelper = (ctx: Context, row: Array[Float64], i: Int32) -> {
            Array.foreachWithIndex((_, j) -> (ctx.matrix)[i][j] = Int32.toFloat64(if (j < i) 2 * (j + 1) else 2 * (i + 1)), row);
            (ctx.b)[i] = Int32.toFloat64(i);
            (ctx.swap)[i] = i
        };
        Array.foreachWithIndex((row, i) -> initMatrixHelper(context, row, i), context.matrix)
    };

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

def computeGauss(context: Context): () & Impure =
    let computeGaussILoop = (row: Array[Float64], i: Int32) -> {
        // getPivot
        Array.foreachWithIndex((_, _) -> computeGaussJLoop(context, i, i + 1), row)
    };
    Array.foreachWithIndex((row, i) -> computeGaussILoop(row, i), context.matrix)    

def main(): String & Impure =
    // The following variables would be command line arguments, but Flix doesn't
    // have the ability to recieve command line arguments
    let matrixSize = 16;
    /// let verify = true;
    /// let numThreads = 4;
    
    let context = buildContext(matrixSize);

    computeGauss(context);

    contextToString(context)
    
