
type alias Context = {
    matrix: Array[Array[Float64]], 
    b: Array[Float64], 
    v: Array[Float64], 
    c: Array[Float64],
    swap: Array[Int32]
}


def buildMatrix(matrix: Array[Array[Float64]], i: Int32, matrixSize: Int32): Array[Array[Float64]] & Impure = 
    if (i == matrixSize) matrix
    else {
        matrix[i] = [0.0; matrixSize];
        buildMatrix(matrix, i+1, matrixSize)
    } 


pub def buildContext(matrixSize: Int32): Context & Impure = 
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