/* -------------------------------------------------------------------------- */
/*           Gaussian Elimination using Java channel implementation           */
/* -------------------------------------------------------------------------- */

namespace GaussParallelJava {

    def computeGaussKLoop(context: Context, pivotValue: Float64, i: Int32, j: Int32, k: Int32): () & Impure =
        if (k < Array.length((context.matrix)[i])) {
            (context.matrix)[j][k] = (pivotValue * (context.matrix)[i][k]) - (context.matrix)[j][k];
            computeGaussKLoop(context, pivotValue, i, j, k + 1)
        } else ()

    def computeGaussJLoopParallel(context: Context, i: Int32, channel: Channel[Int32], latch: ##java.util.concurrent.CountDownLatch, processId: Int32): () & Impure =
        import java.util.concurrent.CountDownLatch.countDown();

        select {
            case j <- channel => {
                latch.countDown();
                let pivotValue = (context.matrix)[j][i];
                (context.matrix)[j][i] = 0.0;
                computeGaussKLoop(context, pivotValue, i, j, i + 1);
                (context.b)[j] = (pivotValue * (context.b)[i]) - (context.b)[j];
                computeGaussJLoopParallel(context, i, channel, latch, processId)
            }
            case _ => ()
        }


    def fillChannel(context: Context, channel: Channel[Int32], j: Int32, matrixSize: Int32): () & Impure =
        if (j == matrixSize) ()
        else {
            channel <- j;
            fillChannel(context, channel, j+1, matrixSize)
        }

    def startThreads(context: Context, i: Int32, channel: Channel[Int32], numThreads: Int32, latch: ##java.util.concurrent.CountDownLatch, processId: Int32): () & Impure =
        if (numThreads == 0) ()
        else {
            spawn computeGaussJLoopParallel(context, i, channel, latch, processId);
            startThreads(context, i, channel, numThreads-1, latch, processId+1)
        }


    pub def computeGauss(context: Context, matrixSize: Int32, numThreads: Int32): () & Impure =
        let computeGaussILoop = (i: Int32) -> {
            getPivot(context, i, matrixSize);   

            let rowsLeft = matrixSize - i - 1;
            let c = Channel.newWithCapacity(rowsLeft);        
            fillChannel(context, c, i+1, matrixSize);

            import new java.util.concurrent.CountDownLatch(Int32) as latch;
            import java.util.concurrent.CountDownLatch.await();
            let signal = latch(rowsLeft);
            startThreads(context, i, c, numThreads, signal, 0);
            signal.await()
        };
        Array.foreachWithIndex((_, i) -> computeGaussILoop(i), context.matrix)    


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

    pub def solveGauss(context: Context, matrixSize: Int32): Context & Impure =
        (context.v)[matrixSize-1] = (context.b)[matrixSize-1];
        ref context := solveGaussVILoop(context, matrixSize-2, matrixSize);
        solveGaussCILoop(context, 0, matrixSize)

    pub def verifyLoop(context: Context, i: Int32, matrixSize: Int32): Unit & Impure =
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

}
