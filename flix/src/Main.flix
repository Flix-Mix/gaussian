def printf[a, e](toString: a -> String & e, str: String, val: a): Unit & Impure =
    let stringBuilder = StringBuilder.new();
    let convertedStr = toString(val);
    StringBuilder.appendString!(stringBuilder, str);
    StringBuilder.appendString!(stringBuilder, convertedStr);
    let fullStr = StringBuilder.toString(stringBuilder);
    Console.printLine(fullStr)

// The main entry point.
def main(_args: Array[String]): Int32 & Impure =
  // The following variables would be command line arguments, but Flix doesn't
  // have the ability to recieve command line arguments
  import java.lang.System:currentTimeMillis();
  //let matrixSize = 128;
  let verify = false;
  let printContext = false;
  //let numThreads = 1;
  let matrixSizes = [128, 256, 512, 1024, 2048];
  let threadCounts = [1, 2, 4, 8, 16, 32, 64, 72];
  // let context = buildContext(128);
  Array.foreach((n) -> {
      Array.foreach((t) -> {
      let context = buildContext(n);
      let start = currentTimeMillis();
      //printf(Int32.toString, "Matrix size: ", n);
      //printf(Int32.toString, "Number of threads: ", t);
      computeGauss(context, n, t);
      let finish = currentTimeMillis();
      let timeElapsed = Int64.toFloat64(finish - start);
      printf(Float64.toString, "", timeElapsed);

    if (verify) ref context := solveGauss(context, 128)
    else ();

      if (printContext) Console.printLine(contextToString(context))
    else ();

      if (verify) {
        verifyLoop(context, 0, 128)
    } else ()
  }, threadCounts)}, matrixSizes);

  /// let timeElapsed_sec = Int64.toFloat64((finish - start))/1000.0;

  0 // exit code
