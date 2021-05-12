
def printf[a, e](toString: a -> String & e, str: String, val: a): Unit & Impure =
    let stringBuilder = StringBuilder.new();
    let convertedStr = toString(val);
    StringBuilder.appendString!(stringBuilder, str);
    StringBuilder.appendString!(stringBuilder, convertedStr);
    let fullStr = StringBuilder.toString(stringBuilder);
    Console.printLine(fullStr)


  enum Implementation {
      case Seq,
      case ParallelJava,
      case ParallelFlix
  }


// The main entry point.
def main(): Int32 & Impure =
  // The following variables would be command line arguments, but Flix doesn't
  // have the ability to recieve command line arguments
  import java.lang.System:currentTimeMillis();
  let matrixSize = 128;
  let verify = false;
  let printContext = false;
  let numThreads = 4;
  let context = buildContext(matrixSize);
  let start = currentTimeMillis();

  let impl = Implementation.ParallelFlix;

  match impl {
    case Seq => GaussSequential.computeGauss(context, matrixSize)
    case ParallelJava => GaussParallelJava.computeGauss(context, matrixSize, numThreads)
    case ParallelFlix => GaussParallelFlix.computeGauss(context, matrixSize, numThreads)
  };

  let finish = currentTimeMillis();
  let timeElapsed = Int64.toFloat64(finish - start);
  printf(Int32.toString, "Matrix size: ", matrixSize);
  printf(Int32.toString, "Number of threads: ", numThreads);
  printf(Float64.toString, "Time elapsed [ms]: ", timeElapsed);

	if (verify) {
    ref context := match impl {
      case Seq => GaussSequential.solveGauss(context, matrixSize)
      case ParallelJava => GaussParallelJava.solveGauss(context, matrixSize)
      case ParallelFlix => GaussParallelFlix.solveGauss(context, matrixSize)
    }
  }
	else ();

    if (printContext) Console.printLine(contextToString(context))
	else ();

	if (verify) {
    match impl {
      case Seq => GaussSequential.verifyLoop(context, 0, matrixSize)
      case ParallelJava => GaussParallelJava.verifyLoop(context, 0, matrixSize)
      case ParallelFlix => GaussParallelFlix.verifyLoop(context, 0, matrixSize)
    }
	} else ();

  // Exit code
  0