type alias Context = {
    matrix: Array[Array[Float64]], 
    b: Array[Float64], 
    v: Array[Float64], 
    c: Array[Float64],
    swap: Array[Int32]
}

def buildContext(matrixSize: Int32): Context & Impure = 
    {
        matrix = [[0; matrixSize]; matrixSize],
        b = [0; matrixSize],
        v = [0; matrixSize],
        c = [0; matrixSize],
        swap = [0; matrixSize]
    }

def initMatrix(context: Context): () & Impure = 
    Array.foreachWithIndex((array, i) -> 
        Array.foreachWithIndex((element, j) -> 
            element := match j < i {
                case true => 2 * (j + 1)
                case false => 2 * (i + 1)
            };
            context.b := i;
            context.swap := i, 
        array), 
    context.matrix)

        
/// def initializeMatrix(matrixSize: Int32): 

/// /// Meow's `n` times on channel `c`.
/// def meow(c: Channel[String], n: Int): Unit & Impure =
///     match n {
///         case 0 => ()
///         case x => c <- "Meow!"; meow(c, x - 1)
///     }

/// /// Hiss'es `n` times on channel `c`.
/// def hiss(c: Channel[String], n: Int): Unit & Impure =
///     match n {
///         case 0 => ()
///         case x => c <- "Hiss!"; hiss(c, x - 1)
///     }

/// Start the animal farm...
def main(): String & Impure =
    // The following variables would be command line arguments, but Flix doesn't
    // have the ability to recieve command line arguments
    let matrixSize = 1024;
    let verify = true;
    let numThreads = 4;
    
    let context = buildContext(matrixSize);

    "Hey Flix!"

    