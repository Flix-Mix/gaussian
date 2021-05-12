def printf[a, e](toString: a -> String & e, str: String, val: a): Unit & Impure =
    let stringBuilder = StringBuilder.new();
    let convertedStr = toString(val);
    StringBuilder.appendString!(stringBuilder, str);
    StringBuilder.appendString!(stringBuilder, convertedStr);
    let fullStr = StringBuilder.toString(stringBuilder);
    Console.printLine(fullStr)


pub enum CircularBuffer[a] {
    case CircularBuffer(
        Ref[Array[a]], //buf: array storing elements of the circular buffer 
        Ref[Int], //startInd: ptr to the start of stored data
        Ref[Int], //endInd: ptr to the end of stored data
        Ref[Int] //count: current number of elements
    )
}


namespace CircularBuffer {
    ///
    /// Returns a new empty mutable list with a default capacity.
    ///
    pub def new(capacity: Int): CircularBuffer[a] & Impure = CircularBuffer(ref [default; capacity], ref 0, ref 0, ref 0)

    ///
    /// Returns the capacity of `cb`.
    ///
    pub def capacity(cb: CircularBuffer[a]): Int & Impure =
        let CircularBuffer(a, _, _, _) = cb;
        Array.length(deref a)


    ///
    /// Returns the first-in element in the circular buffer `cb`.
    ///
    pub def get(cb: CircularBuffer[a]): a & Impure =
        let CircularBuffer(a, si, _, l) = cb;
        let len = deref l;
        let startInd = deref si;
        if (isEmpty(cb)) {
            Console.printLine("Index out of bounds: the buffer is empty");
            (deref a)[-1] /// deliberately trigger index-out-of-bounds
        }
        else {
            l := len - 1;
            // Move the start index to start if at the end
            if (startInd == capacity(cb)-1) si := 0
            else si := startInd + 1;
            Array.get(deref a, startInd)
        }


    ///
    /// Returns `true` if the given circular buffer `cb` is empty.
    ///
    pub def isEmpty(cb: CircularBuffer[a]): Bool & Impure = length(cb) == 0


    ///
    /// Returns `true` if the given circular buffer `cb` is full.
    ///
    pub def isFull(cb: CircularBuffer[a]): Bool & Impure = length(cb) == capacity(cb)


    ///
    /// Returns the number of elements in the given circular buffer `cb`.
    ///
    pub def length(cb: CircularBuffer[a]): Int & Impure =
        let CircularBuffer(_, _, _, l) = cb;
        deref l


    ///
    /// Inserts the given element `x` at the end of the given circular buffer `cb` and returns true if successful.
    ///
    pub def put(x: a, cb: mut CircularBuffer[a]): Bool & Impure =
        let CircularBuffer(a, _, ei, l) = cb;
        let len = deref l;
        let endInd = deref ei;

        if (isFull(cb)) false
        else {
            (deref a)[endInd] = x;
            l := len + 1;
            // Move the end index to start if at the end
            if (endInd == capacity(cb)-1) ei := 0
            else ei := endInd + 1;
            true
        }

    ///
    /// Returns the first-in element in the circular buffer `cb` as an Optional type.
    ///
    pub def tryGet(cb: CircularBuffer[a]): Option[a] & Impure =
        let CircularBuffer(a, si, _, l) = cb;
        let len = deref l;
        let startInd = deref si;
        if (isEmpty(cb)) {
            None
        }
        else {
            l := len - 1;
            si := startInd + 1;
            Some(Array.get(deref a, startInd))
        }
}

def main(): String & Impure =

    let elementQueue = CircularBuffer.new(3);
    Console.printLine("Creating buffer of size 3");

    if (CircularBuffer.isEmpty(elementQueue)) Console.printLine("Buffer is empty") else Console.printLine("Buffer is not empty");

    let result = CircularBuffer.put(1, elementQueue);
    Console.printLine("");
    Console.printLine("Adding 1");
    if (result) Console.printLine("Put was a success") else Console.printLine("Put failed");

    let result2 = CircularBuffer.put(2, elementQueue);
    Console.printLine("");
    Console.printLine("Adding 2");
    if (result2) Console.printLine("Put was a success") else Console.printLine("Put failed");

    let result3 = CircularBuffer.put(3, elementQueue);
    Console.printLine("");
    Console.printLine("Adding 3");
    if (result3) Console.printLine("Put was a success") else Console.printLine("Put failed");
    if (CircularBuffer.isFull(elementQueue)) Console.printLine("Buffer is full") else Console.printLine("Buffer is not full");

    let element1 = CircularBuffer.get(elementQueue);
    printf(Int32.toString, "Getting element ", element1);

    let result4 = CircularBuffer.put(4, elementQueue);
    Console.printLine("");
    Console.printLine("Adding 4");
    if (result4) Console.printLine("Put was a success") else Console.printLine("Put failed");

    let element2 = CircularBuffer.get(elementQueue);
    Console.printLine("");
    printf(Int32.toString, "Getting element ", element2);
    let element3 = CircularBuffer.get(elementQueue);
    Console.printLine("");
    printf(Int32.toString, "Getting element ", element3);
    if (CircularBuffer.isEmpty(elementQueue)) Console.printLine("Buffer is empty") else Console.printLine("Buffer is not empty");
    let element4 = CircularBuffer.get(elementQueue);
    Console.printLine("");
    printf(Int32.toString, "Getting element ", element4);
    ""
