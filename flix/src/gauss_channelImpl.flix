/*
 * Copyright 2019 Magnus Madsen
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/**
* Print formatted string (str+val)
* @Params
*  Type.toString
*  string
*  var of Type
*/
def printf[a, e](toString: a -> String & e, str: String, val: a): Unit & Impure =
    let stringBuilder = StringBuilder.new();
    let convertedStr = toString(val);
    StringBuilder.appendString!(stringBuilder, str);
    StringBuilder.appendString!(stringBuilder, convertedStr);
    let fullStr = StringBuilder.toString(stringBuilder);
    Console.printLine(fullStr)




/* -------------------------------------------------------------------------- */
/*                       Circular Buffer implementation                       */
/* -------------------------------------------------------------------------- */

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




/* -------------------------------------------------------------------------- */
/*                           Channel implementation                           */
/* -------------------------------------------------------------------------- */


/**
* The tuple corresponds to (id, channelLock, bufferSize, unbuffered, elementQueue, waitingGetters, waitingSetters)
*   - id
*   - channelLock is the channelLock of this channel.
*   - bufferSize is the size of a channel. If you try to put an element in a channel that's full, you wait until there's space.
*   - a flag for whether the channel is unbuffered or not
*   - elementQueue is a circular buffer of elements in the list.
*   - waitingGetters is a set of conditions that is waiting for get. This set is cleared after each new element.
*   - waitingSetters is a condition that can notify threads of available space in the elementQueue.
*/

/* ---------------------------- Type declarations --------------------------- */

pub enum ChannelImpl[a] {
    case Mpmc(
            Int64,
            ReentrantLock,
            Int,
            Bool,
            CircularBuffer[a],
            MutList[(ReentrantLock, Condition)],
            Condition
        )
}

type alias Condition = ##java.util.concurrent.locks.Condition

type alias ReentrantLock = ##java.util.concurrent.locks.ReentrantLock

/* ----------------------------- Implementation ----------------------------- */

namespace ChannelImpl {


/* ------------------------------ Java wrappers ----------------------------- */

    def newReentrantLock(): ReentrantLock & Impure =
        import new java.util.concurrent.locks.ReentrantLock() as javaNewReentrantLock;
        javaNewReentrantLock()

    def newCondition(reentrantLock: ReentrantLock): Condition & Impure =
        import java.util.concurrent.locks.ReentrantLock.newCondition() as javaNewCondition;
        javaNewCondition(reentrantLock)

    def java_lock(reentrantLock: ReentrantLock): Unit & Impure =
        import java.util.concurrent.locks.ReentrantLock.lock() as javaLock;
        javaLock(reentrantLock)

    def java_unlock(reentrantLock: ReentrantLock): Unit & Impure =
        import java.util.concurrent.locks.ReentrantLock.unlock() as javaUnlock;
        javaUnlock(reentrantLock)

    def signalAll(cond: Condition): Unit & Impure =
        import java.util.concurrent.locks.Condition.signalAll() as javaSignalAll;
        javaSignalAll(cond)

    def await(condition: Condition): Unit & Impure =
        import java.util.concurrent.locks.Condition.await() as javaAwait;
        javaAwait(condition)

    def lock(chan: ChannelImpl[a]): Unit & Impure =
        let Mpmc(_, lock, _, _, _, _, _) = chan;
        java_lock(lock)

    def unlock(chan: ChannelImpl[a]): Unit & Impure =
        let Mpmc(_, lock, _, _, _, _, _) = chan;
        java_unlock(lock)

    def lockAllChannels(channels: Array[ChannelImpl[a]]): Unit & Impure =
        Array.foreach(lock, channels)

    def unlockAllChannels(channels: Array[ChannelImpl[a]]): Unit & Impure =
        Array.foreach(unlock, channels)

    def threadInterrupted(): Bool & Impure =
        import java.lang.Thread:interrupted();
        let isInterrupted = interrupted();
        if (isInterrupted)
            true
        else
            false

    // TODO : implement equals, hashCode, Order, Eq

    /// @Override
    /// pub def equals(o: a): Bool & Impure =
    ///     if (this == o) return true;
    ///     if (o == null || getClass() != o.getClass()) return false;
    ///     ChannelIndexPair that = (ChannelIndexPair) o;
    ///     return index == that.index && Objects.equals(channel, that.channel);
    /// }

    /// @Override
    /// pub def hashCode(cip: ChannelIndexPair): Int =
    ///     import java.util.Objects.hash as hash;
    ///     hash(getChannel(cip), getIndex(cip))

    /// instance Order[ChannelImpl[a]] {
    ///     pub def compare(x: Channel, y: Channel): Comparison =
    ///         compare(x.id, y.id)
    /// }

    /// instance Eq[ChannelImpl[a]] {
    ///     pub def eq(x: Channel, y: Channel[b]): Bool =
    ///         channelId(x) == channelId(y)
    /// }


    // TODO implement global counter for id
    def newId(): Int64 =
        /// import library.Channel.GlobalCounter:newId() as javaNewId;
        /// javaNewId()
        1i64


    def bug!(message: String): Unit & Impure = 
        Console.printLine(message)


    ///
    /// Returns a fresh buffered channel with internal capacity `c`.
    ///
    pub def new(bufferSize: Int): ChannelImpl[a] & Impure =
        let (bufSize, unbuffered) =
            if (bufferSize < 0) {
                bug!("Channel's buffer size must be non-negative.");
                (1, false)
            } else if (bufferSize == 0)
                (1, true)
            else
                (bufferSize, false);

        let reentrantLock = newReentrantLock();

        Mpmc(
            newId(),
            reentrantLock,
            bufSize,
            unbuffered,
            CircularBuffer.new(bufferSize),
            MutList.new(), 
            newCondition(reentrantLock)
        )


    ///
    /// Sends the element `x` on the channel `c`.
    ///
    /// Equivalent to the expression `c <- x`.
    ///
    pub def put_dup(c: ChannelImpl[a], x: a): ChannelImpl[a] & Impure =
        lock(c);

        // Check if the channel is full
        putHelperWaitOnFullChannel(c);

        // There was space to put another element in the channel
        let Mpmc(_, _, _, unbuffered, queue, waitingGetters, waitingSetters) = c;
        CircularBuffer.put(x, queue);

        // Signal waitingGetters that there is an element available
        MutList.foreach(
            (lockConditionPair) -> {
                let (conditionLock, condition) = lockConditionPair;
                java_lock(conditionLock);
                signalAll(condition);
                java_unlock(conditionLock)
            },
            waitingGetters);

        // Clear waitingGetters.
        // If a waitingGetter does not receive an element, it can add itself again
        MutList.clear!(waitingGetters);

        // If the channel is unbuffered, wait for the element to be handed off before continuing
        if (unbuffered) 
            await(waitingSetters)
        else 
            ();

        unlock(c);

        // Return the channel
        c


    def putHelperWaitOnFullChannel(c: ChannelImpl[a]): () & Impure =
        let Mpmc(_, _, bufferSize, _, queue, _, waitingSetters) = c;
        if (CircularBuffer.length(queue) >= bufferSize) {
            await(waitingSetters);
            putHelperWaitOnFullChannel(c)
        } else {
            ()
        }


    ///
    /// Receives an element from the channel `c`.
    ///
    /// Equivalent to the expression `<- c`.
    ///
    pub def get_dup(chan: ChannelImpl[a]): a & Impure =
        lock(chan);

        let Mpmc(_, _, _, _, queue, _, waitingSetters) = chan;
        let optionalElement = CircularBuffer.tryGet(queue);
        let element = getHelperWaitOnEmptyChannel(chan, optionalElement);

        // Signal waiting setters that the channel has space
        signalAll(waitingSetters);

        unlock(chan);

        // Return the element from the channel
        element


    /*
    * Recursive helper function for get, loops and awaits on empty channel
    */
    def getHelperWaitOnEmptyChannel(chan: ChannelImpl[a], element: Option[a]): a & Impure =
        match element {
            case None => {
                let Mpmc(_, _, _, _, queue, waitingGetters, _) = chan;
                // No element was found

                // Create a new Lock and Condition
                let conditionLock = newReentrantLock();
                java_lock(conditionLock);

                let condition = newCondition(conditionLock);
                // Add LockConditionPair to the channel
                MutList.push!((conditionLock, condition), waitingGetters);

                // Temporarily unlock the channel while waiting. This is necessary as the Condition comes from a different Lock.
                unlock(chan);
                await(condition);
                lock(chan);

                // Someone signalled that an element was put in the channel.
                // Try to get the element (which could already be taken by someone else)
                let optionalElement = CircularBuffer.tryGet(queue);

                java_unlock(conditionLock);

                getHelperWaitOnEmptyChannel(chan, optionalElement)
            }
            case Some(e) => e
        }


    pub def tryGet(chan: ChannelImpl[a]): Option[a] & Impure =
        lock(chan);

        let Mpmc(_, _, _, _, queue, _, waitingSetters) = chan;
        // Try to get an element from the channel
        let element = CircularBuffer.tryGet(queue);

        // If there was an element, signal waiting setters
        match element {
            case None => ()
            case _ => {
                let ws = waitingSetters;
                signalAll(ws)
            }
        };

        unlock(chan);
        // Return the element from the channel, or None if channel was empty
        element


    def cmp(a: ChannelImpl[a], b: ChannelImpl[a]): Int32 =
        let Mpmc(a_id, _, _, _, _, _, _) = a;
        let Mpmc(b_id, _, _, _, _, _, _) = b;
        if (a_id == b_id)
            (0)
        else if (a_id > b_id)
            (1)
        else
            (-1)


    def sortChannels(channels: Array[ChannelImpl[a]]): Array[ChannelImpl[a]] & Impure =
        Array.sortWith(cmp, channels)


  /**
   * Given a array of channels, returns the first channel that has an element
   * and return the index of that channel and the retrieved element in a
   * SelectChoice object.
   *
   * @param channels the channels to select on
   * @return the channel index of the channel with an element and the element
   */
    pub def selectImpl(channels: Array[ChannelImpl[a]], hasDefault: Bool): Option[(Int, a)] & Impure =
        // Create new Condition and channelLock the current thread
        let selectLock = newReentrantLock();
        let condition = newCondition(selectLock);

        // Sort channels to avoid deadlock when locking
        let sortedChannels = sortChannels(channels);
        selectHelper(channels, sortedChannels, selectLock, condition, hasDefault)

        
    def selectHelper(channels: Array[ChannelImpl[a]], sortedChannels: Array[ChannelImpl[a]], selectLock: ReentrantLock, condition: Condition, hasDefault: Bool): Option[(Int, a)] & Impure =
        if (not threadInterrupted()) {
            // Lock all channels in sorted order
            lockAllChannels(sortedChannels);

            // Lock the select lock after the channels
            java_lock(selectLock);

            // Find channels with waiting elements
            let selectChoicesWithWaitingElements = MutList.new();
            Array.foreachWithIndex(
                (c, i) -> {
                    match tryGet(c) {
                        case None => ()
                        case Some(element) => {
                            // There is a waiting element in this channel.
                            // Return the element and the branchNumber of this channel

                            // If an element found, add a new selectChoice to the MutList in order to later return
                            MutList.push!((i, element), selectChoicesWithWaitingElements)
                        }
                    }
                },
                channels
                );

            if (not MutList.isEmpty(selectChoicesWithWaitingElements)) {
                // Element found
                Some(MutList.get(0, selectChoicesWithWaitingElements))
            }
            else if (hasDefault) {
                // No element was found and default case exists
                None
            } else {
                // No element was found.
                // Add our condition to all channels to get notified when a new element is added
                Array.foreach(
                    (c) -> {
                        addGetter(c, selectLock, condition)
                    },
                    channels);

                // Unlock all channels in sorted order, so other threads may input elements
                unlockAllChannels(sortedChannels);

                // Wait for an element to be added to any of the channels
                await(condition);

                // Unlock the selectLock, which is relevant when a different thread wants to put
                // an element into a channel that was not selected from the select.
                // This other channel will then signal the condition from selectLock (in the put method),
                // so it needs the lock.
                java_unlock(selectLock);

                selectHelper(channels, sortedChannels, selectLock, condition, hasDefault)
            }
        } else {
            bug!("Thread is interrupted");
            None
        }


    /**
    * Adds the given condition to the list of conditions waiting to
    * retrieve elements from the elementQueue.
    *
    * @param condition the condition to add
    */
    def addGetter(channel: ChannelImpl[a], conditionLock: ReentrantLock, condition: Condition): Unit & Impure =
        lock(channel);
        let pair = (conditionLock, condition);

        // Add LockConditionPair to the channel
        let Mpmc(_, _, _, _, _, waitingGetters, _) = channel;
        MutList.push!(pair, waitingGetters);
        unlock(channel)
}





/* -------------------------------------------------------------------------- */
/*                            Gaussian Elimination                            */
/* -------------------------------------------------------------------------- */



/* --------------------- Context and formatting methods --------------------- */

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



/* ----------------------------- Implementation ----------------------------- */

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


def computeGaussKLoop(context: Context, pivotValue: Float64, i: Int32, j: Int32, k: Int32): () & Impure =
    if (k < Array.length((context.matrix)[i])) {
        (context.matrix)[j][k] = (pivotValue * (context.matrix)[i][k]) - (context.matrix)[j][k];
        computeGaussKLoop(context, pivotValue, i, j, k + 1)
    } else ()


def computeGaussJLoopParallel(context: Context, i: Int32, channel: ChannelImpl[Int], latch: ##java.util.concurrent.CountDownLatch, processId: Int32): () & Impure =
    import java.util.concurrent.CountDownLatch.countDown();

    match ChannelImpl.tryGet(channel) {
        case Some(j)   => {
            latch.countDown();
            let pivotValue = (context.matrix)[j][i];
            (context.matrix)[j][i] = 0.0;
            computeGaussKLoop(context, pivotValue, i, j, i + 1);
            (context.b)[j] = (pivotValue * (context.b)[i]) - (context.b)[j];
            computeGaussJLoopParallel(context, i, channel, latch, processId)
        }
        case None => ()
    }


def fillChannel(context: Context, channel: ChannelImpl[Int], j: Int32, matrixSize: Int32): () & Impure =
    if (j == matrixSize) ()
    else {
        ChannelImpl.put_dup(channel, j);
        fillChannel(context, channel, j+1, matrixSize)
    }


def startThreads(context: Context, i: Int32, channel: ChannelImpl[Int], numThreads: Int32, latch: ##java.util.concurrent.CountDownLatch, processId: Int32): () & Impure =
    if (numThreads == 0) ()
    else {
        spawn computeGaussJLoopParallel(context, i, channel, latch, processId);
        startThreads(context, i, channel, numThreads-1, latch, processId+1)
    }


def computeGauss(context: Context, matrixSize: Int32, numThreads: Int32): () & Impure =
    let computeGaussILoop = (i: Int32) -> {
        getPivot(context, i, matrixSize);   

		let rowsLeft = matrixSize - i - 1;
        let c = ChannelImpl.new(rowsLeft);        
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
        /// Console.printLine("The matrix is singular\n");
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





/* -------------------------------------------------------------------------- */
/*                                    Main                                    */
/* -------------------------------------------------------------------------- */

def main(): String & Impure =

    import java.lang.System:currentTimeMillis();
    let matrixSize = 128;
    let verify = true;
    let printContext = false;
    let numThreads = 4;
    let context = buildContext(matrixSize);
    let start = currentTimeMillis();

    computeGauss(context, matrixSize, numThreads);

    let finish = currentTimeMillis();
    let timeElapsed = Int64.toFloat64(finish - start);

    printf(Int32.toString, "Matrix size: ", matrixSize);
    printf(Int32.toString, "Number of threads: ", numThreads);
    printf(Float64.toString, "Time elapsed [ms]: ", timeElapsed);

	if (verify) ref context := solveGauss(context, matrixSize)
	else ();

    if (printContext) Console.printLine(contextToString(context))
	else ();

	if (verify) {
		verifyLoop(context, 0, matrixSize)
	} else ();
    ""