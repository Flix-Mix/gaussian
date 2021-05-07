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

    // todo: remove this
    def exit(status: Int32): () & Impure =
        import java.lang.Runtime:getRuntime();
        let runtime = getRuntime();
        import java.lang.Runtime.exit(Int32);
        runtime.exit(status)

    /**
    * The tuple corresponds to (id, channelLock, bufferSize, unbuffered, elementQueue, waitingGetters, waitingSetters)
    *   - id
    *   - channelLock is the channelLock of this channel.
    *   - bufferSize is the size of a channel. If you try to put an
    *   element in a channel that's full, you wait until there's space.
    *   - a flag for whether the channel is unbuffered or not
    *   - elementQueue is the elementQueue of elements in the list.
    *   - waitingGetters is a set of conditions that is waiting for get.
    *   This set is cleared after each new element.
    *   - waitingSetters is a condition that can notify threads of
    *   available space in the elementQueue.
    */    
    pub opaque type ChannelImpl[a] = (
        Int64,
        ReentrantLock,
        Int,
        Bool,
        MutList[a],
        MutList[LockConditionPair],
        Condition
    )


    // channel, index
    pub opaque type ChannelIndexPair[a] = (
        ChannelImpl[a],
        Int
        )

    pub def newChannelIndexPair(channel: ChannelImpl[a], index: Int): ChannelIndexPair[a] =
        ChannelIndexPair(
            channel,
            index
        )

    pub def getChannel(cip: ChannelIndexPair[a]): ChannelImpl[a] = match cip {
            case ChannelIndexPair(c, _)   => c
        }

    pub def getIndex(cip: ChannelIndexPair[a]): Int = match cip {
            case ChannelIndexPair(_, ind)   => ind
        }

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


    pub def getId(c: ChannelImpl[a]): Int64 = match c {
            case ChannelImpl(id, _, _, _, _, _, _)   => id
        }

    pub def getChannelLock(c: ChannelImpl[a]): ReentrantLock = match c {
            case ChannelImpl(_, channelLock, _, _, _, _, _)   => channelLock
        }

    pub def getBufferSize(c: ChannelImpl[a]): Int = match c {
            case ChannelImpl(_, _, bufferSize, _, _, _, _)   => bufferSize
        }

    pub def getUnbuffered(c: ChannelImpl[a]): Bool = match c {
            case ChannelImpl(_, _, _, unbuffered, _, _, _)   => unbuffered
        }

    pub def getElementQueue(c: ChannelImpl[a]): MutList[a] = match c {
            case ChannelImpl(_, _, _, _, elementQueue, _, _)   => elementQueue
        }

    // TODO : make it a MutList or MutSet??
    pub def getWaitingGetters(c: ChannelImpl[a]): MutList[LockConditionPair] = match c {
            case ChannelImpl(_, _, _, _, _, waitingGetters, _)   => waitingGetters
        }

    pub def getWaitingSetters(c: ChannelImpl[a]): Condition = match c {
            case ChannelImpl(_, _, _, _, _, _, waitingSetters)   => waitingSetters
        }


    /*
    instance Order[Channel] {
        pub def compare(x: Channel, y: Channel): Comparison =
            compare(x.id, y.id)
    }
    */

    /// instance Eq[Channel] {
    ///     pub def eq(x: Channel[a], y: Channel[b]): Bool =
    ///         channelId(x) == channelId(y)
    /// }

    /// def channelId(channel: Channel[a]): Int =
    ///     let Channel(c) = channel;
    ///     getId(c)

    type alias ReentrantLock = ##java.util.concurrent.locks.ReentrantLock

    /// type alias Lock = ##java.util.concurrent.locks.ReentrantLock

    type alias LockConditionPair = (ReentrantLock, Condition) 

    pub def getLock(pair: LockConditionPair): ReentrantLock = match pair {
            case (lock, _)   => lock
        }

    pub def getCondition(pair: LockConditionPair): Condition = match pair {
            case (_, cond)   => cond
        }

    def java_signalAll(cond: Condition): Unit & Impure =
        import java.util.concurrent.locks.Condition.signalAll() as signalAll;
        signalAll(cond);
        ()
    
    def java_lock(reentrantLock: ReentrantLock): Unit & Impure =
        import java.util.concurrent.locks.ReentrantLock.lock() as lock;
        lock(reentrantLock);
        ()

    def java_unlock(reentrantLock: ReentrantLock): Unit & Impure =
        import java.util.concurrent.locks.ReentrantLock.unlock() as unlock;
        unlock(reentrantLock);
        ()
        
    def java_await(condition: Condition): Unit & Impure =
        import java.util.concurrent.locks.Condition.await();
        condition.await()

    def newReentrantLock(): ReentrantLock & Impure =
        import new java.util.concurrent.locks.ReentrantLock() as javaNewReentrantLock;
        javaNewReentrantLock()

    type alias Condition = ##java.util.concurrent.locks.Condition

    def newCondition(reentrantLock: ReentrantLock): Condition & Impure =
        import java.util.concurrent.locks.ReentrantLock.newCondition() as javaNewCondition;
        javaNewCondition(reentrantLock)

    /// def newId(): Int64 & Impure =
    ///     import library.Channel.GlobalCounter:newId() as javaNewId;
    ///     javaNewId()

    ///
    /// Returns a fresh unbuffered channel.
    ///
    /// An unbuffered channel has an internal capacity of `0`.
    ///
    pub def new(): ChannelImpl[a] & Impure = newWithCapacity(0)

    ///
    /// Returns a fresh buffered channel with internal capacity `c`.
    ///
    pub def newWithCapacity(bufferSize: Int): ChannelImpl[a] & Impure =

///   TODO  This does not do what you think it does :)
/// Here is what you want:
/// let unbuffered = ref false;
/// ...
/// unbuffered := true

        let unbuffered = false;
        if (bufferSize < 0) {
            // TODO bug!("foo")
            Console.printLine("Channel bufferSize must be positive");
            exit(1)
        } else if (bufferSize == 0) {
            ref bufferSize := 1;
            ref unbuffered := true
        }
        else
            ref unbuffered := false;

        let mutList = MutList.new();
        MutList.reserve!(bufferSize, mutList);

        let reentrantLock = newReentrantLock();

        ChannelImpl(
            0i64,
            reentrantLock,
            bufferSize,
            unbuffered,
            mutList,
            MutList.new(), // MutSet.empty()
            newCondition(reentrantLock)
        )
    ///
    /// Receives an element from the channel `c`.
    ///
    /// Equivalent to the expression `<- c`.
    ///
    /// pub def get(c: Channel[a]): a & Impure = <- c

    pub def put_waitOnFullChannel(c: ChannelImpl[a], bufferSize: Int): () & Impure = 
        // If channel not full, return
        if (MutList.length(getElementQueue(c)) < bufferSize) ()
        // Else, wait
        else {
        java_await(getWaitingSetters(c));
        put_waitOnFullChannel(c, bufferSize)
        }


    ///
    /// Sends the element `x` on the channel `c`.
    ///
    /// Equivalent to the expression `c <- x`.
    ///
    pub def put(c: ChannelImpl[a], x: a): ChannelImpl[a] & Impure =
        lock(c);

        // Check if the channel is full
        put_waitOnFullChannel(c, getBufferSize(c));

        // There was space to put another element in the channel
        MutList.push!(x, getElementQueue(c));

        // Signal waitingGetters that there is an element available
        MutList.foreach(
            (pair) -> {
                let conditionLock = getLock(pair);
                java_lock(conditionLock);
                signalAll(pair);
                java_unlock(conditionLock)
            },
            getWaitingGetters(c));

        // Clear waitingGetters.
        // If a waitingGetter does not receive an element, it can add itself again
        MutList.clear!(getWaitingGetters(c));

        // If the channel is unbuffered, wait for the element to be handed off before continuing
        if (getUnbuffered(c)) java_await(getWaitingSetters(c))
        else ();

        unlock(c);
        c

    pub def threadInterrupted(): Bool & Impure = 
        /// import java.lang.Thread:currentThread();
        /// let thread = currentThread();
        import java.lang.Thread:interrupted();
        let isInterrupted = interrupted();
        if (isInterrupted) true
        else false


    pub def select_helper(channels: Array[ChannelImpl[a]], sortedChannels: Array[ChannelImpl[a]], selectLock: ReentrantLock, condition: Condition, hasDefault: Bool): Option[(Int, a, Bool)] & Impure =
        // MutList used to persist the Select Choice from the inner scope
        let selectChoices = MutList.new();

        if (threadInterrupted()) {Console.printLine("Thread is interrupted"); None}
        else {
            // Lock all channels in sorted order
            lockAllChannels(sortedChannels);

            // Lock the select lock after the channels
            java_lock(selectLock);

            // Find channels with waiting elements in a random order to prevent backpressure.
            // Build list mapping a channel to it's branchNumber (the index of the 'channels' array)
            let channelIndexPairs = MutList.new();
            Array.foreachWithIndex(
                (c, i) -> MutList.push!(newChannelIndexPair(c, i), channelIndexPairs),
                channels
                );  


            // TODO : Randomize the order channels are looked at.
            // This prevents backpressure from building up on one channel.
            /// Collections.shuffle(channelIndexPairs, ThreadLocalRandom.current());

            // Find channels with waiting elements
            MutList.foreachWithIndex(
                (cip, _) -> {
                    let e_opt = tryGet(getChannel(cip));

                    match e_opt {
                        case None => ()
                        case Some(element) => {
                            // There is a waiting element in this channel.
                            // Return the element and the branchNumber of this channel

                            // If an element found, add a new selectChoice to the MutList in order to later return
                            let selectChoice = Some((getIndex(cip), element, false));
                            MutList.push!(selectChoice, selectChoices)
                        }
                    }
                },
                channelIndexPairs
                );


            // Element found
            if (!MutList.isEmpty(selectChoices)) {
                /// Console.printLine("Selected an element");
                
                MutList.get(0, selectChoices)}
            // No element was found but if there is a default case, choose this
            else if (hasDefault) {
                /// Console.printLine("Default selection");
                // TODO what exactly to return?
                None
            } 
            // No element was found.
            else {
                /// Console.printLine("Nothing selected");
                // Add our condition to all channels to get notified when a new element is added
                Array.foreach(
                    (c) -> {
                        addGetter(c, selectLock, condition)
                    },
                    channels);

                // Unlock all channels in sorted order, so other threads may input elements
                unlockAllChannels(sortedChannels);

                // Wait for an element to be added to any of the channels
                java_await(condition);

                // Unlock the selectLock, which is relevant when a different thread wants to put
                // an element into a channel that was not selected from the select.
                // This other channel will then signal the condition from selectLock (in the put method),
                // so it needs the lock.
                java_unlock(selectLock);

                select_helper(channels, sortedChannels, selectLock, condition, hasDefault)
            }
        }

  /**
   * Given a array of channels, returns the first channel that has an element
   * and return the index of that channel and the retrieved element in a
   * SelectChoice object.
   *
   * @param channels the channels to select on
   * @return the channel index of the channel with an element and the element
   */
   // TODO : resolve the function name issue
    pub def selectImpl(channels: Array[ChannelImpl[a]], hasDefault: Bool): Option[(Int, a, Bool)] & Impure =
        // Create new Condition and channelLock the current thread
        let selectLock = newReentrantLock();
        let condition = newCondition(selectLock);

        // Sort channels to avoid deadlock when locking
        let sortedChannels = sortChannels(channels);
        select_helper(channels, sortedChannels, selectLock, condition, hasDefault)


    def cmp(a: ChannelImpl[a], b: ChannelImpl[a]): Int32=
        let a_id = getId(a);
        let b_id = getId(b);
        if (a_id == b_id) (0)
        else if (a_id > b_id) (1)
        else (-1)

    def sortChannels(channels: Array[ChannelImpl[a]]): Array[ChannelImpl[a]] & Impure=
        Array.sortWith(cmp, channels)

    def lock(chan: ChannelImpl[a]): Unit & Impure=
        let lock = getChannelLock(chan);
        java_lock(lock);
        ()

    def unlock(chan: ChannelImpl[a]): Unit & Impure= 
        let lock = getChannelLock(chan);
        java_unlock(lock);
        ()


    def signalAll(pair: LockConditionPair): Unit & Impure=
        let condition = getCondition(pair);
        Console.printLine("java_signalAll called");
        java_signalAll(condition)

    def lockAllChannels(channels: Array[ChannelImpl[a]]): Unit & Impure=
        Array.foreach(lock,channels);
        ()

    def unlockAllChannels(channels: Array[ChannelImpl[a]]): Unit & Impure=
        Array.foreach(unlock,channels);
        ()


    /*
    * Recursive helper function for get, loops and awaits on empty channel
    */
    def get_waitOnEmptyChannel(chan: ChannelImpl[a], e: Option[a]): Option[a] & Impure=
        let queue = getElementQueue(chan);

        match e {
            case None => { // Null
                Console.printLine("Empty channel");

                // Create a new Lock and Condition
                let conditionLock = newReentrantLock();
                java_lock(conditionLock);
                let condition = newCondition(conditionLock);

                // Create a new LockConditionPair
                let pair = (conditionLock, condition);
                let getters = getWaitingGetters(chan);

                // Add LockConditionPair to the channel
                MutList.push!(pair, getters);

                // Temporarily unlock the channel while waiting. This is necessary as the Condition comes from a different Lock.
                unlock(chan);
                java_await(condition);
                lock(chan);

                // TODO :  will always pop LIFO, should pop randomly?
                // Someone signalled that an element was put in the channel.
                // Try to get the element (which could already be taken by someone else)
                ref e := MutList.pop!(queue);

                java_unlock(conditionLock); 
                get_waitOnEmptyChannel(chan, e);

                // Return the element to get
                e
            }
            case _ => { //Not Null
                // Return the element to get
                e
            }
        }
            
    

    pub def get(chan: ChannelImpl[a]): a & Impure=
        lock(chan);
        let queue = getElementQueue(chan);

        let e_opt = MutList.pop!(queue);
        get_waitOnEmptyChannel(chan, e_opt);

        // Unwrap the Optional type element
        let e_mutList = List.toMut(Option.toList(e_opt));
        let e = MutList.get(0, e_mutList);

        // Signal waiting setters that the channel has space
        let waitingSetters = getWaitingSetters(chan);
        java_signalAll(waitingSetters);
        unlock(chan);

        // Return the element from the channel
        e
    // TODO : test putting on full channel



    def tryGet(chan: ChannelImpl[a]): Option[a] & Impure=
        lock(chan);
        let queue = getElementQueue(chan);
        let element = MutList.pop!(queue);

        match element {
            case None => ()
            case _ => {
                let ws = getWaitingSetters(chan);
                java_signalAll(ws)
            }
        };
        
        unlock(chan);
        element


    /**
    * Adds the given condition to the list of conditions waiting to
    * retrieve elements from the elementQueue.
    *
    * @param condition the condition to add
    */
    pub def addGetter(channel: ChannelImpl[a], conditionLock: ReentrantLock, condition: Condition): Unit & Impure=
        lock(channel);
        let pair = (conditionLock, condition);
        // Add LockConditionPair to the channel
        MutList.push!(pair, getWaitingGetters(channel));
        unlock(channel)



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

/**
* Print formatted string (str+val)
* @Params
*  Type.toString
*  string
*  var of Type
*/
/// def printf[a, e](toString: a -> String & e, str: String, val: a): Unit & Impure =
///     let stringBuilder = StringBuilder.new();
///     let convertedStr = toString(val);
///     StringBuilder.appendString!(stringBuilder, str);
///     StringBuilder.appendString!(stringBuilder, convertedStr);
///     let fullStr = StringBuilder.toString(stringBuilder);
///     Console.printLine(fullStr)

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

def computeGaussJLoopParallel(context: Context, i: Int32, channel: ChannelImpl[Int32], latch: ##java.util.concurrent.CountDownLatch, processId: Int32): () & Impure =
    /// import java.util.concurrent.CountDownLatch.getCount();
    import java.util.concurrent.CountDownLatch.countDown();

    let channels = [channel];

    match selectImpl(channels, true) {
        case Some((_, j, _))   => {
            latch.countDown();
            printf(Int32.toString, "j ", j);
            let pivotValue = (context.matrix)[j][i];
            (context.matrix)[j][i] = 0.0;
            computeGaussKLoop(context, pivotValue, i, j, i + 1);
            (context.b)[j] = (pivotValue * (context.b)[i]) - (context.b)[j];
            /// printf(Int64.toString, "Counted down, now ", latch.getCount());
            computeGaussJLoopParallel(context, i, channel, latch, processId)
        }
            ///     1 => {printf(Int32.toString, "Selected channel ", cInd); el}
            ///     2 => {printf(Int32.toString, "Selected channel ", cInd); el}
            ///     3 => {printf(Int32.toString, "Selected channel ", cInd); el}
            /// }
        case None => ()
    }


    /// select {
    ///     case j <- channel => {
    ///         latch.countDown();
    ///         /// printf(Int32.toString, "j ", j);
    ///         let pivotValue = (context.matrix)[j][i];
    ///         (context.matrix)[j][i] = 0.0;
    ///         computeGaussKLoop(context, pivotValue, i, j, i + 1);
    ///         (context.b)[j] = (pivotValue * (context.b)[i]) - (context.b)[j];
    ///         /// printf(Int64.toString, "Counted down, now ", latch.getCount());
    ///         computeGaussJLoopParallel(context, i, channel, latch, processId)
    ///     }
    ///     case _ => ()
    ///     /// printf(Int32.toString, "In default of process ", processId)
    /// }


def fillChannel(context: Context, channel: ChannelImpl[Int32], j: Int32, matrixSize: Int32): () & Impure =
    if (j == matrixSize) ()
    else {
        put(channel, j);
        /// channel <- j;
        fillChannel(context, channel, j+1, matrixSize)
    }

def startThreads(context: Context, i: Int32, channel: ChannelImpl[Int32], numThreads: Int32, latch: ##java.util.concurrent.CountDownLatch, processId: Int32): () & Impure =
    if (numThreads == 0) ()
    else {
        spawn computeGaussJLoopParallel(context, i, channel, latch, processId);
        startThreads(context, i, channel, numThreads-1, latch, processId+1)
    }


def computeGauss(context: Context, matrixSize: Int32, numThreads: Int32): () & Impure =
    let computeGaussILoop = (i: Int32) -> {
        getPivot(context, i, matrixSize);   

		let rowsLeft = matrixSize - i - 1;
        let c = newWithCapacity(rowsLeft);        
        fillChannel(context, c, i+1, matrixSize);

        import new java.util.concurrent.CountDownLatch(Int32) as latch;
        import java.util.concurrent.CountDownLatch.await();
        let signal = latch(rowsLeft);
        startThreads(context, i, c, numThreads, signal, 0);
        signal.await()
        /// Console.printLine("Barrier. Next iteration");
        /// Console.printLine("")
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
    import java.lang.System:currentTimeMillis();
    let matrixSize = 16;
    let verify = true;
    let printContext = true;
    let numThreads = 2;
    let context = buildContext(matrixSize);
    let start = currentTimeMillis();

    computeGauss(context, matrixSize, numThreads);

    let finish = currentTimeMillis();
    let timeElapsed = Int64.toFloat64(finish - start);
    /// let timeElapsed_sec = Int64.toFloat64((finish - start))/1000.0;
    /// printf(Int32.toString, "Matrix size: ", matrixSize);
    /// printf(Int32.toString, "Number of threads: ", numThreads);
    /// printf(Float64.toString, "", timeElapsed);

	if (verify) ref context := solveGauss(context, matrixSize)
	else ();

    if (printContext) Console.printLine(contextToString(context))
	else ();

	if (verify) {
		verifyLoop(context, 0, matrixSize)
	} else ();
    

    let stringBuilder = StringBuilder.new();
    let convertedStr = Float64.toString(timeElapsed);
    StringBuilder.appendString!(stringBuilder, convertedStr);
    let fullStr = StringBuilder.toString(stringBuilder);
    fullStr
