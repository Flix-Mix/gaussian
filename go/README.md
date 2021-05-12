# go

To compile and run from CL:
    $ go build GaussSequential.go
    $ ./GaussSequential -s 1024 -v


    $ go build GaussParallel.go
    $ ./GaussParallel -s 1024 -n 8 -v

Flags:
    -v      verify
    -s Int  size of matrix
    -n Int  number of threads