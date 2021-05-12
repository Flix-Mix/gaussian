package main

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"time"
)

var matrix [][]float64
var B []float64
var V []float64
var C []float64

var swap_arr []int

func printArray(array []float64, nsize int) {
	fmt.Println("\nPrinting array:")
	for i := 0; i < nsize; i++ {
		s := fmt.Sprintf("%6.5f ", array[i])
		io.WriteString(os.Stdout, s)
	}
	fmt.Println()
}

func printMatrix(nsize int) {
	fmt.Println("\nPrinting matrix:")
	for i := 0; i < nsize; i++ {
		for j := 0; j < nsize; j++ {
			s := fmt.Sprintf("%6.5f ", matrix[i][j])
			io.WriteString(os.Stdout, s)
		}
		fmt.Println()
	}
}

func initMatrix(nsize int) {
	matrix = make([][]float64, nsize)
	for i := 0; i < nsize; i++ {
		matrix[i] = make([]float64, nsize)
	}
	B = make([]float64, nsize)
	V = make([]float64, nsize)
	C = make([]float64, nsize)
	swap_arr = make([]int, nsize)

	for i := 0; i < nsize; i++ {
		for j := 0; j < nsize; j++ {
			// fmt.Println("%d \n", i)
			// fmt.Println("%d \n", j)
			if j < i {
				matrix[i][j] = float64(2 * (j + 1))
			} else {
				matrix[i][j] = float64(2 * (i + 1))
			}
		}
		B[i] = float64(i)
		swap_arr[i] = i
	}

}

func getPivot(nsize int, currow int) {
	big := matrix[currow][currow]
	irow := currow

	if big == 0.0 {
		for i := currow; i < nsize; i++ {
			tmp := matrix[i][currow]
			if tmp != 0.0 {
				big = tmp
				irow = i
				break
			}
		}
	}

	if big == 0.0 {
		fmt.Printf("The matrix is singular\n")
		os.Exit(3)
	}

	if irow != currow {
		for i := currow; i < nsize; i++ {
			matrix[irow][i], matrix[currow][i] = matrix[currow][i], matrix[irow][i]
		}
		B[irow], B[currow] = B[currow], B[irow]
		swap_arr[irow], swap_arr[currow] = swap_arr[currow], swap_arr[irow]
	}

	{
		pivotVal := matrix[currow][currow]

		if pivotVal != 1.0 {
			matrix[currow][currow] = 1.0
			for i := currow + 1; i < nsize; i++ {
				matrix[currow][i] /= pivotVal
			}
			B[currow] /= pivotVal
		}
	}
}

func computeGauss(nsize int) {

	for i := 0; i < nsize; i++ {
		getPivot(nsize, i)

		pivotVal := matrix[i][i]

		for j := i + 1; j < nsize; j++ {
			pivotVal = matrix[j][i]
			matrix[j][i] = 0.0
			for k := i + 1; k < nsize; k++ {
				matrix[j][k] -= pivotVal * matrix[i][k]
			}
			B[j] -= pivotVal * B[i]
		}
	}

}

func solveGauss(nsize int) {

	V[nsize-1] = B[nsize-1]
	for i := nsize - 2; i >= 0; i-- {
		V[i] = B[i]
		for j := nsize - 1; j > i; j-- {
			V[i] -= matrix[i][j] * V[j]
		}
	}

	for i := 0; i < nsize; i++ {
		C[i] = V[i]
	}

}

func main() {
	args := os.Args[1:]
	nsize := 1024
	verify := false
	print_matrix := false

	for i, arg := range args {
		if arg == "-v" {
			verify = true
		}
		if arg == "-p" {
			print_matrix = true
		}
		if arg == "-s" {
			s, err := strconv.Atoi(args[i+1])
			if err != nil {
				// handle error
				fmt.Println(err)
				os.Exit(2)
			}
			if s > 0 {
				fmt.Printf("s = %d\n", s)
				nsize = s
			} else {
				fmt.Printf("Entered size is negative, hence using the default (%d)\n", nsize)
			}
		}
	}

	start := time.Now()
	initMatrix(nsize)
	computeGauss(nsize)
	if verify {
		solveGauss(nsize)
	}

	duration := time.Since(start)

	fmt.Print("Application time: ")
	fmt.Println(duration)

	if print_matrix {
		printMatrix(nsize)
		printArray(B, nsize)
		printArray(V, nsize)
		printArray(C, nsize)
	}

	if verify {
		for i := 0; i < nsize; i++ {
			s := fmt.Sprintf("%6.5f %5.5f\n", B[i], C[i])
			io.WriteString(os.Stdout, s)
		}
	}
}
