package pbrt

import (
	"testing"
	"fmt"
)

func Swap(v []*MortonPrimitive) {
	t := make([]*MortonPrimitive, len(v))
	for i := 0; i < len(v); i++ {
		t[i] = &MortonPrimitive{5, 20}
	}
	v = t
}

func PrintIt(v []*MortonPrimitive) {
	for _, mp := range v {
		fmt.Printf("%+v, ", mp)
	}
	fmt.Print("\n")
}

func TestRadixSortInPlace(t *testing.T) {
	var mps []*MortonPrimitive
	for i := 0; i < 10; i++ {
		mps = append(mps, &MortonPrimitive{primitiveIndex:i, mortonCode:uint64(10-i)})
	}

	PrintIt(mps)
	RadixSortInPlace(&mps)
	PrintIt(mps)
	t.Fail()
}