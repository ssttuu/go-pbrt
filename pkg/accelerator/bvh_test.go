package accelerator

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func Swap(v []*MortonPrimitive) {
	t := make([]*MortonPrimitive, len(v))
	for i := 0; i < len(v); i++ {
		t[i] = &MortonPrimitive{5, 20}
	}
	v = t
}

func TestRadixSortInPlace(t *testing.T) {
	var mps []*MortonPrimitive
	for i := 0; i < 10; i++ {
		mps = append(mps, &MortonPrimitive{i, uint64(10 - (i + 1))})
	}

	RadixSortInPlace(&mps)

	expected := []*MortonPrimitive{
		{primitiveIndex: 9, mortonCode: 0},
		{primitiveIndex: 8, mortonCode: 1},
		{primitiveIndex: 7, mortonCode: 2},
		{primitiveIndex: 6, mortonCode: 3},
		{primitiveIndex: 5, mortonCode: 4},
		{primitiveIndex: 4, mortonCode: 5},
		{primitiveIndex: 3, mortonCode: 6},
		{primitiveIndex: 2, mortonCode: 7},
		{primitiveIndex: 1, mortonCode: 8},
		{primitiveIndex: 0, mortonCode: 9},
	}
	assert.Equal(t, expected, mps)
}
