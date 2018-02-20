package geometry

import "testing"

func TestVector2_Add(t *testing.T) {
	v1 := VectorXY{1, 2}
	v2 := VectorXY{2, 3}

	result := v1.Add(&v2)
	if !result.Equals(&VectorXY{3,5}) {
		t.Fail()
	}

	resultV3 := v1.Add(&VectorXYZ{4, 5, 6})
	if !resultV3.Equals(&VectorXY{5, 7}) {
		t.Fail()
	}
}
