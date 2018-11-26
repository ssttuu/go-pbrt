package efloat

import (
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestEFloat_Add(t *testing.T) {
	f := New(1.0, 0.0)
	ret := f.Add(New(1.0, 0.0))
	assert.Equal(t, &EFloat{Value: 2.0, Low: 1.9999999999999998, High: 2.0000000000000004}, ret)
}
