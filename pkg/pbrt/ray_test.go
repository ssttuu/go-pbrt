package pbrt

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stupschwartz/go-pbrt/pkg/math"
)

func TestOffsetRayOrigin(t *testing.T) {
	out := OffsetRayOrigin(
		&Point3f{},
		&Vector3f{math.MachineEpsilon, math.MachineEpsilon, math.MachineEpsilon},
		&Normal3f{1, 1, 1},
		&Vector3f{1, 1, 1},
	)

	assert.Equal(t, &Point3f{1.5183e-320, 1.5183e-320, 1.5183e-320}, out)
}
