package pbrt

import (
	"testing"
	"github.com/stretchr/testify/assert"
)

func TestMatchesFlags(t *testing.T) {
	assert.True(t, MatchesFlags(BSDFDiffuse, BSDFDiffuse))
	assert.True(t, MatchesFlags(BSDFDiffuse, BSDFDiffuse|BSDFReflection))
	assert.True(t, MatchesFlags(BSDFReflection, BSDFDiffuse|BSDFReflection))
	assert.False(t, MatchesFlags(BSDFReflection, BSDFDiffuse))
	assert.True(t, MatchesFlags(BSDFReflection, BSDFAll))
}
