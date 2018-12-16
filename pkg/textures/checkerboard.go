package textures

import (
	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type AntiAliasingMethod int

const (
	none AntiAliasingMethod = iota
	closedForm
)

func NewCheckerboard2D(m pbrt.TextureMapping2D, tex1, tex2 pbrt.SpectrumTexture) *Checkerboard2D {
	return &Checkerboard2D{
		mapping: m,
		tex1: tex1,
		tex2: tex2,
		antiAliasingMethod: none,
	}
}

type Checkerboard2D struct {
	mapping            pbrt.TextureMapping2D
	tex1, tex2         pbrt.SpectrumTexture
	antiAliasingMethod AntiAliasingMethod
}

func (c *Checkerboard2D) Evaluate(si *pbrt.SurfaceInteraction) pbrt.Spectrum {
	st, _, _ := c.mapping.Map(si)
	if c.antiAliasingMethod == none {
		if int(math.Floor(st.X) + math.Floor(st.Y)) % 2 == 0 {
			return c.tex1.Evaluate(si)
		}
		return c.tex2.Evaluate(si)
	}
	// TODO:
	return pbrt.NewSpectrum(1.0)
}

