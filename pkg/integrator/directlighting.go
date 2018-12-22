package integrator

import (
	"context"
	"log"

	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

type LightStrategy uint8

const (
	UniformSampleAll LightStrategy = iota + 1
	UniformSampleOne
)

type DirectLighting struct {
	camera      pbrt.Camera
	sampler     pbrt.Sampler
	pixelBounds *pbrt.Bounds2i

	strategy      LightStrategy
	maxDepth      int
	nLightSamples []int32
}

func NewDirectLighting(strategy LightStrategy, maxDepth int, camera pbrt.Camera, sampler pbrt.Sampler, pixelBounds *pbrt.Bounds2i) *DirectLighting {
	return &DirectLighting{
		camera:      camera,
		sampler:     sampler,
		pixelBounds: pixelBounds,
		strategy:    strategy,
		maxDepth:    maxDepth,
	}
}

func (dli *DirectLighting) GetSampler() pbrt.Sampler {
	return dli.sampler
}

func (dli *DirectLighting) GetCamera() pbrt.Camera {
	return dli.camera
}

func (dli *DirectLighting) Preprocess(scene pbrt.Scene, sampler pbrt.Sampler) {
	if dli.strategy == UniformSampleAll {
		// compute number of samples to use for each light
		for i := 0; i < len(scene.Lights()); i++ {
			dli.nLightSamples = append(dli.nLightSamples, sampler.RoundCount(scene.Light(i).GetSamples()))
		}

		// request samples for sampling all lights
		for i := 0; i < dli.maxDepth; i++ {
			for j := 0; j < len(scene.Lights()); j++ {
				sampler.Request2DArray(dli.nLightSamples[j])
				sampler.Request2DArray(dli.nLightSamples[j])
			}
		}
	}
}

func (dli *DirectLighting) Li(ctx context.Context, ray *pbrt.Ray, scene pbrt.Scene, sampler pbrt.Sampler, depth int) pbrt.Spectrum {
	L := pbrt.NewSpectrum(0)

	// find closest ray intersection or return background radiance
	si := pbrt.NewSurfaceInteraction()
	intersects := scene.Intersect(ray, si)
	if !intersects {
		for i := 0; i < len(scene.Lights()); i++ {
			L.AddAssign(scene.Light(i).Le(pbrt.NewRayDifferentialFromRay(ray)))
		}
		return L
	}

	// compute scattering function for surface interaction
	si.ComputeScatteringFunctions(ray, false, pbrt.Radiance)
	if si.BSDF == nil {
		return dli.Li(ctx, pbrt.NewRayDifferentialFromRay(si.SpawnRay(ray.Direction)), scene, sampler, depth)
	}

	//compute emitted light if ray hit an area light source
	L.AddAssign(si.Le(si.Wo))

	if len(scene.Lights()) > 0 {
		// compute direct lighting for DirectLighting integrator
		switch dli.strategy {
		case UniformSampleAll:
			L.AddAssign(pbrt.UniformSampleAllLights(si, scene, sampler, dli.nLightSamples, false))
		case UniformSampleOne:
			L.AddAssign(pbrt.UniformSampleOneLight(si, scene, sampler, false, nil))
		default:
			// TODO: don't panic
			log.Panic("unknown lighting strategy")
		}
	}

	if depth+1 < dli.maxDepth {
		// trace rays for specular reflection and refraction
		L.AddAssign(dli.SpecularReflect(ctx, ray, si, scene, sampler, depth+1))
		L.AddAssign(dli.SpecularTransmit(ctx, ray, si, scene, sampler, depth+1))
	}

	return L
}

func (dli *DirectLighting) SpecularReflect(ctx context.Context, ray *pbrt.Ray, si *pbrt.SurfaceInteraction, scene pbrt.Scene, sampler pbrt.Sampler, depth int) pbrt.Spectrum {
	return pbrt.SamplerIntegratorSpecularReflect(dli, ctx, ray, si, scene, sampler, depth)
}

func (dli *DirectLighting) SpecularTransmit(ctx context.Context, ray *pbrt.Ray, si *pbrt.SurfaceInteraction, scene pbrt.Scene, sampler pbrt.Sampler, depth int) pbrt.Spectrum {
	return pbrt.SamplerIntegratorSpecularTransmit(dli, ctx, ray, si, scene, sampler, depth)
}
