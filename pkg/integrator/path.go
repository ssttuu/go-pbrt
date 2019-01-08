package integrator

import (
	"context"

	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
)

func NewPath(maxDepth int32, camera pbrt.Camera, sampler pbrt.Sampler, pixelBounds *pbrt.Bounds2i, rrThreshold float64, lightSampleStrategy pbrt.LightSampleStrategy) pbrt.Integrator {
	return &Path{
		SamplerIntegrator:   pbrt.NewSamplerIntegrator(camera, sampler, pixelBounds),
		maxDepth:            maxDepth,
		rrThreshold:         rrThreshold,
		lightSampleStrategy: lightSampleStrategy,
	}
}

type Path struct {
	*pbrt.SamplerIntegrator

	maxDepth            int32
	rrThreshold         float64
	lightSampleStrategy pbrt.LightSampleStrategy
	lightDistribution   pbrt.LightDistribution
}

func (p *Path) Preprocess(scene pbrt.Scene, sampler pbrt.Sampler) {
	p.lightDistribution = pbrt.CreateLightSampleDistribution(p.lightSampleStrategy, scene)
}

func (p *Path) Li(ctx context.Context, r *pbrt.Ray, scene pbrt.Scene, sampler pbrt.Sampler, depth int) pbrt.Spectrum {
	L := pbrt.NewSpectrum(0.0)
	beta := pbrt.NewSpectrum(1.0)
	ray := pbrt.NewRayDifferentialFromRay(r)
	specularBounce := false
	bounces := int32(0)
	etaScale := 1.0

	for {
		bounces++

		// intersect ray with scene and store intersection in isect
		isect := pbrt.NewSurfaceInteraction()
		foundIntersection := scene.Intersect(ray, isect)

		// possibly add emitted light at intersection
		if bounces == 0 || specularBounce {
			// add emitted light at path vertex or from the environment
			if foundIntersection {
				Le := isect.Le(ray.Direction.MulScalar(-1))
				betaLe := beta.Mul(Le)
				L.AddAssign(betaLe)
				p.check(L)
			} else {
				for _, light := range scene.InfiniteLights() {
					Le := light.Le(ray)
					betaMul := beta.Mul(Le)
					L.AddAssign(betaMul)
					p.check(L)
				}
			}
		}

		// terminate path if ray escaped or maxDepth was reached
		if !foundIntersection || bounces >= p.maxDepth {
			p.check(L)
			break
		}

		// Compute scattering functions and skip over medium boundaries
		isect.ComputeScatteringFunctions(ray, true, pbrt.Radiance)
		if isect.BSDF == nil {
			// TODO: double check
			ray = isect.SpawnRay(ray.Direction)
			bounces--
			continue
		}

		distrib := p.lightDistribution.Lookup(isect.Point)

		//sample illumination from lights to find path contribution
		// - skip this for perfectly specular BSDFs
		if isect.BSDF.NumComponents(pbrt.BSDFAll&^pbrt.BSDFSpecular) > 0 {
			Ld := beta.Mul(pbrt.UniformSampleOneLight(isect, scene, sampler, false, distrib))
			L.AddAssign(Ld)
			p.check(L)
		}

		// sample BSDF to get new path direction
		wo := ray.Direction
		f, wi, pdf, flags := isect.BSDF.SampleF(wo, sampler.Get2D(), pbrt.BSDFAll)
		if f.IsBlack() || pdf == 0.0 {
			p.check(L)
			break
		}
		wiAbsDot := wi.AbsDot(isect.Shading.Normal)
		wiAbsDotPdf := wiAbsDot / pdf
		fMul := f.MulScalar(wiAbsDotPdf)
		beta.MulAssign(fMul)
		p.check(beta)

		specularBounce = flags&pbrt.BSDFSpecular != 0

		if flags&pbrt.BSDFSpecular > 0 && flags&pbrt.BSDFTransmission > 0 {
			eta := isect.BSDF.Eta
			// Update the term that tracks radiance scaling for refraction
			// depending on whether the ray is entering or leaving the
			// medium.
			if wo.Dot(isect.Normal) > 0 {
				etaScale *= eta * eta
			} else {
				etaScale *= 1 / (eta * eta)
			}
		}

		ray = isect.SpawnRay(wi)

		// Account for subsurface scattering, if applicable
		if isect.BSSRDF != nil && flags&pbrt.BSDFTransmission > 0 {
			// importance sample the BSSRDF
			pi := pbrt.NewSurfaceInteraction()
			S, pdf := isect.BSSRDF.SampleS(scene, sampler.Get1D(), sampler.Get2D(), pi)
			if S.IsBlack() || pdf == 0 {
				break
			}
			beta = S.DivScalar(pdf)

			// account for the direct subsurface scattering component
			L.AddAssign(beta.Mul(pbrt.UniformSampleOneLight(pi, scene, sampler, false, p.lightDistribution.Lookup(pi.Point))))

			// account for the indirect subsurface scattering component
			f, wi, pdf, _ := pi.BSDF.SampleF(pi.Wo, sampler.Get2D(), pbrt.BSDFAll)
			if f.IsBlack() || pdf == 0 {
				break
			}
			beta.MulAssign(f.MulScalar(wi.AbsDot(pi.Shading.Normal) / pdf))
			specularBounce = flags & pbrt.BSDFSpecular != 0
			ray = pi.SpawnRay(wi)

		}

		// possibly terminate the path with Russian roulette.
		// factor out radiance scaling due to refraction in rrBeta
		rrBeta := beta.MulScalar(etaScale)
		if rrBeta.MaxComponentValue() < p.rrThreshold && bounces > 3 {
			q := math.Max(0.05, 1-rrBeta.MaxComponentValue())
			if sampler.Get1D() < q {
				p.check(L)
				break
			}
			beta = beta.DivScalar(1 - q)
		}
	}

	return L
}

func (p *Path) check(L pbrt.Spectrum) {
	//if L.MaxComponentValue() > 1 {
	//	fmt.Println(L)
	//}
}

func (p *Path) SpecularReflect(ctx context.Context, ray *pbrt.Ray, si *pbrt.SurfaceInteraction, scene pbrt.Scene, sampler pbrt.Sampler, depth int) pbrt.Spectrum {
	return pbrt.SamplerIntegratorSpecularReflect(p, ctx, ray, si, scene, sampler, depth)
}
func (p *Path) SpecularTransmit(ctx context.Context, ray *pbrt.Ray, si *pbrt.SurfaceInteraction, scene pbrt.Scene, sampler pbrt.Sampler, depth int) pbrt.Spectrum {
	return pbrt.SamplerIntegratorSpecularTransmit(p, ctx, ray, si, scene, sampler, depth)
}
