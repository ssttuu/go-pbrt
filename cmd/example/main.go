package main

import (
	"context"
	"flag"
	"fmt"
	"log"
	"os"
	"runtime/pprof"
	"time"

	"github.com/stupschwartz/go-pbrt/pkg/accelerator"
	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
)

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to `file`")

func main() {

	flag.Parse()
	if *cpuprofile != "" {
		f, err := os.Create(*cpuprofile)
		if err != nil {
			log.Fatal("could not create CPU profile: ", err)
		}
		if err := pprof.StartCPUProfile(f); err != nil {
			log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}

	// START

	var primitives []pbrt.Primitive

	n := 8

	for k := 0; k <= n; k++ {

		for i := 0; i < 3; i++ {
			var x, y, z float64
			var color pbrt.Spectrum

			switch i {
			case 0:
				x = float64(k) / float64(n) * 100
				color = pbrt.NewRGBSpectrum(1, 0, 0)
			case 1:
				y = float64(k) / float64(n) * 100
				color = pbrt.NewRGBSpectrum(0, 1, 0)
			case 2:
				z = float64(k) / float64(n) * 100
				color = pbrt.NewRGBSpectrum(0, 0, 1)
			}

			radius := 2.0
			sphere := pbrt.NewSphereShape(
				fmt.Sprintf("Sphere: %v, %v, %v - MatteMaterial", x, y, z),
				pbrt.Translate(new(pbrt.Vector3f)), true, radius)

			xform := pbrt.Translate(&pbrt.Vector3f{x, y, z})
			kd := pbrt.NewConstantSpectrumTexture(color)
			sigma := pbrt.NewConstantFloatTexture(0.0)
			geoPrim := pbrt.NewGeometricPrimitive(sphere, pbrt.NewMatteMaterial(kd, sigma, nil))
			xformedPrimitive := pbrt.NewTransformedPrimitive(geoPrim, pbrt.NewAnimatedTransform(xform, xform, 0, 1))

			primitives = append(primitives, xformedPrimitive)
		}
	}

	xformLocal := pbrt.Translate(new(pbrt.Vector3f))
	sphere := pbrt.NewSphereShape("Green Sphere", xformLocal, true, 5.0)

	sigma := pbrt.NewConstantFloatTexture(0.0)
	geoPrim := pbrt.NewGeometricPrimitive(sphere, pbrt.NewMatteMaterial(pbrt.NewConstantSpectrumTexture(pbrt.NewRGBSpectrum(1, 1, 1)), sigma, nil))
	xformPrim := pbrt.NewTransformedPrimitive(geoPrim, pbrt.NewAnimatedTransform(xformLocal, xformLocal, 0, 1))
	primitives = append(primitives, xformPrim)

	agg := accelerator.NewSimpleAggregate(primitives)

	lights := []pbrt.Light{
		pbrt.NewPointLight(
			pbrt.Translate(&pbrt.Vector3f{300, 0, 0}),
			nil,
			pbrt.NewSpectrum(200000),
		),
	}

	scene := pbrt.NewScene(agg, lights, []pbrt.Light{})

	shutterOpen, shutterClose := 0.0, 1.0

	resolution := &pbrt.Point2i{X: 1920, Y: 1080}
	resolution = resolution.DivScalar(2)

	cropBounds := &pbrt.Bounds2f{Min: &pbrt.Point2f{X: 0, Y: 0}, Max: &pbrt.Point2f{X: 1, Y: 1}}
	boxFilter := pbrt.NewBoxFilter(&pbrt.Point2f{X: 1, Y: 1})

	pixelBounds := &pbrt.Bounds2i{Min: &pbrt.Point2i{X: 0, Y: 0}, Max: resolution}
	sampler := pbrt.NewRandomSampler(32, 4)

	film := pbrt.NewFilm(fmt.Sprintf("build/render-%s.png", time.Now().Format(time.RFC3339)), resolution, cropBounds, boxFilter, 100.0, 1.0, 1.0)

	camXform, err := pbrt.LookAt(&pbrt.Point3f{X: 150, Y: 150, Z: 150}, &pbrt.Point3f{X: 0, Y: 0, Z: 0}, &pbrt.Vector3f{0, 1, 0})
	if err != nil {
		log.Fatal(err)
	}
	camXform = camXform.Mul(pbrt.RotateY(-30)).Mul(pbrt.RotateX(-30))

	camAnimXform := pbrt.NewAnimatedTransform(camXform, camXform, 0, 1)
	camera := pbrt.NewPerspectiveCamera(camAnimXform, cropBounds, shutterOpen, shutterClose, 0, 20, 100, film, nil)

	integrator := pbrt.NewDirectLightingIntegrator(pbrt.UniformSampleAll, 10, camera, sampler, pixelBounds)

	ctx := context.Background()
	err = pbrt.Render(ctx, integrator, scene, 16)
	if err != nil {
		log.Fatal(err)
	}

}