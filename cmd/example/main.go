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

	var primitives []*pbrt.TransformedPrimitive

	n := 8

	kd := pbrt.NewConstantSpectrumTexture(pbrt.NewRGBSpectrum(1.0, 0.0, 0.0))
	sigma := pbrt.NewConstantFloatTexture(0.0)

	//for i := 0; i <= n; i++ {
	//	for j := 0; j <= n; j++ {
	for k := 0; k <= n; k++ {
		x := (float64(1) / float64(n) * 200) - 100
		y := (float64(1) / float64(n) * 200) - 100
		z := (float64(k) / float64(n) * 200) - 100
		radius := 2.0

		xformLocal := pbrt.Translate(new(pbrt.Vector3f))
		sphere := pbrt.NewSphereShape(
			fmt.Sprintf("Sphere: %v, %v, %v - MatteMaterial", x, y, z),
			xformLocal, xformLocal.Inverse(), false, radius)

		var geoPrim pbrt.Primitive
		//if (i + j + k % 2) == 0 {
		//	geoPrim = pbrt.NewGeometricPrimitive(sphere, material.NewMirror())
		//} else {
		xform := pbrt.Translate(&pbrt.Vector3f{x, y, z})
		geoPrim = pbrt.NewGeometricPrimitive(sphere, pbrt.NewMatteMaterial(kd, sigma, nil))
		xformedPrimitive := pbrt.NewTransformedPrimitive(geoPrim, pbrt.NewAnimatedTransform(xform, xform, 0, 1))
		//}

		primitives = append(primitives, xformedPrimitive)
	}
	//}
	//}

	xformLocal := pbrt.Translate(new(pbrt.Vector3f))
	sphere := pbrt.NewSphereShape("Green Sphere", xformLocal, xformLocal.Inverse(), false, 5.0)

	xform := pbrt.Translate(&pbrt.Vector3f{10, 0, 0})
	geoPrim := pbrt.NewGeometricPrimitive(sphere, pbrt.NewMatteMaterial(pbrt.NewConstantSpectrumTexture(pbrt.NewRGBSpectrum(0.0, 2.0, 0.0)), sigma, nil))
	xformPrim := pbrt.NewTransformedPrimitive(geoPrim, pbrt.NewAnimatedTransform(xform, xform, 0, 1))
	primitives = append(primitives, xformPrim)

	agg := accelerator.NewSimpleAggregate(primitives)

	lights := []pbrt.Light{
		pbrt.NewPointLight(
			pbrt.Translate(&pbrt.Vector3f{300, 0, 0}),
			nil,
			pbrt.NewSpectrum(10000000),
		),
		pbrt.NewPointLight(
			pbrt.Translate(&pbrt.Vector3f{-300, 0, 0}),
			nil,
			pbrt.NewSpectrum(10000000),
		),
		pbrt.NewPointLight(
			pbrt.Translate(&pbrt.Vector3f{300, 300, 300}),
			nil,
			pbrt.NewSpectrum(10000000),
		),
	}

	scene := pbrt.NewScene(agg, lights, []pbrt.Light{})

	shutterOpen, shutterClose := 0.0, 1.0

	resolution := &pbrt.Point2i{X: 1920, Y: 1080}
	resolution = resolution.DivScalar(1)

	cropBounds := &pbrt.Bounds2f{Min: &pbrt.Point2f{X: 0, Y: 0}, Max: &pbrt.Point2f{X: 1, Y: 1}}
	boxFilter := pbrt.NewBoxFilter(&pbrt.Point2f{X: 1, Y: 1})

	pixelBounds := &pbrt.Bounds2i{Min: &pbrt.Point2i{X: 0, Y: 0}, Max: resolution}
	sampler := pbrt.NewRandomSampler(10, 4)

	film := pbrt.NewFilm(fmt.Sprintf("build/render-%s.png", time.Now().Format(time.RFC3339)), resolution, cropBounds, boxFilter, 100.0, 1.0, 1.0)

	camXform, err := pbrt.LookAt(&pbrt.Point3f{X: 200, Y: 200, Z: 200}, &pbrt.Point3f{X: 0, Y: 0, Z: 0}, &pbrt.Vector3f{0, 0, 1})
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
