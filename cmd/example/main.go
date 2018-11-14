package main

import (
	"github.com/stupschwartz/go-pbrt/pkg/pbrt"
	"context"
	"log"
	"flag"
	"os"
	"runtime/pprof"
	"fmt"
	"time"
	"github.com/stupschwartz/go-pbrt/pkg/accelerator"
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

	kd := pbrt.NewConstantSpectrumTexture(pbrt.NewRGBSpectrum(1.0, 0.0, 0.0))
	sigma := pbrt.NewConstantFloatTexture(0.0)

	for i := 0; i <= n; i++ {
		for j := 0; j <= n; j++ {
			for k := 0; k <= n; k++ {
				x := (float64(i) / float64(n) * 200) - 100
				y := (float64(j) / float64(n) * 200) - 100
				z := (float64(k) / float64(n) * 200) - 100
				radius := 2.0

				xform := pbrt.Translate(&pbrt.Vector3f{x, y, z})
				sphere := pbrt.NewSphereShape(
					fmt.Sprintf("Sphere: %v, %v, %v - MatteMaterial", x, y, z),
					xform, xform.Inverse(), false, radius)

				var geoPrim pbrt.Primitive
				//if (i + j + k % 2) == 0 {
				//	geoPrim = pbrt.NewGeometricPrimitive(sphere, material.NewMirror())
				//} else {
				geoPrim = pbrt.NewGeometricPrimitive(sphere, pbrt.NewMatteMaterial(kd, sigma, nil))
				//}

				primitives = append(primitives, geoPrim)
			}
		}
	}


	xform := pbrt.Translate(&pbrt.Vector3f{10, 0, 0})
	sphere := pbrt.NewSphereShape("Green Sphere", xform, xform.Inverse(), false, 5.0)

	primitives = append(primitives, pbrt.NewGeometricPrimitive(sphere, pbrt.NewMatteMaterial(pbrt.NewConstantSpectrumTexture(pbrt.NewRGBSpectrum(0.0, 2.0, 0.0)), sigma, nil)))

	agg := accelerator.NewSimpleAggregate(primitives)

	//areaLightXform := pbrt.Translate(&pbrt.Vector3f{0, 0, 0})
	//areaLightShape := pbrt.NewSphereShape(areaLightXform, areaLightXform.Inverse(), false, 10.0)
	//areaLight := pbrt.NewDiffuseAreaLight(areaLightXform, nil, pbrt.NewSpectrum(100000.0), 10, areaLightShape, false)

	pointLight := pbrt.NewPointLight(
		pbrt.Translate(&pbrt.Vector3f{0, 0, 0}),
		nil,
		pbrt.NewSpectrum(pbrt.Pi*100000),
	)

	lights := []pbrt.Light{
		pointLight,
		//areaLight,
	}

	scene := &pbrt.Scene{
		Aggregate:  agg,
		Lights:     lights,
		WorldBound: agg.WorldBound(),
	}

	shutterOpen, shutterClose := 0.0, 1.0

	resolution := &pbrt.Point2i{X: 1920, Y: 1080}

	cropBounds := &pbrt.Bounds2f{Min: &pbrt.Point2f{X: 0, Y: 0}, Max: &pbrt.Point2f{X: 1, Y: 1}}
	boxFilter := pbrt.NewBoxFilter(&pbrt.Point2f{X: 1,Y: 1})

	pixelBounds := &pbrt.Bounds2i{Min: &pbrt.Point2i{X:0, Y: 0}, Max: resolution}
	sampler := pbrt.NewRandomSampler(10, 4)

	film := pbrt.NewFilm(fmt.Sprintf("build/render-%s.png", time.Now().Format(time.RFC3339)), resolution, cropBounds, boxFilter, 100.0, 1.0, 1.0)

	camXform, err := pbrt.LookAt(&pbrt.Point3f{X: 200, Y: 200, Z: 200}, &pbrt.Point3f{X: 0,Y: 0, Z: 0}, &pbrt.Vector3f{0, 0, 1})
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
