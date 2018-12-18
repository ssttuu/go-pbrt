package render

import (
	"context"
	"fmt"
	"github.com/pkg/errors"
	"github.com/ssttuu/go-pbrt/pkg/accelerator"
	"github.com/ssttuu/go-pbrt/pkg/integrator"
	"github.com/ssttuu/go-pbrt/pkg/lights"
	"github.com/ssttuu/go-pbrt/pkg/materials"
	"github.com/ssttuu/go-pbrt/pkg/math"
	"github.com/ssttuu/go-pbrt/pkg/pbrt"
	"github.com/ssttuu/go-pbrt/pkg/proto/render"
	"github.com/ssttuu/go-pbrt/pkg/shapes"
	"github.com/ssttuu/go-pbrt/pkg/textures"
	"google.golang.org/grpc"
	"os"
	"time"
)

func RegisterServer(s *grpc.Server) error {
	render.RegisterRenderServer(s, &server{})
	return nil
}

type server struct{}

func (s *server) Render(ctx context.Context, req *render.RenderRequest) (*render.RenderResponse, error) {
	var primitives []pbrt.Primitive

	n := 8

	for k := 1; k < n; k++ {
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
			y = math.Max(y, radius/2)

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

	xformLocal := pbrt.Translate(&pbrt.Vector3f{50, 2.5, 50})
	sphere := pbrt.NewSphereShape("Green Sphere", pbrt.Translate(new(pbrt.Vector3f)), true, 5.0)

	glass := materials.NewGlass(
		pbrt.NewConstantSpectrumTexture(pbrt.NewSpectrum(1)),
		pbrt.NewConstantSpectrumTexture(pbrt.NewSpectrum(1)),
		pbrt.NewConstantFloatTexture(0),
		pbrt.NewConstantFloatTexture(0),
		pbrt.NewConstantFloatTexture(1.5),
		nil,
	)

	geoPrim := pbrt.NewGeometricPrimitive(sphere, glass)
	xformPrim := pbrt.NewTransformedPrimitive(geoPrim, pbrt.NewAnimatedTransform(xformLocal, xformLocal, 0, 1))
	primitives = append(primitives, xformPrim)

	sigma := pbrt.NewConstantFloatTexture(0.0)
	diskXform := pbrt.Translate(new(pbrt.Vector3f)).Mul(pbrt.RotateX(90))
	disk := shapes.NewDisk(diskXform, 0.01, 10000, 0, 360)
	checkerboard := textures.NewCheckerboard2D(
		pbrt.NewPlanarMapping2D(&pbrt.Vector3f{.2, 0, 0}, &pbrt.Vector3f{0, 0, .2}, 0, 0),
		pbrt.NewConstantSpectrumTexture(pbrt.NewSpectrum(1.0)),
		pbrt.NewConstantSpectrumTexture(pbrt.NewSpectrum(0.18)),
	)
	diskPrim := pbrt.NewGeometricPrimitive(disk, pbrt.NewMatteMaterial(checkerboard, sigma, nil))
	primitives = append(primitives, diskPrim)

	agg := accelerator.NewBVH(primitives, 2, accelerator.SplitSAH)

	ls := []pbrt.Light{
		lights.NewDistant(
			pbrt.Translate(&pbrt.Vector3f{-100, 100, 100}),
			pbrt.NewSpectrum(0.05),
			&pbrt.Vector3f{-1, 1, 1},
		),
		lights.NewPoint(
			pbrt.Translate(&pbrt.Vector3f{50, 20, 50}),
			nil,
			pbrt.NewSpectrum(100),
		),
		lights.NewPoint(
			pbrt.Translate(&pbrt.Vector3f{-50, 30, -50}),
			nil,
			pbrt.NewSpectrum(50),
		),
	}

	scene := pbrt.NewScene(agg, ls)

	shutterOpen, shutterClose := 0.0, 1.0

	resolution := &pbrt.Point2i{X: int64(req.Width), Y: int64(req.Height)}

	cropBounds := &pbrt.Bounds2f{Min: &pbrt.Point2f{X: 0, Y: 0}, Max: &pbrt.Point2f{X: 1, Y: 1}}
	boxFilter := pbrt.NewBoxFilter(&pbrt.Point2f{X: 1, Y: 1})

	pixelBounds := &pbrt.Bounds2i{Min: &pbrt.Point2i{X: 0, Y: 0}, Max: resolution}
	sampler := pbrt.NewRandomSampler(8, 4)

	if _, err := os.Stat("build"); os.IsNotExist(err) {
		if err := os.Mkdir("build", 0644); err != nil {
			return nil, errors.Wrap(err, "making directory")
		}
	}
	path := fmt.Sprintf("build/render-%s.png", time.Now().Format(time.RFC3339))
	film := pbrt.NewFilm(path, resolution, cropBounds, boxFilter, 100.0, 1.0, 1.0)

	camXform, err := pbrt.LookAt(&pbrt.Point3f{X: 150, Y: 150, Z: 150}, &pbrt.Point3f{X: 0, Y: 0, Z: 0}, &pbrt.Vector3f{0, 1, 0})
	if err != nil {
		return nil, errors.Wrap(err, "look at")
	}
	camXform = camXform.Mul(pbrt.RotateY(-30)).Mul(pbrt.RotateX(-30))

	camAnimXform := pbrt.NewAnimatedTransform(camXform, camXform, 0, 1)
	camera := pbrt.NewPerspectiveCamera(camAnimXform, cropBounds, shutterOpen, shutterClose, 0, 20, 100, film, nil)

	dli := integrator.NewDirectLighting(integrator.UniformSampleAll, 10, camera, sampler, pixelBounds)

	err = pbrt.Render(ctx, dli, scene, 16)
	if err != nil {
		return nil, errors.Wrap(err, "rendering")
	}

	return &render.RenderResponse{
		Path: path,
	}, nil
}
