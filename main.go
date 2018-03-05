package main

import (
	"github.com/stupschwartz/go-pbrt/core"
	"math"
	"image"
	"image/color"
	"os"
	"image/png"
	"log"
	"fmt"

	"github.com/uber/jaeger-client-go/config"
	"github.com/uber/jaeger-lib/metrics/prometheus"
	"github.com/opentracing/opentracing-go"
	"io"
	"flag"
	"runtime/pprof"
)

type TraceData struct {
	x, y  int
	depth int
}

type TraceResult struct {
	x, y  int
	value *pbrt.Vector3
}

const INFINITY float64 = 10e8

func trace(origin, direction *pbrt.Vector3, spheres []*pbrt.Sphere, depth int) *pbrt.Vector3 {
	tNear := INFINITY

	//var sphere pbrt.Sphere
	foundSphere := false

	// Find sphere intersection
	var t0, t1 float64
	var intersects bool

	for i := 0; i < len(spheres); i++ {
		t0, t1, intersects = spheres[i].Intersects(origin, direction)
		if (intersects) {
			if t0 < 0 {
				t0 = t1
			}
			if t0 < tNear {
				tNear = t0
				//sphere := spheres[i]
				foundSphere = true
			}
		}
	}

	if !foundSphere {
		return pbrt.NewVector3(0.1, 0.1, 0.1)
	}

	return pbrt.NewVector3(1, 1, 1)
}

func worker(parentSpan opentracing.Span, id int, camera CameraSettings, spheres []*pbrt.Sphere, traceQueue <-chan *TraceData, traceResults chan <- *TraceResult, done chan <- bool) {
	span := parentSpan.Tracer().StartSpan("worker",
		opentracing.ChildOf(parentSpan.Context()),
	)
	defer span.Finish()

	//previousSpanContext := span.Context()

	fmt.Println("Worker", id, "started")
	itemProcessing := 0

	origin := pbrt.NewVector3(0, 0, 0)
	rayDir := pbrt.NewVector3(0, 0, -1)

	var xx, yy float64

	for traceItem := range traceQueue {
		xx = (2.0 * ((float64(traceItem.x) + 0.5) * camera.InverseWidth) - 1.0) * camera.Angle * camera.AspectRatio
		yy = (1.0 - 2.0 * ((float64(traceItem.y) + 0.5) * camera.InverseHeight)) * camera.Angle

		rayDir.SetX(xx)
		rayDir.SetY(yy)
		rayDir.Normalize()

		pixel := trace(origin, rayDir, spheres, traceItem.depth)
		traceResults <- &TraceResult{
			x: traceItem.x,
			y: traceItem.y,
			value: pixel,
		}
		itemProcessing += 1
	}
	done <- true
}

type CameraSettings struct {
	Width         int
	Height        int
	InverseWidth  float64
	InverseHeight float64
	FOV           float64
	AspectRatio   float64
	Angle         float64
}

func NewCameraSettings(width, height int) CameraSettings {
	fov := 30.0
	return CameraSettings{
		Width: width,
		Height: height,
		InverseWidth: 1 / float64(width),
		InverseHeight: 1 / float64(width),
		FOV: fov,
		AspectRatio: float64(width) / float64(height),
		Angle: math.Tan(math.Pi * 0.5 * fov / 180.0),

	}
}

func render(parentSpan opentracing.Span, spheres []*pbrt.Sphere) {
	renderSpan := parentSpan.Tracer().StartSpan("render",
		opentracing.ChildOf(parentSpan.Context()),
	)
	defer renderSpan.Finish()

	camera := NewCameraSettings(1920, 1080)

	imgPng := image.NewNRGBA(image.Rect(0, 0, camera.Width, camera.Height))

	traceQueue := make(chan *TraceData, 2048)
	traceResults := make(chan *TraceResult, 2048)
	nWorkers := 8
	done := make(chan bool, nWorkers)
	for i := 0; i < nWorkers; i++ {
		go worker(renderSpan, i, camera, spheres, traceQueue, traceResults, done)
	}

	go func() {
		populateTraceQueue := parentSpan.Tracer().StartSpan("populate-trace",
			opentracing.ChildOf(renderSpan.Context()),
		)
		for y := 0; y < camera.Height; y++ {
			for x := 0; x < camera.Width; x++ {

				traceQueue <- &TraceData{
					x: x,
					y: y,
					depth: 0,
				}
			}
		}
		close(traceQueue)
		populateTraceQueue.Finish()
	}()

	writeImageTrace := parentSpan.Tracer().StartSpan("write-image",
		opentracing.ChildOf(renderSpan.Context()),
	)
	nWorkersCompleted := 0
	for {
		if nWorkersCompleted == nWorkers {
			break
		}

		select {
		case traceResult := <-traceResults:
			imgPng.Set(traceResult.x, traceResult.y, color.NRGBA{
				R: uint8(traceResult.value.GetX() * 255),
				G: uint8(traceResult.value.GetY() * 255),
				B: uint8(traceResult.value.GetZ() * 255),
				A: 255,
			})
		case <-done:
			nWorkersCompleted += 1
		}
	}

	close(traceResults)

	f, err := os.Create("image.png")
	if err != nil {
		log.Fatal(err)
	}

	if err := png.Encode(f, imgPng); err != nil {
		f.Close()
		log.Fatal(err)
	}

	if err := f.Close(); err != nil {
		log.Fatal(err)
	}

	writeImageTrace.Finish()

}

func getTracer() (opentracing.Tracer, io.Closer) {
	metricsFactory := prometheus.New()
	cfg := config.Configuration{
		Sampler: &config.SamplerConfig{
			Type: "const",
			Param: 1,
		},
		Reporter: &config.ReporterConfig{
			LocalAgentHostPort: "jaeger:6831",
		},
	}
	tracer, closer, err := cfg.New(
		"pbrt",
		config.Metrics(metricsFactory),
	)
	if err != nil {
		log.Fatal("Failed to start tracer", err)
	}

	return tracer, closer

}

var cpuprofile = flag.String("cpuprofile", "", "write cpu profile to file")

func main() {
	flag.Parse()
    if *cpuprofile != "" {
        f, err := os.Create(*cpuprofile)
        if err != nil {
            log.Fatal(err)
        }
        pprof.StartCPUProfile(f)
        defer pprof.StopCPUProfile()
    }

	//tracer, closer := getTracer()
	//defer closer.Close()

	// START

	spheres := []*pbrt.Sphere{}

	n := 9

	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				x := (float64(i) / float64(n) * 200) - 100
				y := (float64(j) / float64(n) * 200) - 100
				z := (float64(k) / float64(n) * 200) - 100
				radius := 1.0
				if x > 0 {
					radius = 2.0
				}
				if z < 0 {
					radius *= 2
				}
				spheres = append(spheres, &pbrt.Sphere{
					Center: pbrt.NewVector3(x, y, z),
					Radius: radius,
					SurfaceColor: pbrt.NewVector3(1.0, 1.0, 1.0),
					Reflection: 0,
					Transparency: 0,
				}, )
			}
		}
	}

	span := opentracing.StartSpan("render")
	defer span.Finish()

	render(span, spheres)
}
