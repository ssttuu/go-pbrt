package pbrt

import (
	"context"
	"fmt"
	"sync/atomic"
	"time"
)

type ProgressReporter interface {
	Start(ctx context.Context)
	Step()
}

func NewProgress(steps int64) *StdoutProgress {
	p := &StdoutProgress{
		total: steps,
		step: 0,
		steps: make(chan int64, steps),
	}
	return p
}

type StdoutProgress struct {
	total int64
	step int64
	steps chan int64

	startTime, endTime time.Time
}

func (p *StdoutProgress) Start(ctx context.Context) {
	p.startTime = time.Now()
	fmt.Printf("Started: %s\n", p.startTime.Format(time.RFC3339))

	go func() {
		defer func() {
			p.endTime = time.Now()
			fmt.Printf("\nCompleted: %s\n", p.endTime.Format(time.RFC3339))
			fmt.Printf("Duration: %s\n", p.endTime.Sub(p.startTime))
		}()

		defer close(p.steps)
		for step := range p.steps {
			select {
			case <-ctx.Done():
				return
			default:
				fmt.Printf("\rProgress: %3.2f%%", float64(step)/float64(p.total) * 100)
				if step == p.total {
					return
				}
			}
		}
	}()
}

func (p *StdoutProgress) Step() {
	atomic.AddInt64(&p.step, 1)
	p.steps <- p.step
}
