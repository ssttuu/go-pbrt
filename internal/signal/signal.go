package signal

import (
	"context"
	"github.com/pkg/errors"
	"os"
	"os/signal"
	"syscall"
)

func Func(ctx context.Context) func() error {
	return func() error {
		sigs := make(chan os.Signal, 1)

		signal.Notify(sigs, syscall.SIGINT, syscall.SIGTERM, syscall.SIGKILL)

		select {
		case <-ctx.Done():
			return ctx.Err()
		case sig := <-sigs:
			signal.Reset(syscall.SIGINT, syscall.SIGTERM, syscall.SIGKILL)
			return errors.Errorf("received signal: %v", sig)
		}
	}
}
