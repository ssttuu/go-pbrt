

generate:
	go generate ./pkg/pbrt-gen/

build: generate
	go build -o run


run: build
	./run

