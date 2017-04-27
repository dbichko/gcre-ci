DOCKER = docker
REPO = gcre/$(notdir $(CURDIR))

build: .built

push: build
	$(DOCKER) push $(REPO)

.built: Dockerfile
	$(DOCKER) build -t $(REPO) .
	touch .built

run:
	$(DOCKER) run -it $(REPO)

clean:
	rm -f .built
