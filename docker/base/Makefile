DOCKER = docker
REPO = gcre/$(notdir $(CURDIR))	

info:
	echo $(REPO)

run: build
	$(DOCKER) run -it $(REPO)

build:
	$(DOCKER) build --no-cache --force-rm -t $(REPO) .

push: build
	$(DOCKER) push $(REPO)
