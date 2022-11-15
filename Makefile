EXECUTABLE_BINARY_PATH := "./target/release/vis-a-vis"
EXECUTABLE_BINARY_VERSION := $(shell $(EXECUTABLE_BINARY_PATH) --version)
DEPLOY_DEST_DIR_NAME := "visavis-${EXECUTABLE_BINARY_VERSION}"
DEPLOY_DEST_DIR_PATH := "deploy/${DEPLOY_DEST_DIR_NAME}"
DEPLOY_DEST_ARCHIVE_NAME := "visavis-${EXECUTABLE_BINARY_VERSION}.tar.bz2"


DEFAULT_PARAMETERS := "parameters/default.json"
DEFAULT_PROTOCOL := "protocols/default.protocol"

ifeq ($(PARAMETERS),)
PARAMETERS=${DEFAULT_PARAMETERS}
endif

ifeq ($(PROTOCOL),)
PROTOCOL=${DEFAULT_PROTOCOL}
endif


all: compile


compile:
	@cargo build --release


run:
	@bash -c "if [ x${PARAMETERS} == x${DEFAULT_PARAMETERS} ]; then echo Note: Assuming default parameter values from ${PARAMETERS}.; fi"
	@bash -c "if [ x${PROTOCOL}   == x${DEFAULT_PROTOCOL}   ]; then echo Note: Assuming default protocol from ${PROTOCOL}.; fi"
	@cargo run --release $(PARAMETERS) $(PROTOCOL) -- --images


archive:
	@bash -c 'if [ -d visavis ]; then rm -fR deploy/visavis; fi'

	@mkdir -p                  deploy/visavis
	@cp -a  ReadMe.md          deploy/visavis/
	@cp -a  License.txt        deploy/visavis/
	@cp -ar protocols          deploy/visavis/

	@mkdir                     deploy/visavis/parameters
	@cp -ar parameters/*.json  deploy/visavis/parameters

	@cp -a  Cargo.toml         deploy/visavis/
	@mkdir                     deploy/visavis/src
	@cp -a  src/*.rs           deploy/visavis/src/

	@mkdir                         deploy/visavis/extra
	@cp -a  extra/*.py             deploy/visavis/extra/
	@mkdir                         deploy/visavis/extra/parameters
	@cp -a  extra/parameters/*.py  deploy/visavis/extra/parameters/

	@rm -fR $(DEPLOY_DEST_DIR_PATH)
	@mv deploy/visavis $(DEPLOY_DEST_DIR_PATH)
	@( cd deploy && \
	   tar cfj $(DEPLOY_DEST_ARCHIVE_NAME) $(DEPLOY_DEST_DIR_NAME) ) && \
	 echo "Created archive deploy/${DEPLOY_DEST_ARCHIVE_NAME} (source code with parameter sets and protocols)."


.PHONY:      \
	all      \
	compile  \
	run      \
	archive
