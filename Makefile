
all:
	env RDBASE=$(RDBASE) ninja -C build install

test: RDBASE:=$(shell pwd)
test:
	cd build; env RDBASE=$(RDBASE) ctest
