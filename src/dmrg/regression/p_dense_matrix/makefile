BOOST_PATH = "/Users/akosenko/workspace/build/boost_1_45_0"
BOOST_LIBS_PATH = "/opt/boost-1.44-binary/lib"
ALPS_PATH = "/Users/akosenko/workspace/build/alps_14.12.2010"
HDF5_PATH = "/Users/akosenko/workspace/build/hdf5_1.8.5"
LOCAL_LIBS_PATH = "/opt/local/lib"
SOURCE_PATH = "/Users/akosenko/workspace/dmrg"
BUILD_PATH = "/Users/akosenko/workspace/build/testing/obj"
OUTPUT_BINARY = "p_dense_matrix.test"

all:
	rm -rf $(BUILD_PATH)/*
	cp -R $(SOURCE_PATH)/* $(BUILD_PATH)
	cd ${BUILD_PATH}; mpic++ -c -g -O0 ambient/packets/packet.cpp        -I.
	cd ${BUILD_PATH}; mpic++ -c -g -O0 ambient/packets/packet_t.cpp      -I.
	cd ${BUILD_PATH}; mpic++ -c -g -O0 ambient/groups/packet_manager.cpp -I.
	cd ${BUILD_PATH}; mpic++ -c -g -O0 ambient/groups/multirank.cpp      -I.
	cd ${BUILD_PATH}; mpic++ -c -g -O0 ambient/groups/group.cpp          -I.
	cd ${BUILD_PATH}; mpic++ -c -g -O0 ambient/ambient.cpp               -I.
	cd ${BUILD_PATH}; mpic++ -I$(SOURCE_PATH) -I$(BOOST_PATH) -I$(ALPS_PATH)/include -I$(HDF5_PATH)/include \
        -L$(HDF5_PATH)/lib -L$(ALPS_PATH)/lib -L$(LOCAL_LIBS_PATH) -L$(BOOST_LIBS_PATH) -lboost -lalps -lhdf5 -larpack -lz -lsz -lboost_unit_test_framework -DBOOST_TEST_DYN_LINK -framework Accelerate \
        -DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DCPU_FREQ=2.5e9 -DHAVE_ALPS_HDF5 -DMPI_PARALLEL \
        -g regression/p_dense_matrix/p_dense_matrix_unit_tests.cpp \
           packet_t.o       \
           packet.o         \
           packet_manager.o \
           multirank.o      \
           group.o          \
           ambient.o        \
        -o $(OUTPUT_BINARY)
	
	cd ${BUILD_PATH}; mpiexec -np 2 ./${OUTPUT_BINARY}
