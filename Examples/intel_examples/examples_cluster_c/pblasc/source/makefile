#===============================================================================
# Copyright 2010-2020 Intel Corporation.
#
# This software and the related documents are Intel copyrighted  materials,  and
# your use of  them is  governed by the  express license  under which  they were
# provided to you (License).  Unless the License provides otherwise, you may not
# use, modify, copy, publish, distribute,  disclose or transmit this software or
# the related documents without Intel's prior written permission.
#
# This software and the related documents  are provided as  is,  with no express
# or implied  warranties,  other  than those  that are  expressly stated  in the
# License.
#===============================================================================

##  Content:
##      Intel(R) Math Kernel Library PBLAS C examples creation and run
##******************************************************************************

include ../make.inc



# Object files

pblas1_s_obj = obj$(IDE)/pblas1_s_example.o
pblas2_s_obj = obj$(IDE)/pblas2_s_example.o
pblas3_s_obj = obj$(IDE)/pblas3_s_example.o

pblas1_d_obj = obj$(IDE)/pblas1_d_example.o
pblas2_d_obj = obj$(IDE)/pblas2_d_example.o
pblas3_d_obj = obj$(IDE)/pblas3_d_example.o



# Build PBLAS examples

pblas1_s_example : $(pblas1_s_obj)
		mkdir -p $(RES_DIR)
		$(CCLOADER) $(CCLOADERFLAGS) -o $(RES_DIR)/$@ $^ $(LIBS)
pblas2_s_example : $(pblas2_s_obj)
		mkdir -p $(RES_DIR)
		$(CCLOADER) $(CCLOADERFLAGS) -o $(RES_DIR)/$@ $^ $(LIBS)
pblas3_s_example : $(pblas3_s_obj)
		mkdir -p $(RES_DIR)
		$(CCLOADER) $(CCLOADERFLAGS) -o $(RES_DIR)/$@ $^ $(LIBS)

pblas1_d_example : $(pblas1_d_obj)
		mkdir -p $(RES_DIR)
		$(CCLOADER) $(CCLOADERFLAGS) -o $(RES_DIR)/$@ $^ $(LIBS)
pblas2_d_example : $(pblas2_d_obj)
		mkdir -p $(RES_DIR)
		$(CCLOADER) $(CCLOADERFLAGS) -o $(RES_DIR)/$@ $^ $(LIBS)
pblas3_d_example : $(pblas3_d_obj)
		mkdir -p $(RES_DIR)
		$(CCLOADER) $(CCLOADERFLAGS) -o $(RES_DIR)/$@ $^ $(LIBS)



# Run PBLAS examples

pblas1_s_example.txt : pblas1_s_example
		$(MPIRUN) ./$(RES_DIR)/$< 1>$(RES_DIR)/$@ 2>&1
pblas2_s_example.txt : pblas2_s_example
		$(MPIRUN) ./$(RES_DIR)/$< 1>$(RES_DIR)/$@ 2>&1
pblas3_s_example.txt : pblas3_s_example
		$(MPIRUN) ./$(RES_DIR)/$< 1>$(RES_DIR)/$@ 2>&1

pblas1_d_example.txt : pblas1_d_example
		$(MPIRUN) ./$(RES_DIR)/$< 1>$(RES_DIR)/$@ 2>&1
pblas2_d_example.txt : pblas2_d_example
		$(MPIRUN) ./$(RES_DIR)/$< 1>$(RES_DIR)/$@ 2>&1
pblas3_d_example.txt : pblas3_d_example
		$(MPIRUN) ./$(RES_DIR)/$< 1>$(RES_DIR)/$@ 2>&1



# Compile

obj$(IDE)/%.o : %.c
		mkdir -p obj$(IDE)
		$(CC) $(CCFLAGS) -I$(MKLINC) -c -o $@ $<



# Clean all object files
	
clean: ; rm -f -r obj*
