TARGET_EXEC ?= perflab
TARGET_PLAIN ?= plain

TARGET_DIR ?= ./
BUILD_DIR ?= ./build
SRC_DIRS ?= ./
OPT ?= -O3 -ftree-vectorize -w
GCC ?= gcc

SERVER = lnxsrv09
PORT = 15213
SUBMIT_API = http://${SERVER}.seas.ucla.edu:${PORT}/upload
CHECK_QUEUE_API = http://${SERVER}.seas.ucla.edu:${PORT}/ping-submission-queue
FILE = $(shell basename "$$PWD").tar.gz
TARGETS = kernels.c cookie.txt comment.txt

SRCS := $(shell find $(SRC_DIRS) -maxdepth 1 -name "*.cpp" -or -name "*.c")
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)

DEPS := $(OBJS:.o=.d)

CPPFLAGS ?= -I. -MMD -MP -g $(OPT) -Wall -std=c11 -D_GNU_SOURCE

LDFLAGS ?= -lm checker_lib.a

TOBJS := $(filter-out $(BUILD_DIR)/./plain.c.o,$(OBJS))

all: $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/$(TARGET_PLAIN) 

$(TARGET_DIR)/$(TARGET_EXEC): $(TOBJS) 
	$(GCC) $(TOBJS) -o $@ $(LDFLAGS) 

$(TARGET_DIR)/$(TARGET_PLAIN): $(BUILD_DIR)/plain.c.o $(BUILD_DIR)/kernels.c.o
	$(GCC) $^ -o $@ $(LDFLAGS)

# Add this option "-fopt-info-vec-all" after "$(GCC) to see the compiler's vectorization report!! It's so much "fun".
$(BUILD_DIR)/./kernels.c.o: kernels.c
	$(MKDIR_P) $(dir $@)
	$(GCC) $(CPPFLAGS) -mavx  $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.c.o: %.c
	$(MKDIR_P) $(dir $@)
	$(GCC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

.PHONY: clean all

submit: create-tar
		echo 'sending request to the server...'
		echo "it might take a few seconds..."
		curl -F 'file=@$(FILE)' $(SUBMIT_API)

ping:
		curl $(CHECK_QUEUE_API)
create-tar:
		@tar -czf $(FILE) $(TARGETS) 2>/dev/null

clean:
	$(RM) -r $(BUILD_DIR) $(TARGET_DIR)/$(TARGET_EXEC) $(TARGET_DIR)/$(TARGET_PLAIN)

-include $(DEPS)

MKDIR_P ?= mkdir -p
.SILENT:
