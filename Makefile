NAME=mdlj

EXTRA_SRC=

LIBS=-lm

LDFLAGS=-g

CFLAGS=-Wall -g

OBJS=$(patsubst %.c,%.o,$(NAME).c $(EXTRA_SRC))

$(NAME): $(OBJS)
	${CC} ${LDFLAGS} -o ${NAME} ${OBJS} ${LIBS}

.PHONY: clean
clean:
	rm -f $(NAME) $(OBJS)
