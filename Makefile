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
	rm -f $(NAME) $(OBJS) $(NAME).gcno $(NAME).gcda $(NAME)out.*
	rm -rf Doxyfile html

Doxyfile: Makefile
	doxygen -g - | sed \
	  -e '/PROJECT_NAME/s/".*"/"$(NAME)"/' \
	  -e '/^HAVE_DOT/s/NO/YES/' \
	  -e '/^CALL.*_GRAPH/s/NO/YES/' \
	  -e '/^EXTRACT_/s/ NO/ YES/' \
	  -e '/^INLINE_SOURCES/s/NO/YES/' \
	  -e '/^SOURCE_BROWSER/s/NO/YES/' \
	  -e '/^GENERATE_LATEX/s/YES/NO/' \
	   > $@

.PHONY: doc
doc: html

html: $(NAME).c $(EXTRA_SRC) Doxyfile
	doxygen Doxyfile
