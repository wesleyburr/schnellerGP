/* ------------------------------------------------------------
 This is the file "parameters.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2009
 ------------------------------------------------------------ */

#include "parameters.h"

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int
askforint(const char *question, const char *envname, int deflt)
{
  int       res;
  char      buf[80];
  char     *env;

  env = getenv(envname);
  if (env && sscanf(env, "%d", &res) == 1)
    return res;

  (void) Rprintf("%s (%d)\n", question, deflt);
  res = deflt;
  if (fgets(buf, 80, stdin))
    sscanf(buf, "%d", &res);

  return res;
}

char
askforchar(const char *question, const char *envname, const char *allowed,
	   char deflt)
{
  char      res;
  char      buf[80];
  char     *env;
  const char *c;

  env = getenv(envname);
  if (env && sscanf(env, "%c", &res) == 1)
    return tolower(res);

  do {
    (void) Rprintf("%s (", question);
    res = tolower(deflt);
    for (c = allowed; *c; c++)
      if (*c == res)
	(void) Rprintf("%c", toupper(*c));
      else
	(void) Rprintf("%c", tolower(*c));
    (void) Rprintf(")\n");
    if (fgets(buf, 80, stdin))
      sscanf(buf, " %c", &res);
    res = tolower(res);

    for (c = allowed; *c && *c != res; c++);
  } while (!*c);

  return res;
}

real
askforreal(const char *question, const char *envname, real deflt)
{
  real      res;
  char      buf[80];
  char     *env;

  env = getenv(envname);
  if (env && sscanf(env, "%" SCANF_PREFIX "f", &res) == 1)
    return res;

  (void) Rprintf("%s (%.3e)\n", question, deflt);
  res = deflt;
  if (fgets(buf, 80, stdin))
    sscanf(buf, "%" SCANF_PREFIX "f", &res);

  return res;
}

char     *
askforstring(const char *question, const char *envname, const char *deflt,
	     char *buffer, uint bufsize)
{
  char     *env;
  char     *bufptr;

  assert(bufsize > 0);

  env = getenv(envname);
  if (env) {
    strncpy(buffer, env, bufsize);
    buffer[bufsize - 1] = '\0';
    return buffer;
  }

  (void) Rprintf("%s (%s)\n", question, deflt);
  if (fgets(buffer, bufsize, stdin) && buffer[0] != '\n' && buffer[0] != '\0') {
    for (bufptr = buffer; *bufptr != '\n' && *bufptr != '\0'; bufptr++);
    *bufptr = '\0';
    return buffer;
  }

  strncpy(buffer, deflt, bufsize);
  buffer[bufsize - 1] = '\0';
  return buffer;
}
